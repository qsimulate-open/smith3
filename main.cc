//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: main.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the SMITH3 package.
//
// The SMITH3 package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SMITH3 package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SMITH3 package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// This program is supposed to perform Wick's theorem for multireference problems.
// Spin averaged quantities assumed.

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include "equation.h"
#include "tree.h"
#include "residual.h"
#include "energy.h"
#include "density.h"

using namespace std;
using namespace smith;

int main() {

  string theory="CAS_test"; 
  shared_ptr<Operator> proj0(new Op("proj", "x", "c", "x", "x"));
  shared_ptr<Operator> t20(new Op("t2", "x", "x", "x", "c"));
  shared_ptr<Operator> r0(new Op("r", "x", "x", "x", "c"));
  shared_ptr<Operator> f1(new Op("f1", "g", "g"));
  shared_ptr<Operator> v2(new Op("v2", "g", "g", "g", "g"));
  shared_ptr<Operator> h1(new Op("h1", "g", "g"));
  shared_ptr<Operator> proje(new Op("proj"));
  shared_ptr<Operator> t2dagger0(new Op("t2dagger", "x", "c", "x", "x"));
  list<shared_ptr<Operator>> da0 = {proj0, f1, t20};
  list<shared_ptr<Operator>> db0 = {proj0, t20};
  list<shared_ptr<Operator>> dc0 = {proj0, v2};
  list<shared_ptr<Operator>> dd0 = {proj0, h1};
  shared_ptr<Diagram> dda0(new Diagram(da0));
  shared_ptr<Diagram> ddb0(new Diagram(db0, "e0"));
  shared_ptr<Diagram> ddc0(new Diagram(dc0));
  shared_ptr<Diagram> ddd0(new Diagram(dd0));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  shared_ptr<Equation> edc0(new Equation(ddc0, theory));
  shared_ptr<Equation> edd0(new Equation(ddd0, theory));
  eda0->merge(edb0);
  eda0->merge(edc0);
  eda0->merge(edd0);
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Residual(eda0, "residual"));
  tda->sort_gamma();
  list<shared_ptr<Operator>> ea0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> eb0 = {proje, t2dagger0, h1};
  shared_ptr<Diagram> dea0(new Diagram(ea0));
  shared_ptr<Diagram> deb0(new Diagram(eb0));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  eea0->merge(eeb0);
  eea0->duplicates();
  eea0->active();
  shared_ptr<Tree> tea(new Energy(eea0, "energy"));
  tea->sort_gamma(tda->gamma());

  ofstream fs(tda->tree_name() + ".h");
  ofstream es(tda->tree_name() + "_tasks.h");
  list<shared_ptr<Tree>> tea_list = {tea};
  pair<string, string> tmp = tda->generate_task_list(tea_list);
  fs << tmp.first;
  es << tmp.second;
  fs.close();
  es.close();
  cout << std::endl;

  // output
  cout << std::endl << "   *** Residual ***" << std::endl << std::endl;
  tda->print();
  cout << std::endl << "   ***  Energy  ***" << std::endl << std::endl;
  tea->print();
  cout << std::endl << std::endl;

  return 0;
}

