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

using namespace std;
using namespace smith;

int main() {

  string theory="MP2"; 
  shared_ptr<Op> proj0(new Op("proj", "c", "c", "a", "a"));
  shared_ptr<Op> t0(new Op("t", "a", "a", "c", "c"));
  shared_ptr<Op> r0(new Op("r", "a", "a", "c", "c"));
  shared_ptr<Op> f(new Op("f", "g", "g"));
  shared_ptr<Op> H(new Op("H", "g", "g", "g", "g"));
  shared_ptr<Op> proje(new Op("proj"));
  shared_ptr<Op> tdagger0(new Op("tdagger", "c", "c", "a", "a"));
  list<shared_ptr<Op> > da0 = {proj0, f, t0};
  list<shared_ptr<Op> > db0 = {proj0, t0};
  list<shared_ptr<Op> > dc0 = {proj0, H};
  shared_ptr<Diagram> dda0(new Diagram(da0));
  shared_ptr<Diagram> ddb0(new Diagram(db0));
  shared_ptr<Diagram> ddc0(new Diagram(dc0));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  shared_ptr<Equation> edc0(new Equation(ddc0, theory));
  eda0->merge(edb0);
  eda0->merge(edc0);
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Tree(eda0));
  list<shared_ptr<Op> > ea0 = {proje, tdagger0, H};
  list<shared_ptr<Op> > eb0 = {proje, tdagger0, r0};
  shared_ptr<Diagram> dea0(new Diagram(ea0));
  shared_ptr<Diagram> deb0(new Diagram(eb0));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  eea0->merge(eeb0);
  eea0->duplicates();
  eea0->active();
  shared_ptr<Tree> tea(new Tree(eea0));
  tea->sort_gamma(tda->gamma());

  ofstream fs(tda->tree_name() + ".h");
  ofstream es(tda->tree_name() + "_tasks.h");
  pair<string, std::string> tmp = tda->generate_task_list(false, tea);
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

