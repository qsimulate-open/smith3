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

  string theory="CAS_all_active"; 
  shared_ptr<Op> proj0(new Op("proj", "x", "x", "x", "a"));
  shared_ptr<Op> proj1(new Op("proj", "x", "x", "a", "a"));
  shared_ptr<Op> t20(new Op("t2", "x", "a", "x", "x"));
  shared_ptr<Op> t21(new Op("t2", "a", "a", "x", "x"));
  shared_ptr<Op> r0(new Op("r", "x", "a", "x", "x"));
  shared_ptr<Op> r1(new Op("r", "a", "a", "x", "x"));
  shared_ptr<Op> f1(new Op("f1", "g", "g"));
  shared_ptr<Op> v2(new Op("v2", "g", "g", "g", "g"));
  shared_ptr<Op> proje(new Op("proj"));
  shared_ptr<Op> t2dagger0(new Op("t2dagger", "x", "x", "x", "a"));
  shared_ptr<Op> t2dagger1(new Op("t2dagger", "x", "x", "a", "a"));
  list<shared_ptr<Op> > da0 = {proj0, f1, t20};
  list<shared_ptr<Op> > da1 = {proj0, f1, t21};
  list<shared_ptr<Op> > da2 = {proj1, f1, t20};
  list<shared_ptr<Op> > da3 = {proj1, f1, t21};
  list<shared_ptr<Op> > db0 = {proj0, t20};
  list<shared_ptr<Op> > db1 = {proj0, t21};
  list<shared_ptr<Op> > db2 = {proj1, t20};
  list<shared_ptr<Op> > db3 = {proj1, t21};
  list<shared_ptr<Op> > dc0 = {proj0, v2};
  list<shared_ptr<Op> > dc1 = {proj1, v2};
  shared_ptr<Diagram> dda0(new Diagram(da0));
  shared_ptr<Diagram> dda1(new Diagram(da1));
  shared_ptr<Diagram> dda2(new Diagram(da2));
  shared_ptr<Diagram> dda3(new Diagram(da3));
  shared_ptr<Diagram> ddb0(new Diagram(db0, "e0"));
  shared_ptr<Diagram> ddb1(new Diagram(db1, "e0"));
  shared_ptr<Diagram> ddb2(new Diagram(db2, "e0"));
  shared_ptr<Diagram> ddb3(new Diagram(db3, "e0"));
  shared_ptr<Diagram> ddc0(new Diagram(dc0));
  shared_ptr<Diagram> ddc1(new Diagram(dc1));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  shared_ptr<Equation> eda1(new Equation(dda1, theory));
  shared_ptr<Equation> eda2(new Equation(dda2, theory));
  shared_ptr<Equation> eda3(new Equation(dda3, theory));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  shared_ptr<Equation> edb1(new Equation(ddb1, theory));
  shared_ptr<Equation> edb2(new Equation(ddb2, theory));
  shared_ptr<Equation> edb3(new Equation(ddb3, theory));
  shared_ptr<Equation> edc0(new Equation(ddc0, theory));
  shared_ptr<Equation> edc1(new Equation(ddc1, theory));
  eda0->merge(eda1);
  eda0->merge(eda2);
  eda0->merge(eda3);
  eda0->merge(edb0);
  eda0->merge(edb1);
  eda0->merge(edb2);
  eda0->merge(edb3);
  eda0->merge(edc0);
  eda0->merge(edc1);
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Tree(eda0));
  tda->sort_gamma();
  list<shared_ptr<Op> > ea0 = {proje, t2dagger0, v2};
  list<shared_ptr<Op> > ea1 = {proje, t2dagger1, v2};
  list<shared_ptr<Op> > eb0 = {proje, t2dagger0, r0};
  list<shared_ptr<Op> > eb1 = {proje, t2dagger0, r1};
  list<shared_ptr<Op> > eb2 = {proje, t2dagger1, r0};
  list<shared_ptr<Op> > eb3 = {proje, t2dagger1, r1};
  shared_ptr<Diagram> dea0(new Diagram(ea0));
  shared_ptr<Diagram> dea1(new Diagram(ea1));
  shared_ptr<Diagram> deb0(new Diagram(eb0));
  shared_ptr<Diagram> deb1(new Diagram(eb1));
  shared_ptr<Diagram> deb2(new Diagram(eb2));
  shared_ptr<Diagram> deb3(new Diagram(eb3));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eea1(new Equation(dea1, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  shared_ptr<Equation> eeb1(new Equation(deb1, theory));
  shared_ptr<Equation> eeb2(new Equation(deb2, theory));
  shared_ptr<Equation> eeb3(new Equation(deb3, theory));
  eea0->merge(eea1);
  eea0->merge(eeb0);
  eea0->merge(eeb1);
  eea0->merge(eeb2);
  eea0->merge(eeb3);
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

