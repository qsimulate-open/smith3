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
// The SMITH3 package is free software; you can redistribute it and\/or modify
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

// MP2
#if 0
  shared_ptr<Op> proj(new Op("proj", "c", "c", "a", "a")); 
  shared_ptr<Op> t(new Op("t2", "a", "a", "c", "c"));
  shared_ptr<Op> R(new Op("r", "a", "a", "c", "c"));
  shared_ptr<Op> tdagger(new Op("t2dagger", "c", "c", "a", "a"));
  string theory="MP2";
#endif
// simple CASPT2
#if 1
  shared_ptr<Op> proj(new Op("proj", "x", "x", "a", "a")); // test active 1
  shared_ptr<Op> t(new Op("t2", "a", "a", "x", "x"));  // test active 1
  shared_ptr<Op> R(new Op("r", "a", "a", "x", "x"));
  shared_ptr<Op> tdagger(new Op("t2dagger", "x", "x", "a", "a"));
  string theory="CAS_all_active";
#endif
// complicated one 
#if 0
  shared_ptr<Op> tdagger(new Op("t2dagger", "x", "x", "a", "a"));
  shared_ptr<Op> R(new Op("r", "a", "a", "x", "x"));
  shared_ptr<Op> proj(new Op("proj", "c", "x", "x", "x"));
  shared_ptr<Op> t(new Op("t2", "x", "x", "x", "c"));
  string theory="CAS";
#endif

  shared_ptr<Op> dum(new Op("proj"));
  shared_ptr<Op> f(new Op("f1", "g", "g"));
  shared_ptr<Op> H(new Op("v2", "g", "g", "g", "g"));

  list<shared_ptr<Op> > d = {proj, f, t};
  list<shared_ptr<Op> > e = {proj, H};
  list<shared_ptr<Op> > db = {proj, t};

  // amplitude equation
  shared_ptr<Diagram> di(new Diagram(d));
  shared_ptr<Equation> eq(new Equation(di, theory));
  // TODO e0 is actually -e0, and therefore,  the screen output is currently wrong
  shared_ptr<Diagram> dib(new Diagram(db, "e0"));
 
  shared_ptr<Equation> eqb(new Equation(dib, theory));
  shared_ptr<Diagram> dj(new Diagram(e));
  shared_ptr<Equation> eq2(new Equation(dj, theory));
  eq->merge(eqb);
  eq->merge(eq2);
  eq->duplicates();
  eq->active();
  shared_ptr<Tree> res(new Tree(eq));
  res->sort_gamma();

  // energy
  list<shared_ptr<Op> > en0 = {dum, tdagger, H};
  list<shared_ptr<Op> > en1 = {dum, tdagger, R};
  shared_ptr<Diagram> e0(new Diagram(en0));
  shared_ptr<Equation> eneq(new Equation(e0, theory));
  shared_ptr<Diagram> e1(new Diagram(en1));
  shared_ptr<Equation> eneq1(new Equation(e1, theory));
  eneq->merge(eneq1);
  eneq->duplicates();
  eneq->active();
  shared_ptr<Tree> energy(new Tree(eneq));
  energy->sort_gamma(res->gamma());

  ofstream fs(res->tree_name() + ".h");
  ofstream es(res->tree_name() + "_tasks.h");
  pair<string, string> tmp = res->generate_task_list(false,energy);
  fs << tmp.first;
  es << tmp.second;
  fs.close();
  es.close();
  cout << endl;

  // output
  cout << endl << "   *** Residual ***" << endl << endl;
  res->print();
  cout << endl << "   ***  Energy  ***" << endl << endl;
  energy->print();
  cout << endl << endl;

  return 0;
}
