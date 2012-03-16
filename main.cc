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
#include "equation.h"
#include "tree.h"

using namespace std;

int main() {

  shared_ptr<Op> proj(new Op("proj", "c", "c", "a", "a"));
  shared_ptr<Op> f(new Op("f1", "g", "g"));
  shared_ptr<Op> t(new Op("t2", "a", "a", "c", "c"));
  shared_ptr<Op> H(new Op("v2", "g", "g", "g", "g"));

  list<shared_ptr<Op> > d, e;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(t);
  e.push_back(proj);
  e.push_back(H);

  shared_ptr<Diagram> di(new Diagram(d));
  shared_ptr<Equation> eq(new Equation(di, "MP2"));
  shared_ptr<Diagram> dj(new Diagram(e));
  shared_ptr<Equation> eq2(new Equation(dj, "MP2"));
  eq->merge(eq2);
  eq->duplicates();
  eq->active();

  shared_ptr<Tree> res(new Tree(eq));
  res->print();

#if 0
  cout << "-----" << endl;
  cout << res->generate_task_list() << endl;
  cout << "-----" << endl;
#else
  ofstream fs(res->tree_name() + ".h");
  ofstream es(res->tree_name() + "_tasks.h");
  pair<string, string> tmp = res->generate_task_list();
  fs << tmp.first;
  es << tmp.second;
  fs.close();
  es.close();
#endif

  cout << res->generate() << endl;
}
