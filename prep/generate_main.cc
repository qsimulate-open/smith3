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

//////////////////////////////////////////////////////////////////////////////////
//
//  This program generates a simple CASPT2 theory generator
//  
//  Here the CASPT2 test includes two excitation operators, xxaa and xxxa 
//
//  compile: 
//  g++ -std=c++11 generate-main-caspt2-double.cc -o generate-main-caspt2-double
//  run:
//  ./generate-main-caspt2-double > main.cc 
//
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <tuple>
#include <string>
#include <memory>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <sstream>
#include <initializer_list>

#include "indices.h"




using namespace std;

pair<list<shared_ptr<Indices> >, list<shared_ptr<Indices> > > create_proj() {
  list<shared_ptr<Indices> > lp, lt;

  lp.push_back(shared_ptr<Indices>(new Indices{"x", "x", "a", "a"}));
  lp.push_back(shared_ptr<Indices>(new Indices{"x", "x", "x", "a"}));
  lp.push_back(shared_ptr<Indices>(new Indices{"c", "x", "x", "a"}));

  lt.push_back(shared_ptr<Indices>(new Indices{"a", "a", "x", "x"}));
  lt.push_back(shared_ptr<Indices>(new Indices{"a", "x", "x", "x"}));
  lt.push_back(shared_ptr<Indices>(new Indices{"a", "x", "c", "x"}));

  return make_pair(lp, lt);
};


// generators
string double_caspt2(list<shared_ptr<Indices> > dp, list<shared_ptr<Indices> > dt) {
  stringstream mm;

  assert(dp.size() == dt.size());
  int cnt = 0;
  for (auto pp = dp.begin(), tp = dt.begin(); pp != dp.end(); ++pp, ++tp, ++cnt) {
    stringstream ss; ss << cnt;
    mm << (*pp)->generate_tensor("proj" + ss.str());
    mm << (*tp)->generate_tensor("t" + ss.str());
    mm << (*tp)->generate_tensor("r" + ss.str());
  } 
  
  return mm.str();
} 


string generate_body() {
  stringstream mm;
  mm << "  shared_ptr<Op> dum(new Op(\"proj\"));" << endl;
  mm << "  shared_ptr<Op> f(new Op(\"f1\", \"g\", \"g\"));" << endl;
  mm << "  shared_ptr<Op> H(new Op(\"v2\", \"g\", \"g\", \"g\", \"g\"));" << endl;
  mm << "" << endl;
  mm << "  list<shared_ptr<Op> > d = {proj, f, t};" << endl;
  mm << "  list<shared_ptr<Op> > d2 = {proj2, f, t2};" << endl;
  mm << "  list<shared_ptr<Op> > d3 = {proj2, f, t};" << endl;
  mm << "  list<shared_ptr<Op> > d4 = {proj, f, t2};" << endl;
  mm << "  list<shared_ptr<Op> > e = {proj, H};" << endl;
  mm << "  list<shared_ptr<Op> > db = {proj, t};" << endl;
  mm << "" << endl;
  mm << "  // amplitude equation" << endl;
  mm << "  shared_ptr<Diagram> di(new Diagram(d));" << endl;
  mm << "  shared_ptr<Diagram> di2(new Diagram(d2));" << endl;
  mm << "  shared_ptr<Diagram> di3(new Diagram(d3));" << endl;
  mm << "  shared_ptr<Diagram> di4(new Diagram(d4));" << endl;
  mm << "  shared_ptr<Equation> eqn(new Equation(di, theory));" << endl;
  mm << "  shared_ptr<Equation> eqn2(new Equation(di2, theory));" << endl;
  mm << "  shared_ptr<Equation> eqn3(new Equation(di3, theory));" << endl;
  mm << "  shared_ptr<Equation> eqn4(new Equation(di4, theory));" << endl;
  mm << "  // TODO e0 is actually -e0, and therefore,  the screen output is currently wrong" << endl;
  mm << "  shared_ptr<Diagram> dib(new Diagram(db, \"e0\"));" << endl;
  mm << "" << endl;
  mm << "  shared_ptr<Equation> eqb(new Equation(dib, theory));" << endl;
  mm << "  shared_ptr<Diagram> dj(new Diagram(e));" << endl;
  mm << "  shared_ptr<Equation> eq2(new Equation(dj, theory));" << endl;
  mm << "  eqn->merge(eqn2);" << endl;
  mm << "  eqn->merge(eqn3);" << endl;
  mm << "  eqn->merge(eqn4);" << endl;
  mm << "  eqn->merge(eqb);" << endl;
  mm << "  eqn->merge(eq2);" << endl;
  mm << "  eqn->duplicates();" << endl;
  mm << "  eqn->active();" << endl;
  mm << "  shared_ptr<Tree> res(new Tree(eqn));" << endl;
  mm << "  res->sort_gamma();" << endl;
  return mm.str();
}


string header() {
  stringstream mm;
  mm << "//" << endl;
  mm << "// SMITH3 - generates spin-free multireference electron correlation programs." << endl;
  mm << "// Filename: main.cc" << endl;
  mm << "// Copyright (C) 2012 Toru Shiozaki" << endl;
  mm << "//" << endl;
  mm << "// Author: Toru Shiozaki <shiozaki@northwestern.edu>" << endl;
  mm << "// Maintainer: Shiozaki group" << endl;
  mm << "//" << endl;
  mm << "// This file is part of the SMITH3 package." << endl;
  mm << "//" << endl;
  mm << "// The SMITH3 package is free software; you can redistribute it and/or modify" << endl;
  mm << "// it under the terms of the GNU Library General Public License as published by" << endl;
  mm << "// the Free Software Foundation; either version 2, or (at your option)" << endl;
  mm << "// any later version." << endl;
  mm << "//" << endl;
  mm << "// The SMITH3 package is distributed in the hope that it will be useful," << endl;
  mm << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
  mm << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
  mm << "// GNU Library General Public License for more details." << endl;
  mm << "//" << endl;
  mm << "// You should have received a copy of the GNU Library General Public License" << endl;
  mm << "// along with the SMITH3 package; see COPYING.  If not, write to" << endl;
  mm << "// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA." << endl;
  mm << "//" << endl;
  mm << "" << endl;
  mm << "" << endl;
  mm << "// This program is supposed to perform Wick's theorem for multireference problems." << endl;
  mm << "// Spin averaged quantities assumed." << endl;
  mm << "" << endl;
  mm << "#include <iostream>" << endl;
  mm << "#include <fstream>" << endl;
  mm << "#include <list>" << endl;
  mm << "#include <string>" << endl;
  mm << "#include \"equation.h\"" << endl;
  mm << "#include \"tree.h\"" << endl;
  mm << "" << endl;
  mm << "using namespace std;" << endl;
  mm << "using namespace smith;" << endl;
  mm << "" << endl;
  mm << "int main() {" << endl;
  return mm.str();
};


string footer() {
  stringstream mm;
  mm << "  // energy" << endl;
  mm << "  list<shared_ptr<Op> > en0 = {dum, tdagger, H};" << endl;
  mm << "  list<shared_ptr<Op> > en1 = {dum, tdagger, R};" << endl;
  mm << "  shared_ptr<Diagram> e0(new Diagram(en0));" << endl;
  mm << "  shared_ptr<Equation> eneq(new Equation(e0, theory));" << endl;
  mm << "  shared_ptr<Diagram> e1(new Diagram(en1));" << endl;
  mm << "  shared_ptr<Equation> eneq1(new Equation(e1, theory));" << endl;
  mm << "  eneq->merge(eneq1);" << endl;
  mm << "  eneq->duplicates();" << endl;
  mm << "  eneq->active();" << endl;
  mm << "  shared_ptr<Tree> energy(new Tree(eneq));" << endl;
  mm << "  energy->sort_gamma(res->gamma());" << endl;
  mm << "" <<  endl; 
  mm << "  ofstream fs(res->tree_name() + \".h\");" << endl;
  mm << "  ofstream es(res->tree_name() + \"_tasks.h\");" << endl;
  mm << "  pair<string, string> tmp = res->generate_task_list(false, energy);" << endl;
  mm << "  fs << tmp.first;" << endl;
  mm << "  es << tmp.second;" << endl;
  mm << "  fs.close();" << endl;
  mm << "  es.close();" << endl;
  mm << "  cout << endl;" << endl;
  mm << "" <<  endl; 
  mm << "  // output" << endl;
  mm << "  cout << endl << \"   *** Residual ***\" << endl << endl;" << endl;
  mm << "  res->print();" << endl;
  mm << "  cout << endl << \"   ***  Energy  ***\" << endl << endl;" << endl;
  mm << "  energy->print();" << endl;
  mm << "  cout << endl << endl;" << endl;
  mm << "" <<  endl; 
  mm << "  return 0;" << endl;
  mm << "}" << endl;
  return mm.str();
};


int main() {

  // generate common header
  cout << header() << endl;

  list<shared_ptr<Indices> > proj_list, t_list;
  tie(proj_list, t_list) = create_proj();

  // generate caspt2 test
  cout << double_caspt2(proj_list, t_list) << endl;

  cout << generate_body() << endl;

  // done generate the footer
  cout << footer() << endl;


  return 0;
}


