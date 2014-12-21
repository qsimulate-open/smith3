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
#include "forest.h"
#include "tree.h"
#include "residual.h"
#include "energy.h"
#include "dedci.h"
#include "density.h"
#include "density1.h"
#include "density2.h"
#include "correction.h"

using namespace std;
using namespace smith;

int main() {

  string theory="CAS_test";

  shared_ptr<Operator> ex_0 = make_shared<Ex>("x", "x", "a", "a");
  shared_ptr<Operator> t20 = make_shared<Op>("t2", "a", "a", "x", "x");
  shared_ptr<Operator> r0 = make_shared<Op>("r", "a", "a", "x", "x");
  shared_ptr<Operator> f1 = make_shared<Op>("f1", "g", "g");
  shared_ptr<Operator> v2 = make_shared<Op>("v2", "g", "g", "g", "g");
  shared_ptr<Operator> h1 = make_shared<Op>("h1", "g", "g");
  shared_ptr<Operator> proje = make_shared<Op>("proj");
  shared_ptr<Operator> t2dagger0 = make_shared<Op>("t2dagger", "x", "x", "a", "a");
  shared_ptr<Operator> ex_1b = make_shared<Ex>("g", "g");

  list<shared_ptr<Operator>> ra0 = {proje, ex_0, f1, t20};
  list<shared_ptr<Operator>> rb0 = {proje, ex_0, t20};
  list<shared_ptr<Operator>> rc0 = {proje, ex_0, v2};
  list<shared_ptr<Operator>> rd0 = {proje, ex_0, h1};
  auto dra0 = make_shared<Diagram>(ra0);
  auto drb0 = make_shared<Diagram>(rb0, -1, "e0");
  auto drc0 = make_shared<Diagram>(rc0);
  auto drd0 = make_shared<Diagram>(rd0);
  auto era0 = make_shared<Equation>(dra0, theory);
  auto erb0 = make_shared<Equation>(drb0, theory);
  auto erc0 = make_shared<Equation>(drc0, theory);
  auto erd0 = make_shared<Equation>(drd0, theory);
  era0->merge(erb0);
  era0->merge(erc0);
  era0->merge(erd0);
  era0->duplicates();
  era0->active();
  auto tra = make_shared<Residual>(era0, "residual");

  list<shared_ptr<Operator>> ea0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> eb0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> ec0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> ed0 = {proje, t2dagger0, h1};
  auto dea0 = make_shared<Diagram>(ea0, 0.25);
  auto deb0 = make_shared<Diagram>(eb0, -0.25, "e0");
  auto dec0 = make_shared<Diagram>(ec0, 0.5);
  auto ded0 = make_shared<Diagram>(ed0, 0.5);
  auto eea0 = make_shared<Equation>(dea0, theory);
  auto eeb0 = make_shared<Equation>(deb0, theory);
  auto eec0 = make_shared<Equation>(dec0, theory);
  auto eed0 = make_shared<Equation>(ded0, theory);
  eea0->merge(eeb0);
  eea0->merge(eec0);
  eea0->merge(eed0);
  eea0->duplicates();
  eea0->active();
  auto tea = make_shared<Energy>(eea0, "energy");

  list<shared_ptr<Operator>> ca0 = {proje, t2dagger0, t20};
  auto dca0 = make_shared<Diagram>(ca0, 0.25);
  auto eca0 = make_shared<Equation>(dca0, theory);
  eca0->duplicates();
  eca0->active();
  auto tca = make_shared<Correction>(eca0, "correction");

  list<shared_ptr<Operator>> da0 = {proje, t2dagger0, ex_1b, t20};
  auto dda0 = make_shared<Diagram>(da0, 0.25);
  auto eda0 = make_shared<Equation>(dda0, theory);
  eda0->duplicates();
  eda0->active();
  auto tda = make_shared<Density>(eda0, "density");

  list<shared_ptr<Operator>> db0 = {proje, ex_1b, t20};
  auto ddb0 = make_shared<Diagram>(db0);
  auto edb0 = make_shared<Equation>(ddb0, theory);
  edb0->duplicates();
  edb0->active();
  auto tdb = make_shared<Density1>(edb0, "density1");

  list<shared_ptr<Operator>> d2a0 = {proje, ex_0, t20};
  auto dd2a0 = make_shared<Diagram>(d2a0, 0.5);
  auto ed2a0 = make_shared<Equation>(dd2a0, theory);
  ed2a0->duplicates();
  ed2a0->active();
  auto td2a = make_shared<Density2>(ed2a0, "density2");

  list<shared_ptr<Operator>> dedcia0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcib0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcic0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcid0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcie0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcif0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcig0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> dedcih0 = {proje, t2dagger0, h1};
  auto ddedcia0 = make_shared<Diagram>(dedcia0, 0.25, make_pair(true, false));
  auto ddedcib0 = make_shared<Diagram>(dedcib0, 0.25, make_pair(false, true));
  auto ddedcic0 = make_shared<Diagram>(dedcic0, -0.25, "e0", make_pair(true, false));
  auto ddedcid0 = make_shared<Diagram>(dedcid0, -0.25, "e0", make_pair(false, true));
  auto ddedcie0 = make_shared<Diagram>(dedcie0, 0.5, make_pair(true, false));
  auto ddedcif0 = make_shared<Diagram>(dedcif0, 0.5, make_pair(false, true));
  auto ddedcig0 = make_shared<Diagram>(dedcig0, 0.5, make_pair(true, false));
  auto ddedcih0 = make_shared<Diagram>(dedcih0, 0.5, make_pair(false, true));
  auto ededcia0 = make_shared<Equation>(ddedcia0, theory);
  auto ededcib0 = make_shared<Equation>(ddedcib0, theory);
  auto ededcic0 = make_shared<Equation>(ddedcic0, theory);
  auto ededcid0 = make_shared<Equation>(ddedcid0, theory);
  auto ededcie0 = make_shared<Equation>(ddedcie0, theory);
  auto ededcif0 = make_shared<Equation>(ddedcif0, theory);
  auto ededcig0 = make_shared<Equation>(ddedcig0, theory);
  auto ededcih0 = make_shared<Equation>(ddedcih0, theory);
  ededcia0->merge(ededcib0);
  ededcia0->merge(ededcic0);
  ededcia0->merge(ededcid0);
  ededcia0->merge(ededcie0);
  ededcia0->merge(ededcif0);
  ededcia0->merge(ededcig0);
  ededcia0->merge(ededcih0);
  ededcia0->absorb_ket();
  ededcia0->duplicates();
  ededcia0->active();
  auto tdedcia = make_shared<Dedci>(ededcia0, "dedci");

  list<shared_ptr<Tree>> trees = {tra, tea, tca, tda, tdb, td2a, tdedcia};
  auto fr = make_shared<Forest>(trees);

  fr->filter_gamma();
  list<shared_ptr<Tensor>> gm = fr->gamma();
  const list<shared_ptr<Tensor>> gamma = gm;

  auto tmp = fr->generate_code();

  ofstream fs(fr->name() + ".h");
  ofstream es(fr->name() + "_tasks.h");
  ofstream cs(fr->name() + "_gen.cc");
  ofstream ds(fr->name() + "_tasks.cc");
  fs << tmp.ss.str();
  es << tmp.tt.str();
  cs << tmp.cc.str();
  ds << tmp.dd.str();
  fs.close();
  es.close();
  cs.close();
  ds.close();
  cout << std::endl;

  // output
  cout << std::endl << "   ***  Residual  ***" << std::endl << std::endl;
  tra->print();
  cout << std::endl << "   ***  Energy E2 ***" << std::endl << std::endl;
  tea->print();
  cout << std::endl << "   ***  Correlated norm <1|1> ***" << std::endl << std::endl;
  tca->print();
  cout << std::endl << "   ***  Correlated one-body density matrix d2 ***" << std::endl << std::endl;
  tda->print();
  cout << std::endl << "   ***  Correlated one-body density matrix d1 ***" << std::endl << std::endl;
  tdb->print();
  cout << std::endl << "   ***  Correlated two-body density matrix D1 ***" << std::endl << std::endl;
  td2a->print();
  cout << std::endl << "   ***  CI derivative  ***" << std::endl << std::endl;
  tdedcia->print();
  cout << std::endl << std::endl;

  return 0;
}

