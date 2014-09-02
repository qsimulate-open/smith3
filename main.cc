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

  shared_ptr<Operator> ex_0(new Ex("c", "c", "a", "a"));
  shared_ptr<Operator> t20(new Op("t2", "a", "a", "c", "c"));
  shared_ptr<Operator> r0(new Op("r", "a", "a", "c", "c"));
  shared_ptr<Operator> f1(new Op("f1", "g", "g"));
  shared_ptr<Operator> v2(new Op("v2", "g", "g", "g", "g"));
  shared_ptr<Operator> h1(new Op("h1", "g", "g"));
  shared_ptr<Operator> proje(new Op("proj"));
  shared_ptr<Operator> t2dagger0(new Op("t2dagger", "c", "c", "a", "a"));
  shared_ptr<Operator> ex_1b(new Ex("g", "g"));

  list<shared_ptr<Operator>> ra0 = {proje, ex_0, f1, t20};
  list<shared_ptr<Operator>> rb0 = {proje, ex_0, t20};
  list<shared_ptr<Operator>> rc0 = {proje, ex_0, v2};
  list<shared_ptr<Operator>> rd0 = {proje, ex_0, h1};
  shared_ptr<Diagram> dra0(new Diagram(ra0));
  shared_ptr<Diagram> drb0(new Diagram(rb0, -1, "e0"));
  shared_ptr<Diagram> drc0(new Diagram(rc0));
  shared_ptr<Diagram> drd0(new Diagram(rd0));
  shared_ptr<Equation> era0(new Equation(dra0, theory));
  shared_ptr<Equation> erb0(new Equation(drb0, theory));
  shared_ptr<Equation> erc0(new Equation(drc0, theory));
  shared_ptr<Equation> erd0(new Equation(drd0, theory));
  era0->merge(erb0);
  era0->merge(erc0);
  era0->merge(erd0);
  era0->duplicates();
  era0->active();
  shared_ptr<Tree> tra(new Residual(era0, "residual"));

  list<shared_ptr<Operator>> ea0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> eb0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> ec0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> ed0 = {proje, t2dagger0, h1};
  shared_ptr<Diagram> dea0(new Diagram(ea0, 0.25));
  shared_ptr<Diagram> deb0(new Diagram(eb0, -0.25, "e0"));
  shared_ptr<Diagram> dec0(new Diagram(ec0, 0.5));
  shared_ptr<Diagram> ded0(new Diagram(ed0, 0.5));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  shared_ptr<Equation> eec0(new Equation(dec0, theory));
  shared_ptr<Equation> eed0(new Equation(ded0, theory));
  eea0->merge(eeb0);
  eea0->merge(eec0);
  eea0->merge(eed0);
  eea0->duplicates();
  eea0->active();
  shared_ptr<Tree> tea(new Energy(eea0, "energy"));

  list<shared_ptr<Operator>> ca0 = {proje, t2dagger0, t20};
  shared_ptr<Diagram> dca0(new Diagram(ca0, 0.25));
  shared_ptr<Equation> eca0(new Equation(dca0, theory));
  eca0->duplicates();
  eca0->active();
  shared_ptr<Tree> tca(new Correction(eca0, "correction"));

  list<shared_ptr<Operator>> da0 = {proje, t2dagger0, ex_1b, t20};
  shared_ptr<Diagram> dda0(new Diagram(da0, 0.25));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Density(eda0, "density"));

  list<shared_ptr<Operator>> db0 = {proje, ex_1b, t20};
  shared_ptr<Diagram> ddb0(new Diagram(db0));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  edb0->duplicates();
  edb0->active();
  shared_ptr<Tree> tdb(new Density1(edb0, "density1"));

  list<shared_ptr<Operator>> d2a0 = {proje, ex_0, t20};
  shared_ptr<Diagram> dd2a0(new Diagram(d2a0, 0.5));
  shared_ptr<Equation> ed2a0(new Equation(dd2a0, theory));
  ed2a0->duplicates();
  ed2a0->active();
  shared_ptr<Tree> td2a(new Density2(ed2a0, "density2"));

  list<shared_ptr<Operator>> dedcia0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcib0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcic0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcid0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcie0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcif0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcig0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> dedcih0 = {proje, t2dagger0, h1};
  shared_ptr<Diagram> ddedcia0(new Diagram(dedcia0, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcib0(new Diagram(dedcib0, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcic0(new Diagram(dedcic0, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcid0(new Diagram(dedcid0, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcie0(new Diagram(dedcie0, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcif0(new Diagram(dedcif0, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcig0(new Diagram(dedcig0, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcih0(new Diagram(dedcih0, 0.5, std::make_pair(false, true)));
  shared_ptr<Equation> ededcia0(new Equation(ddedcia0, theory));
  shared_ptr<Equation> ededcib0(new Equation(ddedcib0, theory));
  shared_ptr<Equation> ededcic0(new Equation(ddedcic0, theory));
  shared_ptr<Equation> ededcid0(new Equation(ddedcid0, theory));
  shared_ptr<Equation> ededcie0(new Equation(ddedcie0, theory));
  shared_ptr<Equation> ededcif0(new Equation(ddedcif0, theory));
  shared_ptr<Equation> ededcig0(new Equation(ddedcig0, theory));
  shared_ptr<Equation> ededcih0(new Equation(ddedcih0, theory));
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
  shared_ptr<Tree> tdedcia(new Dedci(ededcia0, "dedci"));

  list<shared_ptr<Tree>> trees = {tra, tea, tca, tda, tdb, td2a, tdedcia};
  shared_ptr<Forest> fr(new Forest(trees));

  fr->filter_gamma();
  list<shared_ptr<Tensor>> gm = fr->gamma();
  const list<shared_ptr<Tensor>> gamma = gm;

  pair<string, string> tmp = fr->generate_code();

  ofstream fs(fr->name() + ".h");
  ofstream es(fr->name() + "_tasks.h");
  fs << tmp.first;
  es << tmp.second;
  fs.close();
  es.close();
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

