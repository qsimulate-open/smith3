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
#include "density.h"
#include "density2.h"
#include "correction.h"

using namespace std;
using namespace smith;

int main() {

  string theory="CAS_test";

  shared_ptr<Operator> ex_0(new Ex("c", "c", "a", "a"));
  shared_ptr<Operator> ex_1(new Ex("x", "c", "a", "a"));
  shared_ptr<Operator> ex_2(new Ex("x", "x", "a", "a"));
  shared_ptr<Operator> t20(new Op("t2", "a", "a", "c", "c"));
  shared_ptr<Operator> t21(new Op("t2", "a", "a", "x", "c"));
  shared_ptr<Operator> t22(new Op("t2", "a", "a", "x", "x"));
  shared_ptr<Operator> r0(new Op("r", "a", "a", "c", "c"));
  shared_ptr<Operator> r1(new Op("r", "a", "a", "x", "c"));
  shared_ptr<Operator> r2(new Op("r", "a", "a", "x", "x"));
  shared_ptr<Operator> f1(new Op("f1", "g", "g"));
  shared_ptr<Operator> v2(new Op("v2", "g", "g", "g", "g"));
  shared_ptr<Operator> h1(new Op("h1", "g", "g"));
  shared_ptr<Operator> proje(new Op("proj"));
  shared_ptr<Operator> t2dagger0(new Op("t2dagger", "c", "c", "a", "a"));
  shared_ptr<Operator> t2dagger1(new Op("t2dagger", "x", "c", "a", "a"));
  shared_ptr<Operator> t2dagger2(new Op("t2dagger", "x", "x", "a", "a"));
  shared_ptr<Operator> ex_1b(new Ex("g", "g"));

  list<shared_ptr<Operator>> ra0 = {proje, ex_0, f1, t20};
  list<shared_ptr<Operator>> ra1 = {proje, ex_0, f1, t21};
  list<shared_ptr<Operator>> ra2 = {proje, ex_0, f1, t22};
  list<shared_ptr<Operator>> ra3 = {proje, ex_1, f1, t20};
  list<shared_ptr<Operator>> ra4 = {proje, ex_1, f1, t21};
  list<shared_ptr<Operator>> ra5 = {proje, ex_1, f1, t22};
  list<shared_ptr<Operator>> ra6 = {proje, ex_2, f1, t20};
  list<shared_ptr<Operator>> ra7 = {proje, ex_2, f1, t21};
  list<shared_ptr<Operator>> ra8 = {proje, ex_2, f1, t22};
  list<shared_ptr<Operator>> rb0 = {proje, ex_0, t20};
  list<shared_ptr<Operator>> rb1 = {proje, ex_0, t21};
  list<shared_ptr<Operator>> rb2 = {proje, ex_0, t22};
  list<shared_ptr<Operator>> rb3 = {proje, ex_1, t20};
  list<shared_ptr<Operator>> rb4 = {proje, ex_1, t21};
  list<shared_ptr<Operator>> rb5 = {proje, ex_1, t22};
  list<shared_ptr<Operator>> rb6 = {proje, ex_2, t20};
  list<shared_ptr<Operator>> rb7 = {proje, ex_2, t21};
  list<shared_ptr<Operator>> rb8 = {proje, ex_2, t22};
  list<shared_ptr<Operator>> rc0 = {proje, ex_0, v2};
  list<shared_ptr<Operator>> rc1 = {proje, ex_1, v2};
  list<shared_ptr<Operator>> rc2 = {proje, ex_2, v2};
  list<shared_ptr<Operator>> rd0 = {proje, ex_0, h1};
  list<shared_ptr<Operator>> rd1 = {proje, ex_1, h1};
  list<shared_ptr<Operator>> rd2 = {proje, ex_2, h1};
  shared_ptr<Diagram> dra0(new Diagram(ra0));
  shared_ptr<Diagram> dra1(new Diagram(ra1));
  shared_ptr<Diagram> dra2(new Diagram(ra2));
  shared_ptr<Diagram> dra3(new Diagram(ra3));
  shared_ptr<Diagram> dra4(new Diagram(ra4));
  shared_ptr<Diagram> dra5(new Diagram(ra5));
  shared_ptr<Diagram> dra6(new Diagram(ra6));
  shared_ptr<Diagram> dra7(new Diagram(ra7));
  shared_ptr<Diagram> dra8(new Diagram(ra8));
  shared_ptr<Diagram> drb0(new Diagram(rb0, "e0"));
  shared_ptr<Diagram> drb1(new Diagram(rb1, "e0"));
  shared_ptr<Diagram> drb2(new Diagram(rb2, "e0"));
  shared_ptr<Diagram> drb3(new Diagram(rb3, "e0"));
  shared_ptr<Diagram> drb4(new Diagram(rb4, "e0"));
  shared_ptr<Diagram> drb5(new Diagram(rb5, "e0"));
  shared_ptr<Diagram> drb6(new Diagram(rb6, "e0"));
  shared_ptr<Diagram> drb7(new Diagram(rb7, "e0"));
  shared_ptr<Diagram> drb8(new Diagram(rb8, "e0"));
  shared_ptr<Diagram> drc0(new Diagram(rc0));
  shared_ptr<Diagram> drc1(new Diagram(rc1));
  shared_ptr<Diagram> drc2(new Diagram(rc2));
  shared_ptr<Diagram> drd0(new Diagram(rd0));
  shared_ptr<Diagram> drd1(new Diagram(rd1));
  shared_ptr<Diagram> drd2(new Diagram(rd2));
  shared_ptr<Equation> era0(new Equation(dra0, theory));
  shared_ptr<Equation> era1(new Equation(dra1, theory));
  shared_ptr<Equation> era2(new Equation(dra2, theory));
  shared_ptr<Equation> era3(new Equation(dra3, theory));
  shared_ptr<Equation> era4(new Equation(dra4, theory));
  shared_ptr<Equation> era5(new Equation(dra5, theory));
  shared_ptr<Equation> era6(new Equation(dra6, theory));
  shared_ptr<Equation> era7(new Equation(dra7, theory));
  shared_ptr<Equation> era8(new Equation(dra8, theory));
  shared_ptr<Equation> erb0(new Equation(drb0, theory));
  shared_ptr<Equation> erb1(new Equation(drb1, theory));
  shared_ptr<Equation> erb2(new Equation(drb2, theory));
  shared_ptr<Equation> erb3(new Equation(drb3, theory));
  shared_ptr<Equation> erb4(new Equation(drb4, theory));
  shared_ptr<Equation> erb5(new Equation(drb5, theory));
  shared_ptr<Equation> erb6(new Equation(drb6, theory));
  shared_ptr<Equation> erb7(new Equation(drb7, theory));
  shared_ptr<Equation> erb8(new Equation(drb8, theory));
  shared_ptr<Equation> erc0(new Equation(drc0, theory));
  shared_ptr<Equation> erc1(new Equation(drc1, theory));
  shared_ptr<Equation> erc2(new Equation(drc2, theory));
  shared_ptr<Equation> erd0(new Equation(drd0, theory));
  shared_ptr<Equation> erd1(new Equation(drd1, theory));
  shared_ptr<Equation> erd2(new Equation(drd2, theory));
  era0->merge(era1);
  era0->merge(era2);
  era0->merge(era3);
  era0->merge(era4);
  era0->merge(era5);
  era0->merge(era6);
  era0->merge(era7);
  era0->merge(era8);
  era0->merge(erb0);
  era0->merge(erb1);
  era0->merge(erb2);
  era0->merge(erb3);
  era0->merge(erb4);
  era0->merge(erb5);
  era0->merge(erb6);
  era0->merge(erb7);
  era0->merge(erb8);
  era0->merge(erc0);
  era0->merge(erc1);
  era0->merge(erc2);
  era0->merge(erd0);
  era0->merge(erd1);
  era0->merge(erd2);
  era0->duplicates();
  era0->active();
  shared_ptr<Tree> tra(new Residual(era0, "residual"));

  list<shared_ptr<Operator>> ea0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> ea1 = {proje, t2dagger1, v2};
  list<shared_ptr<Operator>> ea2 = {proje, t2dagger2, v2};
  list<shared_ptr<Operator>> eb0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> eb1 = {proje, t2dagger1, h1};
  list<shared_ptr<Operator>> eb2 = {proje, t2dagger2, h1};
  shared_ptr<Diagram> dea0(new Diagram(ea0));
  shared_ptr<Diagram> dea1(new Diagram(ea1));
  shared_ptr<Diagram> dea2(new Diagram(ea2));
  shared_ptr<Diagram> deb0(new Diagram(eb0));
  shared_ptr<Diagram> deb1(new Diagram(eb1));
  shared_ptr<Diagram> deb2(new Diagram(eb2));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eea1(new Equation(dea1, theory));
  shared_ptr<Equation> eea2(new Equation(dea2, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  shared_ptr<Equation> eeb1(new Equation(deb1, theory));
  shared_ptr<Equation> eeb2(new Equation(deb2, theory));
  eea0->merge(eea1);
  eea0->merge(eea2);
  eea0->merge(eeb0);
  eea0->merge(eeb1);
  eea0->merge(eeb2);
  eea0->duplicates();
  eea0->active();
  shared_ptr<Tree> tea(new Energy(eea0, "energy"));

  list<shared_ptr<Operator>> ca0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> ca1 = {proje, t2dagger0, t21};
  list<shared_ptr<Operator>> ca2 = {proje, t2dagger0, t22};
  list<shared_ptr<Operator>> ca3 = {proje, t2dagger1, t20};
  list<shared_ptr<Operator>> ca4 = {proje, t2dagger1, t21};
  list<shared_ptr<Operator>> ca5 = {proje, t2dagger1, t22};
  list<shared_ptr<Operator>> ca6 = {proje, t2dagger2, t20};
  list<shared_ptr<Operator>> ca7 = {proje, t2dagger2, t21};
  list<shared_ptr<Operator>> ca8 = {proje, t2dagger2, t22};
  shared_ptr<Diagram> dca0(new Diagram(ca0, 0.25));
  shared_ptr<Diagram> dca1(new Diagram(ca1, 0.25));
  shared_ptr<Diagram> dca2(new Diagram(ca2, 0.25));
  shared_ptr<Diagram> dca3(new Diagram(ca3, 0.25));
  shared_ptr<Diagram> dca4(new Diagram(ca4, 0.25));
  shared_ptr<Diagram> dca5(new Diagram(ca5, 0.25));
  shared_ptr<Diagram> dca6(new Diagram(ca6, 0.25));
  shared_ptr<Diagram> dca7(new Diagram(ca7, 0.25));
  shared_ptr<Diagram> dca8(new Diagram(ca8, 0.25));
  shared_ptr<Equation> eca0(new Equation(dca0, theory));
  shared_ptr<Equation> eca1(new Equation(dca1, theory));
  shared_ptr<Equation> eca2(new Equation(dca2, theory));
  shared_ptr<Equation> eca3(new Equation(dca3, theory));
  shared_ptr<Equation> eca4(new Equation(dca4, theory));
  shared_ptr<Equation> eca5(new Equation(dca5, theory));
  shared_ptr<Equation> eca6(new Equation(dca6, theory));
  shared_ptr<Equation> eca7(new Equation(dca7, theory));
  shared_ptr<Equation> eca8(new Equation(dca8, theory));
  eca0->merge(eca1);
  eca0->merge(eca2);
  eca0->merge(eca3);
  eca0->merge(eca4);
  eca0->merge(eca5);
  eca0->merge(eca6);
  eca0->merge(eca7);
  eca0->merge(eca8);
  eca0->duplicates();
  eca0->active();
  shared_ptr<Tree> tca(new Correction(eca0, "correction"));

  list<shared_ptr<Operator>> da0 = {proje, t2dagger0, ex_1b, t20};
  list<shared_ptr<Operator>> da1 = {proje, t2dagger0, ex_1b, t21};
  list<shared_ptr<Operator>> da2 = {proje, t2dagger0, ex_1b, t22};
  list<shared_ptr<Operator>> da3 = {proje, t2dagger1, ex_1b, t20};
  list<shared_ptr<Operator>> da4 = {proje, t2dagger1, ex_1b, t21};
  list<shared_ptr<Operator>> da5 = {proje, t2dagger1, ex_1b, t22};
  list<shared_ptr<Operator>> da6 = {proje, t2dagger2, ex_1b, t20};
  list<shared_ptr<Operator>> da7 = {proje, t2dagger2, ex_1b, t21};
  list<shared_ptr<Operator>> da8 = {proje, t2dagger2, ex_1b, t22};
  list<shared_ptr<Operator>> db0 = {proje, ex_1b, t20};
  list<shared_ptr<Operator>> db1 = {proje, ex_1b, t21};
  list<shared_ptr<Operator>> db2 = {proje, ex_1b, t22};
  shared_ptr<Diagram> dda0(new Diagram(da0, 0.25));
  shared_ptr<Diagram> dda1(new Diagram(da1, 0.25));
  shared_ptr<Diagram> dda2(new Diagram(da2, 0.25));
  shared_ptr<Diagram> dda3(new Diagram(da3, 0.25));
  shared_ptr<Diagram> dda4(new Diagram(da4, 0.25));
  shared_ptr<Diagram> dda5(new Diagram(da5, 0.25));
  shared_ptr<Diagram> dda6(new Diagram(da6, 0.25));
  shared_ptr<Diagram> dda7(new Diagram(da7, 0.25));
  shared_ptr<Diagram> dda8(new Diagram(da8, 0.25));
  shared_ptr<Diagram> ddb0(new Diagram(db0));
  shared_ptr<Diagram> ddb1(new Diagram(db1));
  shared_ptr<Diagram> ddb2(new Diagram(db2));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  shared_ptr<Equation> eda1(new Equation(dda1, theory));
  shared_ptr<Equation> eda2(new Equation(dda2, theory));
  shared_ptr<Equation> eda3(new Equation(dda3, theory));
  shared_ptr<Equation> eda4(new Equation(dda4, theory));
  shared_ptr<Equation> eda5(new Equation(dda5, theory));
  shared_ptr<Equation> eda6(new Equation(dda6, theory));
  shared_ptr<Equation> eda7(new Equation(dda7, theory));
  shared_ptr<Equation> eda8(new Equation(dda8, theory));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  shared_ptr<Equation> edb1(new Equation(ddb1, theory));
  shared_ptr<Equation> edb2(new Equation(ddb2, theory));
  eda0->merge(eda1);
  eda0->merge(eda2);
  eda0->merge(eda3);
  eda0->merge(eda4);
  eda0->merge(eda5);
  eda0->merge(eda6);
  eda0->merge(eda7);
  eda0->merge(eda8);
  eda0->merge(edb0);
  eda0->merge(edb1);
  eda0->merge(edb2);
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Density(eda0, "density"));

  list<shared_ptr<Operator>> d2a0 = {proje, ex_0, t20};
  list<shared_ptr<Operator>> d2a1 = {proje, ex_0, t21};
  list<shared_ptr<Operator>> d2a2 = {proje, ex_0, t22};
  list<shared_ptr<Operator>> d2a3 = {proje, ex_1, t20};
  list<shared_ptr<Operator>> d2a4 = {proje, ex_1, t21};
  list<shared_ptr<Operator>> d2a5 = {proje, ex_1, t22};
  list<shared_ptr<Operator>> d2a6 = {proje, ex_2, t20};
  list<shared_ptr<Operator>> d2a7 = {proje, ex_2, t21};
  list<shared_ptr<Operator>> d2a8 = {proje, ex_2, t22};
  shared_ptr<Diagram> dd2a0(new Diagram(d2a0, 0.5));
  shared_ptr<Diagram> dd2a1(new Diagram(d2a1, 0.5));
  shared_ptr<Diagram> dd2a2(new Diagram(d2a2, 0.5));
  shared_ptr<Diagram> dd2a3(new Diagram(d2a3, 0.5));
  shared_ptr<Diagram> dd2a4(new Diagram(d2a4, 0.5));
  shared_ptr<Diagram> dd2a5(new Diagram(d2a5, 0.5));
  shared_ptr<Diagram> dd2a6(new Diagram(d2a6, 0.5));
  shared_ptr<Diagram> dd2a7(new Diagram(d2a7, 0.5));
  shared_ptr<Diagram> dd2a8(new Diagram(d2a8, 0.5));
  shared_ptr<Equation> ed2a0(new Equation(dd2a0, theory));
  shared_ptr<Equation> ed2a1(new Equation(dd2a1, theory));
  shared_ptr<Equation> ed2a2(new Equation(dd2a2, theory));
  shared_ptr<Equation> ed2a3(new Equation(dd2a3, theory));
  shared_ptr<Equation> ed2a4(new Equation(dd2a4, theory));
  shared_ptr<Equation> ed2a5(new Equation(dd2a5, theory));
  shared_ptr<Equation> ed2a6(new Equation(dd2a6, theory));
  shared_ptr<Equation> ed2a7(new Equation(dd2a7, theory));
  shared_ptr<Equation> ed2a8(new Equation(dd2a8, theory));
  ed2a0->merge(ed2a1);
  ed2a0->merge(ed2a2);
  ed2a0->merge(ed2a3);
  ed2a0->merge(ed2a4);
  ed2a0->merge(ed2a5);
  ed2a0->merge(ed2a6);
  ed2a0->merge(ed2a7);
  ed2a0->merge(ed2a8);
  ed2a0->duplicates();
  ed2a0->active();
  shared_ptr<Tree> td2a(new Density2(ed2a0, "density2"));

  list<shared_ptr<Tree>> trees = {tra, tea, tca, tda, td2a};
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
  cout << std::endl << "   ***  Energy  ***" << std::endl << std::endl;
  tea->print();
  cout << std::endl << "   ***  Correction  ***" << std::endl << std::endl;
  tca->print();
  cout << std::endl << "   ***  One-body Density Matrix  ***" << std::endl << std::endl;
  tda->print();
  cout << std::endl << "   ***  Two-body Density Matrix  ***" << std::endl << std::endl;
  td2a->print();
  cout << std::endl << std::endl;

  return 0;
}

