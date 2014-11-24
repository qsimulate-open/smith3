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
  shared_ptr<Diagram> drb0(new Diagram(rb0, -1, "e0"));
  shared_ptr<Diagram> drb1(new Diagram(rb1, -1, "e0"));
  shared_ptr<Diagram> drb2(new Diagram(rb2, -1, "e0"));
  shared_ptr<Diagram> drb3(new Diagram(rb3, -1, "e0"));
  shared_ptr<Diagram> drb4(new Diagram(rb4, -1, "e0"));
  shared_ptr<Diagram> drb5(new Diagram(rb5, -1, "e0"));
  shared_ptr<Diagram> drb6(new Diagram(rb6, -1, "e0"));
  shared_ptr<Diagram> drb7(new Diagram(rb7, -1, "e0"));
  shared_ptr<Diagram> drb8(new Diagram(rb8, -1, "e0"));
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

  list<shared_ptr<Operator>> ea0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> ea1 = {proje, t2dagger0, f1, t21};
  list<shared_ptr<Operator>> ea2 = {proje, t2dagger0, f1, t22};
  list<shared_ptr<Operator>> ea3 = {proje, t2dagger1, f1, t20};
  list<shared_ptr<Operator>> ea4 = {proje, t2dagger1, f1, t21};
  list<shared_ptr<Operator>> ea5 = {proje, t2dagger1, f1, t22};
  list<shared_ptr<Operator>> ea6 = {proje, t2dagger2, f1, t20};
  list<shared_ptr<Operator>> ea7 = {proje, t2dagger2, f1, t21};
  list<shared_ptr<Operator>> ea8 = {proje, t2dagger2, f1, t22};
  list<shared_ptr<Operator>> eb0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> eb1 = {proje, t2dagger0, t21};
  list<shared_ptr<Operator>> eb2 = {proje, t2dagger0, t22};
  list<shared_ptr<Operator>> eb3 = {proje, t2dagger1, t20};
  list<shared_ptr<Operator>> eb4 = {proje, t2dagger1, t21};
  list<shared_ptr<Operator>> eb5 = {proje, t2dagger1, t22};
  list<shared_ptr<Operator>> eb6 = {proje, t2dagger2, t20};
  list<shared_ptr<Operator>> eb7 = {proje, t2dagger2, t21};
  list<shared_ptr<Operator>> eb8 = {proje, t2dagger2, t22};
  list<shared_ptr<Operator>> ec0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> ec1 = {proje, t2dagger1, v2};
  list<shared_ptr<Operator>> ec2 = {proje, t2dagger2, v2};
  list<shared_ptr<Operator>> ed0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> ed1 = {proje, t2dagger1, h1};
  list<shared_ptr<Operator>> ed2 = {proje, t2dagger2, h1};
  shared_ptr<Diagram> dea0(new Diagram(ea0, 0.25));
  shared_ptr<Diagram> dea1(new Diagram(ea1, 0.25));
  shared_ptr<Diagram> dea2(new Diagram(ea2, 0.25));
  shared_ptr<Diagram> dea3(new Diagram(ea3, 0.25));
  shared_ptr<Diagram> dea4(new Diagram(ea4, 0.25));
  shared_ptr<Diagram> dea5(new Diagram(ea5, 0.25));
  shared_ptr<Diagram> dea6(new Diagram(ea6, 0.25));
  shared_ptr<Diagram> dea7(new Diagram(ea7, 0.25));
  shared_ptr<Diagram> dea8(new Diagram(ea8, 0.25));
  shared_ptr<Diagram> deb0(new Diagram(eb0, -0.25, "e0"));
  shared_ptr<Diagram> deb1(new Diagram(eb1, -0.25, "e0"));
  shared_ptr<Diagram> deb2(new Diagram(eb2, -0.25, "e0"));
  shared_ptr<Diagram> deb3(new Diagram(eb3, -0.25, "e0"));
  shared_ptr<Diagram> deb4(new Diagram(eb4, -0.25, "e0"));
  shared_ptr<Diagram> deb5(new Diagram(eb5, -0.25, "e0"));
  shared_ptr<Diagram> deb6(new Diagram(eb6, -0.25, "e0"));
  shared_ptr<Diagram> deb7(new Diagram(eb7, -0.25, "e0"));
  shared_ptr<Diagram> deb8(new Diagram(eb8, -0.25, "e0"));
  shared_ptr<Diagram> dec0(new Diagram(ec0, 0.5));
  shared_ptr<Diagram> dec1(new Diagram(ec1, 0.5));
  shared_ptr<Diagram> dec2(new Diagram(ec2, 0.5));
  shared_ptr<Diagram> ded0(new Diagram(ed0, 0.5));
  shared_ptr<Diagram> ded1(new Diagram(ed1, 0.5));
  shared_ptr<Diagram> ded2(new Diagram(ed2, 0.5));
  shared_ptr<Equation> eea0(new Equation(dea0, theory));
  shared_ptr<Equation> eea1(new Equation(dea1, theory));
  shared_ptr<Equation> eea2(new Equation(dea2, theory));
  shared_ptr<Equation> eea3(new Equation(dea3, theory));
  shared_ptr<Equation> eea4(new Equation(dea4, theory));
  shared_ptr<Equation> eea5(new Equation(dea5, theory));
  shared_ptr<Equation> eea6(new Equation(dea6, theory));
  shared_ptr<Equation> eea7(new Equation(dea7, theory));
  shared_ptr<Equation> eea8(new Equation(dea8, theory));
  shared_ptr<Equation> eeb0(new Equation(deb0, theory));
  shared_ptr<Equation> eeb1(new Equation(deb1, theory));
  shared_ptr<Equation> eeb2(new Equation(deb2, theory));
  shared_ptr<Equation> eeb3(new Equation(deb3, theory));
  shared_ptr<Equation> eeb4(new Equation(deb4, theory));
  shared_ptr<Equation> eeb5(new Equation(deb5, theory));
  shared_ptr<Equation> eeb6(new Equation(deb6, theory));
  shared_ptr<Equation> eeb7(new Equation(deb7, theory));
  shared_ptr<Equation> eeb8(new Equation(deb8, theory));
  shared_ptr<Equation> eec0(new Equation(dec0, theory));
  shared_ptr<Equation> eec1(new Equation(dec1, theory));
  shared_ptr<Equation> eec2(new Equation(dec2, theory));
  shared_ptr<Equation> eed0(new Equation(ded0, theory));
  shared_ptr<Equation> eed1(new Equation(ded1, theory));
  shared_ptr<Equation> eed2(new Equation(ded2, theory));
  eea0->merge(eea1);
  eea0->merge(eea2);
  eea0->merge(eea3);
  eea0->merge(eea4);
  eea0->merge(eea5);
  eea0->merge(eea6);
  eea0->merge(eea7);
  eea0->merge(eea8);
  eea0->merge(eeb0);
  eea0->merge(eeb1);
  eea0->merge(eeb2);
  eea0->merge(eeb3);
  eea0->merge(eeb4);
  eea0->merge(eeb5);
  eea0->merge(eeb6);
  eea0->merge(eeb7);
  eea0->merge(eeb8);
  eea0->merge(eec0);
  eea0->merge(eec1);
  eea0->merge(eec2);
  eea0->merge(eed0);
  eea0->merge(eed1);
  eea0->merge(eed2);
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
  shared_ptr<Diagram> dda0(new Diagram(da0, 0.25));
  shared_ptr<Diagram> dda1(new Diagram(da1, 0.25));
  shared_ptr<Diagram> dda2(new Diagram(da2, 0.25));
  shared_ptr<Diagram> dda3(new Diagram(da3, 0.25));
  shared_ptr<Diagram> dda4(new Diagram(da4, 0.25));
  shared_ptr<Diagram> dda5(new Diagram(da5, 0.25));
  shared_ptr<Diagram> dda6(new Diagram(da6, 0.25));
  shared_ptr<Diagram> dda7(new Diagram(da7, 0.25));
  shared_ptr<Diagram> dda8(new Diagram(da8, 0.25));
  shared_ptr<Equation> eda0(new Equation(dda0, theory));
  shared_ptr<Equation> eda1(new Equation(dda1, theory));
  shared_ptr<Equation> eda2(new Equation(dda2, theory));
  shared_ptr<Equation> eda3(new Equation(dda3, theory));
  shared_ptr<Equation> eda4(new Equation(dda4, theory));
  shared_ptr<Equation> eda5(new Equation(dda5, theory));
  shared_ptr<Equation> eda6(new Equation(dda6, theory));
  shared_ptr<Equation> eda7(new Equation(dda7, theory));
  shared_ptr<Equation> eda8(new Equation(dda8, theory));
  eda0->merge(eda1);
  eda0->merge(eda2);
  eda0->merge(eda3);
  eda0->merge(eda4);
  eda0->merge(eda5);
  eda0->merge(eda6);
  eda0->merge(eda7);
  eda0->merge(eda8);
  eda0->duplicates();
  eda0->active();
  shared_ptr<Tree> tda(new Density(eda0, "density"));

  list<shared_ptr<Operator>> db0 = {proje, ex_1b, t20};
  list<shared_ptr<Operator>> db1 = {proje, ex_1b, t21};
  list<shared_ptr<Operator>> db2 = {proje, ex_1b, t22};
  shared_ptr<Diagram> ddb0(new Diagram(db0));
  shared_ptr<Diagram> ddb1(new Diagram(db1));
  shared_ptr<Diagram> ddb2(new Diagram(db2));
  shared_ptr<Equation> edb0(new Equation(ddb0, theory));
  shared_ptr<Equation> edb1(new Equation(ddb1, theory));
  shared_ptr<Equation> edb2(new Equation(ddb2, theory));
  edb0->merge(edb1);
  edb0->merge(edb2);
  edb0->duplicates();
  edb0->active();
  shared_ptr<Tree> tdb(new Density1(edb0, "density1"));

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

  list<shared_ptr<Operator>> dedcia0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcia1 = {proje, t2dagger0, f1, t21};
  list<shared_ptr<Operator>> dedcia2 = {proje, t2dagger0, f1, t22};
  list<shared_ptr<Operator>> dedcia3 = {proje, t2dagger1, f1, t20};
  list<shared_ptr<Operator>> dedcia4 = {proje, t2dagger1, f1, t21};
  list<shared_ptr<Operator>> dedcia5 = {proje, t2dagger1, f1, t22};
  list<shared_ptr<Operator>> dedcia6 = {proje, t2dagger2, f1, t20};
  list<shared_ptr<Operator>> dedcia7 = {proje, t2dagger2, f1, t21};
  list<shared_ptr<Operator>> dedcia8 = {proje, t2dagger2, f1, t22};
  list<shared_ptr<Operator>> dedcib0 = {proje, t2dagger0, f1, t20};
  list<shared_ptr<Operator>> dedcib1 = {proje, t2dagger0, f1, t21};
  list<shared_ptr<Operator>> dedcib2 = {proje, t2dagger0, f1, t22};
  list<shared_ptr<Operator>> dedcib3 = {proje, t2dagger1, f1, t20};
  list<shared_ptr<Operator>> dedcib4 = {proje, t2dagger1, f1, t21};
  list<shared_ptr<Operator>> dedcib5 = {proje, t2dagger1, f1, t22};
  list<shared_ptr<Operator>> dedcib6 = {proje, t2dagger2, f1, t20};
  list<shared_ptr<Operator>> dedcib7 = {proje, t2dagger2, f1, t21};
  list<shared_ptr<Operator>> dedcib8 = {proje, t2dagger2, f1, t22};
  list<shared_ptr<Operator>> dedcic0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcic1 = {proje, t2dagger0, t21};
  list<shared_ptr<Operator>> dedcic2 = {proje, t2dagger0, t22};
  list<shared_ptr<Operator>> dedcic3 = {proje, t2dagger1, t20};
  list<shared_ptr<Operator>> dedcic4 = {proje, t2dagger1, t21};
  list<shared_ptr<Operator>> dedcic5 = {proje, t2dagger1, t22};
  list<shared_ptr<Operator>> dedcic6 = {proje, t2dagger2, t20};
  list<shared_ptr<Operator>> dedcic7 = {proje, t2dagger2, t21};
  list<shared_ptr<Operator>> dedcic8 = {proje, t2dagger2, t22};
  list<shared_ptr<Operator>> dedcid0 = {proje, t2dagger0, t20};
  list<shared_ptr<Operator>> dedcid1 = {proje, t2dagger0, t21};
  list<shared_ptr<Operator>> dedcid2 = {proje, t2dagger0, t22};
  list<shared_ptr<Operator>> dedcid3 = {proje, t2dagger1, t20};
  list<shared_ptr<Operator>> dedcid4 = {proje, t2dagger1, t21};
  list<shared_ptr<Operator>> dedcid5 = {proje, t2dagger1, t22};
  list<shared_ptr<Operator>> dedcid6 = {proje, t2dagger2, t20};
  list<shared_ptr<Operator>> dedcid7 = {proje, t2dagger2, t21};
  list<shared_ptr<Operator>> dedcid8 = {proje, t2dagger2, t22};
  list<shared_ptr<Operator>> dedcie0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcie1 = {proje, t2dagger1, v2};
  list<shared_ptr<Operator>> dedcie2 = {proje, t2dagger2, v2};
  list<shared_ptr<Operator>> dedcif0 = {proje, t2dagger0, v2};
  list<shared_ptr<Operator>> dedcif1 = {proje, t2dagger1, v2};
  list<shared_ptr<Operator>> dedcif2 = {proje, t2dagger2, v2};
  list<shared_ptr<Operator>> dedcig0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> dedcig1 = {proje, t2dagger1, h1};
  list<shared_ptr<Operator>> dedcig2 = {proje, t2dagger2, h1};
  list<shared_ptr<Operator>> dedcih0 = {proje, t2dagger0, h1};
  list<shared_ptr<Operator>> dedcih1 = {proje, t2dagger1, h1};
  list<shared_ptr<Operator>> dedcih2 = {proje, t2dagger2, h1};
  shared_ptr<Diagram> ddedcia0(new Diagram(dedcia0, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia1(new Diagram(dedcia1, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia2(new Diagram(dedcia2, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia3(new Diagram(dedcia3, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia4(new Diagram(dedcia4, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia5(new Diagram(dedcia5, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia6(new Diagram(dedcia6, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia7(new Diagram(dedcia7, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcia8(new Diagram(dedcia8, 0.25, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcib0(new Diagram(dedcib0, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib1(new Diagram(dedcib1, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib2(new Diagram(dedcib2, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib3(new Diagram(dedcib3, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib4(new Diagram(dedcib4, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib5(new Diagram(dedcib5, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib6(new Diagram(dedcib6, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib7(new Diagram(dedcib7, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcib8(new Diagram(dedcib8, 0.25, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcic0(new Diagram(dedcic0, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic1(new Diagram(dedcic1, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic2(new Diagram(dedcic2, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic3(new Diagram(dedcic3, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic4(new Diagram(dedcic4, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic5(new Diagram(dedcic5, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic6(new Diagram(dedcic6, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic7(new Diagram(dedcic7, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcic8(new Diagram(dedcic8, -0.25, "e0", std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcid0(new Diagram(dedcid0, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid1(new Diagram(dedcid1, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid2(new Diagram(dedcid2, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid3(new Diagram(dedcid3, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid4(new Diagram(dedcid4, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid5(new Diagram(dedcid5, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid6(new Diagram(dedcid6, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid7(new Diagram(dedcid7, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcid8(new Diagram(dedcid8, -0.25, "e0", std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcie0(new Diagram(dedcie0, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcie1(new Diagram(dedcie1, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcie2(new Diagram(dedcie2, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcif0(new Diagram(dedcif0, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcif1(new Diagram(dedcif1, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcif2(new Diagram(dedcif2, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcig0(new Diagram(dedcig0, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcig1(new Diagram(dedcig1, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcig2(new Diagram(dedcig2, 0.5, std::make_pair(true, false)));
  shared_ptr<Diagram> ddedcih0(new Diagram(dedcih0, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcih1(new Diagram(dedcih1, 0.5, std::make_pair(false, true)));
  shared_ptr<Diagram> ddedcih2(new Diagram(dedcih2, 0.5, std::make_pair(false, true)));
  shared_ptr<Equation> ededcia0(new Equation(ddedcia0, theory));
  shared_ptr<Equation> ededcia1(new Equation(ddedcia1, theory));
  shared_ptr<Equation> ededcia2(new Equation(ddedcia2, theory));
  shared_ptr<Equation> ededcia3(new Equation(ddedcia3, theory));
  shared_ptr<Equation> ededcia4(new Equation(ddedcia4, theory));
  shared_ptr<Equation> ededcia5(new Equation(ddedcia5, theory));
  shared_ptr<Equation> ededcia6(new Equation(ddedcia6, theory));
  shared_ptr<Equation> ededcia7(new Equation(ddedcia7, theory));
  shared_ptr<Equation> ededcia8(new Equation(ddedcia8, theory));
  shared_ptr<Equation> ededcib0(new Equation(ddedcib0, theory));
  shared_ptr<Equation> ededcib1(new Equation(ddedcib1, theory));
  shared_ptr<Equation> ededcib2(new Equation(ddedcib2, theory));
  shared_ptr<Equation> ededcib3(new Equation(ddedcib3, theory));
  shared_ptr<Equation> ededcib4(new Equation(ddedcib4, theory));
  shared_ptr<Equation> ededcib5(new Equation(ddedcib5, theory));
  shared_ptr<Equation> ededcib6(new Equation(ddedcib6, theory));
  shared_ptr<Equation> ededcib7(new Equation(ddedcib7, theory));
  shared_ptr<Equation> ededcib8(new Equation(ddedcib8, theory));
  shared_ptr<Equation> ededcic0(new Equation(ddedcic0, theory));
  shared_ptr<Equation> ededcic1(new Equation(ddedcic1, theory));
  shared_ptr<Equation> ededcic2(new Equation(ddedcic2, theory));
  shared_ptr<Equation> ededcic3(new Equation(ddedcic3, theory));
  shared_ptr<Equation> ededcic4(new Equation(ddedcic4, theory));
  shared_ptr<Equation> ededcic5(new Equation(ddedcic5, theory));
  shared_ptr<Equation> ededcic6(new Equation(ddedcic6, theory));
  shared_ptr<Equation> ededcic7(new Equation(ddedcic7, theory));
  shared_ptr<Equation> ededcic8(new Equation(ddedcic8, theory));
  shared_ptr<Equation> ededcid0(new Equation(ddedcid0, theory));
  shared_ptr<Equation> ededcid1(new Equation(ddedcid1, theory));
  shared_ptr<Equation> ededcid2(new Equation(ddedcid2, theory));
  shared_ptr<Equation> ededcid3(new Equation(ddedcid3, theory));
  shared_ptr<Equation> ededcid4(new Equation(ddedcid4, theory));
  shared_ptr<Equation> ededcid5(new Equation(ddedcid5, theory));
  shared_ptr<Equation> ededcid6(new Equation(ddedcid6, theory));
  shared_ptr<Equation> ededcid7(new Equation(ddedcid7, theory));
  shared_ptr<Equation> ededcid8(new Equation(ddedcid8, theory));
  shared_ptr<Equation> ededcie0(new Equation(ddedcie0, theory));
  shared_ptr<Equation> ededcie1(new Equation(ddedcie1, theory));
  shared_ptr<Equation> ededcie2(new Equation(ddedcie2, theory));
  shared_ptr<Equation> ededcif0(new Equation(ddedcif0, theory));
  shared_ptr<Equation> ededcif1(new Equation(ddedcif1, theory));
  shared_ptr<Equation> ededcif2(new Equation(ddedcif2, theory));
  shared_ptr<Equation> ededcig0(new Equation(ddedcig0, theory));
  shared_ptr<Equation> ededcig1(new Equation(ddedcig1, theory));
  shared_ptr<Equation> ededcig2(new Equation(ddedcig2, theory));
  shared_ptr<Equation> ededcih0(new Equation(ddedcih0, theory));
  shared_ptr<Equation> ededcih1(new Equation(ddedcih1, theory));
  shared_ptr<Equation> ededcih2(new Equation(ddedcih2, theory));
  ededcia0->merge(ededcia1);
  ededcia0->merge(ededcia2);
  ededcia0->merge(ededcia3);
  ededcia0->merge(ededcia4);
  ededcia0->merge(ededcia5);
  ededcia0->merge(ededcia6);
  ededcia0->merge(ededcia7);
  ededcia0->merge(ededcia8);
  ededcia0->merge(ededcib0);
  ededcia0->merge(ededcib1);
  ededcia0->merge(ededcib2);
  ededcia0->merge(ededcib3);
  ededcia0->merge(ededcib4);
  ededcia0->merge(ededcib5);
  ededcia0->merge(ededcib6);
  ededcia0->merge(ededcib7);
  ededcia0->merge(ededcib8);
  ededcia0->merge(ededcic0);
  ededcia0->merge(ededcic1);
  ededcia0->merge(ededcic2);
  ededcia0->merge(ededcic3);
  ededcia0->merge(ededcic4);
  ededcia0->merge(ededcic5);
  ededcia0->merge(ededcic6);
  ededcia0->merge(ededcic7);
  ededcia0->merge(ededcic8);
  ededcia0->merge(ededcid0);
  ededcia0->merge(ededcid1);
  ededcia0->merge(ededcid2);
  ededcia0->merge(ededcid3);
  ededcia0->merge(ededcid4);
  ededcia0->merge(ededcid5);
  ededcia0->merge(ededcid6);
  ededcia0->merge(ededcid7);
  ededcia0->merge(ededcid8);
  ededcia0->merge(ededcie0);
  ededcia0->merge(ededcie1);
  ededcia0->merge(ededcie2);
  ededcia0->merge(ededcif0);
  ededcia0->merge(ededcif1);
  ededcia0->merge(ededcif2);
  ededcia0->merge(ededcig0);
  ededcia0->merge(ededcig1);
  ededcia0->merge(ededcig2);
  ededcia0->merge(ededcih0);
  ededcia0->merge(ededcih1);
  ededcia0->merge(ededcih2);
  ededcia0->absorb_ket();
  ededcia0->duplicates();
  ededcia0->active();
  shared_ptr<Tree> tdedcia(new Dedci(ededcia0, "dedci"));

  list<shared_ptr<Tree>> trees = {tra, tea, tca, tda, tdb, td2a, tdedcia};
  shared_ptr<Forest> fr(new Forest(trees));

  fr->filter_gamma();
  list<shared_ptr<Tensor>> gm = fr->gamma();
  const list<shared_ptr<Tensor>> gamma = gm;

  tuple<string, string, string> tmp = fr->generate_code();

  ofstream fs(fr->name() + ".h");
  ofstream es(fr->name() + "_tasks.h");
  ofstream cs(fr->name() + "_tasks.cc");
  fs << get<0>(tmp);
  es << get<1>(tmp);
  cs << get<2>(tmp);
  fs.close();
  es.close();
  cs.close();
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

