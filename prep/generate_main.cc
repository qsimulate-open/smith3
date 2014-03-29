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

#include <iostream>
#include <tuple>
#include <string>
#include <memory>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <array>
#include <sstream>
#include <initializer_list>

#include "constants.h"
#include "tensor.h"
#include "diagram.h"
#include "equation.h"

using namespace std;

//const string theory = "MP2";
string theory = "CAS_test";

using namespace SMITH3::Prep;

tuple<vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>> create_proj() {
  vector<shared_ptr<Tensor>> lp, lt, ls, td;
  array<string, 3> label = {{"c", "x", "a"}};

  int cnt = 0;
  for (auto& i : label) {
    for (auto& j : label) {
      for (auto& k : label) {
        for (auto& l : label) {
#if 0     // full CASPT2
          if ((l == "c" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "x" && j == "a" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "x") ||
              (l == "c" && k == "x" && j == "x" && i == "a") ||
              (l == "x" && k == "c" && j == "x" && i == "a") ||
              (l == "x" && k == "x" && j == "x" && i == "a") ||
              (l == "x" && k == "c" && j == "x" && i == "x")) {
#else  // turn on one of the following lines
       // if ((l == "c" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "x" && j == "a" && i == "a")) {  // mp2 test ansatz
       // if ((l == "c" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "x" && j == "a" && i == "a") || (l == "c" && k == "c" && j == "x" && i == "a")) {
       // if (l == "c" && k == "x" && j == "x" && i == "a") {  // test cxxa
       // if (l == "x" && k == "c" && j == "x" && i == "a") {  // test xcxa
       // if ((l == "c" && k == "c" && j == "x" && i == "a") || (l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a") || (l == "x" && k == "x" && j == "x" && i == "a")) {  // semi-internal cases
       // if ((l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a") || (l == "x" && k == "x" && j == "x" && i == "a")) {  // y4 internal test
       // if ((l == "c" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "c" && j == "a" && i == "a") || (l == "x" && k == "x" && j == "a" && i == "a") ||  // y5 cases
       //     (l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a") || (l == "x" && k == "x" && j == "x" && i == "a")) {
       // if (l == "x" && k == "c" && j == "x" && i == "a") { // xcxa
// *test single configuration cases*
          if (l == "c" && k == "c" && j == "a" && i == "a") { // ccaa
//        if (l == "x" && k == "c" && j == "a" && i == "a") { // xcaa
//        if (l == "x" && k == "x" && j == "a" && i == "a") { // xxaa
//        if (l == "c" && k == "c" && j == "x" && i == "a") { // ccxa
//        if ((l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a")) { // cxxa or xcxa
//        if (l == "c" && k == "c" && j == "x" && i == "x") { // ccxx
//        if (l == "x" && k == "x" && j == "x" && i == "a") { // xxxa
//        if (l == "x" && k == "c" && j == "x" && i == "x") { // xcxx
// *end test single configuration cases*
#endif
            stringstream ss; ss << cnt;
            lp.push_back(shared_ptr<Tensor>(new Tensor(ss.str(), {l, k, j, i})));
            td.push_back(shared_ptr<Tensor>(new Tensor("t2dagger", ss.str(), {l, k, j, i})));
            lt.push_back(shared_ptr<Tensor>(new Tensor("t2", ss.str(), {j, i, l, k})));
            ls.push_back(shared_ptr<Tensor>(new Tensor("r", ss.str(), {j, i, l, k})));
            ++cnt;
          }
        }
      }
    }
  }

  return tie(lp, lt, ls, td);
};

int main() {

  // generate common header
  cout << header() << endl;

  vector<shared_ptr<Tensor>> proj_list, t_list, r_list, t_dagger;
  tie(proj_list, t_list, r_list, t_dagger) = create_proj();

  // make f and H tensors here
  vector<shared_ptr<Tensor>> f   = {shared_ptr<Tensor>(new Tensor("f1", "", {"g", "g"}))};
  vector<shared_ptr<Tensor>> hc  = {shared_ptr<Tensor>(new Tensor("h1", "", {"g", "g"}))};
  vector<shared_ptr<Tensor>> H   = {shared_ptr<Tensor>(new Tensor("v2", "", {"g", "g", "g", "g"}))};
  vector<shared_ptr<Tensor>> dum = {shared_ptr<Tensor>(new Tensor("proj", "e", {}))};
  vector<shared_ptr<Tensor>> ex1b = {shared_ptr<Tensor>(new Tensor("1b", {"g", "g"}))};

  cout << "  string theory=\"" << theory << "\";" << endl;
  cout << endl;

  for (auto& i : proj_list) cout << i->generate();
  for (auto& i : t_list)    cout << i->generate();
  for (auto& i : r_list)    cout << i->generate();
  for (auto& i : f)         cout << i->generate();
  for (auto& i : H)         cout << i->generate();
  for (auto& i : hc)        cout << i->generate();
  for (auto& i : dum)       cout << i->generate();
  for (auto& i : t_dagger)  cout << i->generate();
  for (auto& i : ex1b)      cout << i->generate();
  cout << endl;

  // residual equations //
  shared_ptr<Equation> eq0(new Equation("ra", {dum, proj_list, f, t_list}));
  shared_ptr<Equation> eq1(new Equation("rb", {dum, proj_list, t_list}, -1.0, "e0"));
  shared_ptr<Equation> eq2(new Equation("rc", {dum, proj_list, H}));
  shared_ptr<Equation> eq2a(new Equation("rd", {dum, proj_list, hc}));
  eq0->merge(eq1);
  eq0->merge(eq2);
  eq0->merge(eq2a);
  eq0->set_tree_type("residual");
  cout << eq0->generate();

  // energy equations //
 // second order energy correction
  shared_ptr<Equation> eq3(new Equation("ea", {dum, t_dagger, H}, 0.25));
  shared_ptr<Equation> eq3a(new Equation("eb", {dum, t_dagger, hc}, 0.25));
  eq3->merge(eq3a);
  eq3->set_tree_type("energy");
  cout << eq3->generate();

  // generate Norm <1|1> to be used in various places, y correction and correction term for one-body density matrix
  shared_ptr<Equation> eq5(new Equation("ca", {dum, t_dagger, t_list}, 0.25));
  eq5->set_tree_type("correction");
  cout << eq5->generate();

  // density matrix equations //
  // one-body contribution
  shared_ptr<Equation> eq6(new Equation("da", {dum, t_dagger, ex1b, t_list}, 0.25));
  shared_ptr<Equation> eq6a(new Equation("db", {dum, ex1b, t_list}, 1.0));
  eq6->merge(eq6a);
  eq6->set_tree_type("density");
  cout << eq6->generate();

  // two-body contribution
  shared_ptr<Equation> eq7(new Equation("d2a", {dum, proj_list, t_list}, 0.5));
  eq7->set_tree_type("density2");
  cout << eq7->generate();

  // cI derivative equations, dedci = dE/dcI  //
  // test energy, NB: for testing use 1/2 correct prefactor (scale). This should equate dedci with energy, otherwise expect 2E.
  double scale = 0.5;
  // test hylleraas eqn:   d/dc( <0|T^+fT|0> -e0<0|T^+T|0> +2<0|T^+V2|0> + 2<0|T^+V2|0>) =>
  //  =   1/2(1/4<I|T^+fT|0> + 1/4<0|T^+fT|I>) - 1/2*(e0/4<I|T^+T|0> + e0/4<0|T^+T|I>) + 2*1/2 (1/4<I|T^+V|0> + 1/4<0|T^+V|I>) + 2*1/2 (1/4<I|T^+h1|0> + 1/4<0|T^+h1|I>)
  shared_ptr<Equation> eq4(new Equation("dedcia", {dum, t_dagger, f, t_list}, scale*0.25, make_pair(true, false)));
  shared_ptr<Equation> eq4a(new Equation("dedcib", {dum, t_dagger, f, t_list}, scale*0.25, make_pair(false, true)));
  shared_ptr<Equation> eq4b(new Equation("dedcic", {dum, t_dagger, t_list}, -1.0*scale*0.25, "e0", make_pair(true, false)));
  shared_ptr<Equation> eq4c(new Equation("dedcid", {dum, t_dagger, t_list}, -1.0*scale*0.25, "e0", make_pair(false, true)));
  shared_ptr<Equation> eq4d(new Equation("dedcie", {dum, t_dagger, H}, scale*0.50, make_pair(true, false)));
  shared_ptr<Equation> eq4e(new Equation("dedcif", {dum, t_dagger, H}, scale*0.50, make_pair(false, true)));
  shared_ptr<Equation> eq4f(new Equation("dedcig", {dum, t_dagger, hc}, scale*0.50, make_pair(true, false)));
  shared_ptr<Equation> eq4g(new Equation("dedcih", {dum, t_dagger, hc}, scale*0.50, make_pair(false, true)));
  eq4->merge(eq4a);
  eq4->merge(eq4b);
  eq4->merge(eq4c);
  eq4->merge(eq4d);
  eq4->merge(eq4e);
  eq4->merge(eq4f);
  eq4->merge(eq4g);
  eq4->set_tree_type("dedci");
  cout << eq4->generate();

  // done. generate the footer
  cout << footer(eq0->tree_label(), eq3->tree_label(), eq5->tree_label(), eq6->tree_label(), eq7->tree_label(), eq4->tree_label()) << endl;


  return 0;
}


