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
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
string theory = "RelCASPT2";

using namespace SMITH3::Prep;

tuple<vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>> create_proj() {
  vector<shared_ptr<Tensor>> lp, lt, ls, td;
  array<string, 3> label = {{"c", "x", "a"}};

  int cnt = 0;
  for (auto& i : label) {
    for (auto& j : label) {
      for (auto& k : label) {
        for (auto& l : label) {
          // full CASPT2
          if (
              // all correct in this block
              (l == "c" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "x" && j == "a" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "x") ||
              (l == "x" && k == "c" && j == "x" && i == "x") ||
              (l == "x" && k == "x" && j == "x" && i == "a") ||
              (l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a")
            ) {
            stringstream ss; ss << cnt;
            lp.push_back(shared_ptr<Tensor>(new Tensor(ss.str(), {l, k, j, i})));
            td.push_back(shared_ptr<Tensor>(new Tensor("t2dagger", ss.str(), {l, k, j, i})));
            lt.push_back(shared_ptr<Tensor>(new Tensor("t2", ss.str(), {j, i, l, k})));
            ++cnt;
          }
        }
      }
    }
  }

  return tie(lp, lt, td);
};

int main() {

  // generate common header
  cout << header() << endl;

  vector<shared_ptr<Tensor>> proj_list, t_list, t_dagger;
  tie(proj_list, t_list, t_dagger) = create_proj();

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
  for (auto& i : f)         cout << i->generate();
  for (auto& i : H)         cout << i->generate();
  for (auto& i : hc)        cout << i->generate();
  for (auto& i : dum)       cout << i->generate();
  for (auto& i : t_dagger)  cout << i->generate();
  for (auto& i : ex1b)      cout << i->generate();
  cout << endl;

  // residual equations //
  shared_ptr<Equation> eq0(new Equation(theory, "ra", {dum, proj_list, f, t_list}));
  shared_ptr<Equation> eq1(new Equation(theory, "rb", {dum, proj_list, t_list}, -1.0, "e0"));
  eq0->merge(eq1);
  eq0->set_tree_type("residual");
  cout << eq0->generate();

  // energy equations //
  // second order energy correction
  // S = <proj|H|0>. <R|T> will be added in bagel
  shared_ptr<Equation> eq3(new Equation(theory, "ec", {dum, proj_list, H}, 0.5));
  shared_ptr<Equation> eq3a(new Equation(theory, "ed", {dum, proj_list, hc}));
  eq3->merge(eq3a);
  eq3->set_tree_type("residual", "source");
  cout << eq3->generate();

  // generate Norm <1|1> to be used in various places
  shared_ptr<Equation> eq5(new Equation(theory, "ca", {dum, proj_list, t_list}));
  eq5->set_tree_type("residual", "norm");
  cout << eq5->generate();

  // done. generate the footer
  cout << footer_ci(eq0->tree_label(), eq3->tree_label(), eq5->tree_label()) << endl;

  return 0;
}


