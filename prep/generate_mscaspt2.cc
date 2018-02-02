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
string theory = "MSCASPT2";

using namespace SMITH3::Prep;

tuple<vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>,
      vector<shared_ptr<Tensor>>, vector<shared_ptr<Tensor>>> create_proj() {
  vector<shared_ptr<Tensor>> lp, lt, ls, td, ltl, tdl;
  array<string, 3> label = {{"c", "x", "a"}};

  int cnt = 0;
  for (auto& i : label) {
    for (auto& j : label) {
      for (auto& k : label) {
        for (auto& l : label) {
          // full CASPT2
          if (
#if 1
              // all correct in this block
              (l == "c" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "c" && j == "a" && i == "a") ||
              (l == "x" && k == "x" && j == "a" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "a") ||
              (l == "c" && k == "c" && j == "x" && i == "x") ||
              (l == "x" && k == "c" && j == "x" && i == "x") ||
              (l == "x" && k == "x" && j == "x" && i == "a") ||
              (l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a")
#endif
            ) {
            stringstream ss; ss << cnt;
            lp.push_back(shared_ptr<Tensor>(new Tensor(ss.str(), {l, k, j, i})));
            td.push_back(shared_ptr<Tensor>(new Tensor("t2dagger", ss.str(), {l, k, j, i})));
            lt.push_back(shared_ptr<Tensor>(new Tensor("t2", ss.str(), {j, i, l, k})));
            tdl.push_back(shared_ptr<Tensor>(new Tensor("l2dagger", ss.str(), {l, k, j, i})));
            ltl.push_back(shared_ptr<Tensor>(new Tensor("l2", ss.str(), {j, i, l, k})));
            ++cnt;
          }
        }
      }
    }
  }

  return tie(lp, lt, td, ltl, tdl);
};

int main() {

  // generate common header
  cout << header() << endl;

  vector<shared_ptr<Tensor>> proj_list, t_list, t_dagger, l_list, l_dagger;
  tie(proj_list, t_list, t_dagger, l_list, l_dagger) = create_proj();

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
  for (auto& i : l_list)    cout << i->generate();
  for (auto& i : f)         cout << i->generate();
  for (auto& i : H)         cout << i->generate();
  for (auto& i : hc)        cout << i->generate();
  for (auto& i : dum)       cout << i->generate();
  for (auto& i : t_dagger)  cout << i->generate();
  for (auto& i : l_dagger)  cout << i->generate();
  for (auto& i : ex1b)      cout << i->generate();
  cout << endl;

  // density matrix equations //
  // one-body contribution d2
  shared_ptr<Equation> eq6(new Equation(theory, "da", {dum, t_dagger, ex1b, l_list}));
  eq6->set_tree_type("residual", "density");
  cout << eq6->generate();
  // one-body contribution d1
  shared_ptr<Equation> eq6a(new Equation(theory, "db", {dum, ex1b, l_list}));
  eq6a->set_tree_type("residual", "density1");
  cout << eq6a->generate();

  // two-body contribution D1
  shared_ptr<Equation> eq7(new Equation(theory, "d2a", {dum, proj_list, l_list}));
  eq7->set_tree_type("residual", "density2");
  cout << eq7->generate();

  // cI derivative equations, dedci = dE/dcI  //
  // test hylleraas eqn:   d/dc( <0|T^+fT|0> -e0<0|T^+T|0> +2<0|T^+h1|0> + 2<0|T^+V2|0>) =>
  //  =   1/2(1/4<I|T^+fT|0> + 1/4<0|T^+fT|I>) - 1/2*(e0/4<I|T^+T|0> + e0/4<0|T^+T|I>) + 2*1/2 (1/4<I|T^+V|0> + 1/4<0|T^+V|I>) + 2*1/2 (1/4<I|T^+h1|0> + 1/4<0|T^+h1|I>)
  // using bracket symmetry in some terms
  shared_ptr<Equation> eq4(new Equation(theory, "dedcia", {dum, l_dagger, f, t_list}, 1.0, make_pair(true, false)));
  eq4->set_tree_type("residual", "deci");
  cout << eq4->generate();

  shared_ptr<Equation> eq4b(new Equation(theory, "dedcic", {dum, l_dagger, t_list}, -1.0, "e0", make_pair(true, false)));
  eq4b->set_tree_type("residual", "deci2");
  cout << eq4b->generate();

  shared_ptr<Equation> eq4d(new Equation(theory, "dedcie", {dum, l_dagger, H}, 0.5, make_pair(true, false)));
  shared_ptr<Equation> eq4f(new Equation(theory, "dedcig", {dum, l_dagger, hc}, 1.0, make_pair(true, false)));
  eq4d->merge(eq4f);
  eq4d->set_tree_type("residual", "deci3");
  cout << eq4d->generate();

  shared_ptr<Equation> eq4e(new Equation(theory, "dedcif", {dum, l_dagger, H}, 0.5, make_pair(false, true)));
  shared_ptr<Equation> eq4g(new Equation(theory, "dedcih", {dum, l_dagger, hc}, 1.0, make_pair(false, true)));
  eq4e->merge(eq4g);
  eq4e->set_tree_type("residual", "deci4");
  cout << eq4e->generate();

  // done. generate the footer
  cout << footer("", "", "", eq6->tree_label(), eq6a->tree_label(), eq7->tree_label(), eq4->tree_label(), eq4b->tree_label(), eq4d->tree_label(), eq4e->tree_label()) << endl;


  return 0;
}


