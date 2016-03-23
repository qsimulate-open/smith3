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

string theory = "SPCASPT2";

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

  vector<shared_ptr<Tensor>> t_list, t_dagger;
  tie(ignore, t_list, t_dagger) = create_proj();

  // make f and H tensors here
  vector<shared_ptr<Tensor>> dum = {shared_ptr<Tensor>(new Tensor("proj", "e", {}))};
  vector<shared_ptr<Tensor>> ex1b = {shared_ptr<Tensor>(new Tensor("1b", {"g", "g"}, /*alpha*/true))};

  cout << "  string theory=\"" << theory << "\";" << endl;
  cout << endl;

  for (auto& i : t_list)    cout << i->generate();
  for (auto& i : dum)       cout << i->generate();
  for (auto& i : t_dagger)  cout << i->generate();
  for (auto& i : ex1b)      cout << i->generate();
  cout << endl;

  // density matrix equations //
  // one-body contribution d2
  shared_ptr<Equation> eq6(new Equation(theory, "da", {dum, t_dagger, ex1b, t_list}));
  eq6->set_tree_type("residual", "density");
  cout << eq6->generate();
  // one-body contribution d1
  shared_ptr<Equation> eq6a(new Equation(theory, "db", {dum, ex1b, t_list}));
  eq6a->set_tree_type("residual", "density1");
  cout << eq6a->generate();

  // done. generate the footer
  cout << footer("", "", "", eq6->tree_label(), eq6a->tree_label(), "", "") << endl;


  return 0;
}


