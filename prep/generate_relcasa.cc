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

string theory = "RelCASA";

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
              (l == "c" && k == "c" && j == "a" && i == "a") ||                                                  // closed->virtual, closed->virtual
              (l == "x" && k == "c" && j == "a" && i == "a") ||                                                  // active->virtual, closed->virtual
              (l == "x" && k == "x" && j == "a" && i == "a") ||                                                  // active->virtual, active->virtual
              (l == "c" && k == "c" && j == "x" && i == "a") ||                                                  // closed->active,  closed->virtual
              (l == "c" && k == "c" && j == "x" && i == "x") ||                                                  // closed->active,  closed->active
              (l == "x" && k == "c" && j == "x" && i == "x") ||                                                  // active->active,  closed->active
              (l == "x" && k == "x" && j == "x" && i == "a") ||                                                  // active->active,  active->virtual
              (l == "c" && k == "x" && j == "x" && i == "a") || (l == "x" && k == "c" && j == "x" && i == "a")   // closed->active,  active->virtual (or equiv)
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
  vector<shared_ptr<Tensor>> f   = {shared_ptr<Tensor>(new Tensor("f1", "0", {"c", "c"})),  // closed,  closed
                                    shared_ptr<Tensor>(new Tensor("f1", "1", {"c", "x"})),  // closed,  active
                                    shared_ptr<Tensor>(new Tensor("f1", "2", {"x", "c"})),  // active,  closed
                                    shared_ptr<Tensor>(new Tensor("f1", "3", {"c", "a"})),  // closed,  virtual
                                    shared_ptr<Tensor>(new Tensor("f1", "4", {"a", "c"})),  // virtual, closed
                                    shared_ptr<Tensor>(new Tensor("f1", "5", {"x", "a"})),  // active,  virtual
                                    shared_ptr<Tensor>(new Tensor("f1", "6", {"a", "x"})),  // virtual, active
                                    shared_ptr<Tensor>(new Tensor("f1", "7", {"a", "a"}))   // virtual, virtual
                                   };

  vector<shared_ptr<Tensor>> hc  = {shared_ptr<Tensor>(new Tensor("h1", "0", {"g", "g"}))};
  vector<shared_ptr<Tensor>> H   = {shared_ptr<Tensor>(new Tensor("v2", "0", {"g", "g", "g", "g"}))};
  vector<shared_ptr<Tensor>> dum = {shared_ptr<Tensor>(new Tensor("proj", "e", {}))};
  vector<shared_ptr<Tensor>> ex1b = {shared_ptr<Tensor>(new Tensor("1b", {"g", "g"}))};

  vector<shared_ptr<Tensor>> hca  = {shared_ptr<Tensor>(new Tensor("h1", "1", {"x", "x"}))};
  vector<shared_ptr<Tensor>> Ha   = {shared_ptr<Tensor>(new Tensor("v2", "1", {"x", "x", "x", "x"}))};

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
  for (auto& i : hca)       cout << i->generate();
  for (auto& i : Ha)        cout << i->generate();
  cout << endl;

  // residual equations //
  // <Omega| F T |0> for all except active-active part 
  shared_ptr<Equation> eq0(new Equation(theory, "ra", {dum, proj_list, f, t_list}));

  // For matching sectors:  -<Omega| T |N> * [E_L^(0) - E_N^(0)] due to the use of commutator to avoid 5RDM
  // For unmatched sectors:  -<Omega| T |N> * E_L^(0), but these terms are zero so we can get away with just one constant "e0"
  shared_ptr<Equation> eq1(new Equation(theory, "rb", {dum, proj_list, t_list}, -1.0, "e0"));
  eq0->merge(eq1);

  // <Omega| H T |0> for active part
  shared_ptr<Equation> eq0x(new Equation(theory, "rax", {dum, proj_list, hca, t_list}));
  shared_ptr<Equation> eq1x(new Equation(theory, "rbx", {dum, proj_list, Ha, t_list}, 0.5));
  eq0->merge(eq0x);
  eq0->merge(eq1x);

  // - <Omega| T F(A) |0> for matching sectors, to generate [F(A), T]  (F(A) = H for active-active part, F for all others)
  for (int i = 0; i != proj_list.size(); ++i) {
    if (i == 3 || i == 4) continue;
    stringstream ss, tt, uu;
    ss << "rax_" << i;
    tt << "rbx_" << i;
    uu << "rcc_" << i;
    shared_ptr<Equation> eq0m(new Equation(theory, ss.str(), {dum, vector<shared_ptr<Tensor>>{proj_list[i]}, vector<shared_ptr<Tensor>>{t_list[i]}, hca}, -1.0));
    shared_ptr<Equation> eq1m(new Equation(theory, tt.str(), {dum, vector<shared_ptr<Tensor>>{proj_list[i]}, vector<shared_ptr<Tensor>>{t_list[i]}, Ha},  -0.5));
    shared_ptr<Equation> eq2m(new Equation(theory, uu.str(), {dum, vector<shared_ptr<Tensor>>{proj_list[i]}, vector<shared_ptr<Tensor>>{t_list[i]}, f}, -1.0));
    eq0->merge(eq0m);
    eq0->merge(eq1m);
    eq0->merge(eq2m);
  }
  shared_ptr<Equation> eq0m(new Equation(theory, "rax_3", {dum, vector<shared_ptr<Tensor>>{proj_list[3], proj_list[4]}, vector<shared_ptr<Tensor>>{t_list[3], t_list[4]}, hca}, -1.0));
  shared_ptr<Equation> eq1m(new Equation(theory, "rbx_3", {dum, vector<shared_ptr<Tensor>>{proj_list[3], proj_list[4]}, vector<shared_ptr<Tensor>>{t_list[3], t_list[4]}, Ha},  -0.5));
  shared_ptr<Equation> eq2m(new Equation(theory, "rcc_3", {dum, vector<shared_ptr<Tensor>>{proj_list[3], proj_list[4]}, vector<shared_ptr<Tensor>>{t_list[3], t_list[4]}, f},  -1.0));
  eq0->merge(eq0m);
  eq0->merge(eq1m);
  eq0->merge(eq2m);

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

  // done. generate the footer
  cout << footer(eq0->tree_label(), eq3->tree_label(), "") << endl;

  return 0;
}


