//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: operator.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the SMITH3 package.
//
// The SMITH3 package is free software; you can redistribute it and\/or modify
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


#include "operator.h"
#include <algorithm>

using namespace std;
using namespace smith;

Operator::Operator(const string& ta, const string& tb)
  : a_(make_shared<Index>(ta,true)), b_(make_shared<Index>(tb,false)) {
  op_.push_back(make_tuple(&a_, ta!="x"?0:2, 0)); // index, dagger, spin
  op_.push_back(make_tuple(&b_, tb!="x"?0:2, 0));
  rho_.push_back(make_shared<Spin>());

  perm_.push_back(0);
}


Operator::Operator(const string& ta, const string& tb, const string& tc, const string& td)
  : a_(make_shared<Index>(ta,true)), b_(make_shared<Index>(tb,true)), c_(make_shared<Index>(tc,false)), d_(make_shared<Index>(td,false)) {
  // accept aa,ii and rearrange it to ai,ai
  op_.push_back(make_tuple(&a_, ta!="x"?0:2, 0)); // index, no-active/active, spin
  op_.push_back(make_tuple(&d_, td!="x"?0:2, 0)); // from historical reasons, it is 0 and 2. -1 when contracted.
  op_.push_back(make_tuple(&b_, tb!="x"?0:2, 1));
  op_.push_back(make_tuple(&c_, tc!="x"?0:2, 1));

  rho_.push_back(make_shared<Spin>());
  rho_.push_back(make_shared<Spin>());

  perm_.push_back(0);
  perm_.push_back(1);

}

// Constructing one-body tensor. No operator is created. Careful...
Operator::Operator(shared_ptr<Index> ta, shared_ptr<Index> tb, shared_ptr<Spin> ts)
 : a_(ta), b_(tb) {
  rho_.push_back(ts); // just to prevent seg fault.
  op_.push_back(make_tuple(&a_, -1, 0));
  op_.push_back(make_tuple(&b_, -1, 0));
}



int Operator::num_nodagger() const {
  int out = 0;
  for (auto& i : op_)
    if (get<1>(i)==0 && !(*get<0>(i))->dagger()) ++out;
  return out;
}


int Operator::num_dagger() const {
  int out = 0;
  for (auto& i : op_)
    if (get<1>(i)==0 && (*get<0>(i))->dagger()) ++out;
  return out;
}


int Operator::num_active_dagger() const {
  int out = 0;
  for (auto& i : op_)
    if (get<1>(i)==2 && (*get<0>(i))->dagger()) ++out;
  return out;
}


int Operator::num_active_nodagger() const {
  int out = 0;
  for (auto& i : op_)
    if (get<1>(i)==2 && !(*get<0>(i))->dagger()) ++out;
  return out;
}


bool Operator::general() const {
  return num_general() != 0;
}


int Operator::num_general() const {
  int out = 0;
  for (auto& i : op_)
    if((*get<0>(i))->label() ==  "g") ++out;
  return out;
}


void Operator::mutate_general(int& in) {
  for (auto& i : op_) {
    if ((*get<0>(i))->label() ==  "g") {
      if (in & 1) {
        (*get<0>(i))->set_label("x");
        get<1>(i) += 2;
      }
      in >>= 1;  // decrease in by one bit
    }
  }
}


bool Operator::contracted() const {
  int out = 0;
  for (auto& i : op_)
    if (get<1>(i) == 0) ++out;
  return out == 0;
}


void Operator::refresh_indices(map<shared_ptr<const Index>, int>& dict,
                         map<shared_ptr<const Index>, int>& done,
                         map<shared_ptr<Spin>, int>& spin) {
  //
  // Note: seperate labeling for those still in the operators and those
  //       already contracted. This is to make it easy to get the minus sign in the
  //       Wick theorem evaluator.
  //
  for (auto& i : op_) {
    // if this is not still contracted
    if (get<1>(i) != -1) {
      auto iter = dict.find(*get<0>(i));
      if (iter == dict.end()) {
        const int c = dict.size();
        dict.insert(make_pair(*get<0>(i), c));
        (*get<0>(i))->set_num(c);
      }
    // if this is already contracted, we use negative values (does not have to be, though - just for print out)
    } else {
      auto iter = done.find(*get<0>(i));
      if (iter == done.end()) {
        const int c = done.size();
        done.insert(make_pair(*get<0>(i), -c-1));
        (*get<0>(i))->set_num(-c-1);
      }
    }

    auto ster = spin.find(rho(get<2>(i)));
    if (get<1>(i) != -1) {
      if (ster == spin.end()) {
        const int c = spin.size();
        spin.insert(make_pair(rho(get<2>(i)), c));
        rho(get<2>(i))->set_num(c);
      }
    }

    // set all the spins into operators
    (*get<0>(i))->set_spin(rho(get<2>(i)));
  }
}


pair<shared_ptr<Index>*, shared_ptr<Spin>* > Operator::first_dagger_noactive() {
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > out;
  auto i = op_.begin();
  for (; i != op_.end(); ++i) {
    if (get<1>(*i)==0 && (*get<0>(*i))->dagger()) { // "x" is active orbitals = 2
      out = make_pair(get<0>(*i), rho_ptr(get<2>(*i)));
      break;
    }
  }
  if (out.first) get<1>(*i) = -1;
  return out;
}


shared_ptr<Index>* Operator::survive(shared_ptr<Index>* a, shared_ptr<Index>* b) {
  string alab = (*a)->label();
  string blab = (*b)->label();
  if (alab == blab) return a;
  else if (alab == "g" && blab != "g") return b;
  else if (alab != "g" && blab == "g") return a;
  else throw logic_error("A strange thing happened in Op::survive");
};


tuple<double, shared_ptr<Spin>, shared_ptr<Spin>>
     Operator::contract(pair<shared_ptr<Index>*, shared_ptr<Spin>* >& dat, const int skip) {
  int cnt = 0;
  auto i = op_.begin();
  double fac = 0.0;
  shared_ptr<Spin> a, b;
  for (; i != op_.end(); ++i) {
    if (get<1>(*i)!=0 || (*get<0>(*i))->dagger()) continue;
    if (contractable((*get<0>(*i))->label(), (*dat.first)->label())) {
      if (cnt == skip) {
        const int n1 = (*dat.first)->num();
        const int n2 = (*get<0>(*i))->num();
        if (n1 == n2) throw logic_error("Should not happen. Op::contract");
        fac = (abs(n2-n1) & 1) ? 1.0 : -1.0;

        *get<0>(*i) = *survive(get<0>(*i), dat.first);
        *dat.first = *get<0>(*i);
        fac *= (*dat.second == rho(get<2>(*i))) ? 2.0 : 1.0;
        a = *dat.second;
        b = rho(get<2>(*i));
        set_rho(get<2>(*i), *dat.second);

        break;
      } else {
        ++cnt;
      }
    }
  }
  get<1>(*i) = -1;
  return make_tuple(fac, a, b);
}




