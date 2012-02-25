//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: op.cc
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


#include "op.h"
#include <algorithm>

using namespace std;

Op::Op(const std::string lab, const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td)
  : label_(lab), a_(new Index(ta,true)), b_(new Index(tb,true)), c_(new Index(tc,false)), d_(new Index(td,false)) {
  op_.push_back(std::make_tuple(&a_, ta!="x"?0:2, 0)); // index, no-active/active, spin
  op_.push_back(std::make_tuple(&d_, td!="x"?0:2, 0)); // from historical reasons, it is 0 and 2. -1 when contracted.
  op_.push_back(std::make_tuple(&b_, tb!="x"?0:2, 1));
  op_.push_back(std::make_tuple(&c_, tc!="x"?0:2, 1));

  std::shared_ptr<Spin> tmp(new Spin());
  rho_.push_back(tmp);

  std::shared_ptr<Spin> tmp2(new Spin());
  rho_.push_back(tmp2);

  perm_.push_back(0);
  perm_.push_back(1);

}


Op::Op(const std::string lab, const std::string& ta, const std::string& tb)
  : label_(lab), a_(new Index(ta,true)), b_(new Index(tb,false)) {
  op_.push_back(std::make_tuple(&a_, ta!="x"?0:2, 0)); // index, dagger, spin
  op_.push_back(std::make_tuple(&b_, tb!="x"?0:2, 0));
  std::shared_ptr<Spin> tmp(new Spin());
  rho_.push_back(tmp);

  perm_.push_back(0);
}


// Constructing one-body tensor. No operator is created. Careful...
Op::Op(const string lab, shared_ptr<Index> ta, shared_ptr<Index> tb, shared_ptr<Spin> ts)
 : label_(lab), a_(ta), b_(tb) {
  rho_.push_back(ts); // just to prevent seg fault.
  op_.push_back(std::make_tuple(&a_, -1, 0));
  op_.push_back(std::make_tuple(&b_, -1, 0));
}


shared_ptr<Op> Op::copy() const {
  // in the case of two-body operators
  if (c_) {
    shared_ptr<Op> tmp(new Op(label_, a_->label(), b_->label(), c_->label(), d_->label()));
    return tmp;
  } else  {
    shared_ptr<Op> tmp(new Op(label_, a_->label(), b_->label()));
    return tmp;
  }
}

int Op::num_nodagger() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)==0 && !(*get<0>(*i))->dagger()) ++out;
  return out;
}


int Op::num_dagger() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)==0 && (*get<0>(*i))->dagger()) ++out;
  return out;
}


int Op::num_active_dagger() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)==2 && (*get<0>(*i))->dagger()) ++out;
  return out;
}


int Op::num_active_nodagger() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)==2 && !(*get<0>(*i))->dagger()) ++out;
  return out;
}


bool Op::general() const {
  return num_general() != 0;
}


int Op::num_general() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if((*get<0>(*i))->label() ==  "g") ++out; 
  return out;
}


void Op::mutate_general(int& in) {
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    if ((*get<0>(*i))->label() ==  "g") {
      if (in & 1) {
        (*get<0>(*i))->set_label("x");
        get<1>(*i) += 2;
      }
      in >>= 1;
    }
  } 
}


bool Op::contracted() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i) == 0) ++out;
  return out == 0;
}


void Op::print() const {
  // tensor
  cout << label_ << "(";
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    cout << (*get<0>(*i))->str(false) << " ";
  }
  cout << ") ";
  // (non active) operator
  if (num_nodagger() + num_dagger() != 0) {
    cout << "{";
    for (auto i = op_.begin(); i != op_.end(); ++i) {
      if (get<1>(*i) == -1 || get<1>(*i) == 2) continue;
      cout << (*get<0>(*i))->str() << rho(get<2>(*i))->str();
    }
    cout << "} ";
  }
}


void Op::refresh_indices(map<shared_ptr<Index>, int>& dict,
                         map<shared_ptr<Index>, int>& done,
                         map<shared_ptr<Spin>, int>& spin) {
  //
  // Note: seperate labeling for those still in the operators and those
  //       already contracted. This is to make it easy to get the minus sign in the
  //       Wick theorem evaluator.
  //
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    // if this is not still contracted
    if (get<1>(*i) != -1) {
      auto iter = dict.find(*get<0>(*i));
      if (iter == dict.end()) {
        const int c = dict.size();
        dict.insert(make_pair(*get<0>(*i), c));
        (*get<0>(*i))->set_num(c);
      }
    // if this is already contracted, we use negative values (does not have to be, though - just for print out)
    } else {
      auto iter = done.find(*get<0>(*i));
      if (iter == done.end()) {
        const int c = done.size();
        done.insert(make_pair(*get<0>(*i), -c-1));
        (*get<0>(*i))->set_num(-c-1);
      }
    }

    auto ster = spin.find(rho(get<2>(*i)));
    if (get<1>(*i) != -1) {
      if (ster == spin.end()) {
        const int c = spin.size();
        spin.insert(make_pair(rho(get<2>(*i)), c));
        rho(get<2>(*i))->set_num(c);
      }
    }

    // set all the spins into operators
    (*get<0>(*i))->set_spin(rho(get<2>(*i)));
  }
}


pair<shared_ptr<Index>*, shared_ptr<Spin>* > Op::first_dagger_noactive() {
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > out;
  auto i = op_.begin();
  for (; i != op_.end(); ++i) {
    if (get<1>(*i)==0 && (*get<0>(*i))->dagger()) { // "x" is active orbitals
      out = make_pair(get<0>(*i), rho_ptr(get<2>(*i)));
      break;
    }
  }
  if (out.first) get<1>(*i) = -1;
  return out;
}


shared_ptr<Index>* Op::survive(shared_ptr<Index>* a, shared_ptr<Index>* b) {
  string alab = (*a)->label();
  string blab = (*b)->label();
  if (alab == blab) return a;
  else if (alab == "g" && blab != "g") return b;
  else if (alab != "g" && blab == "g") return a;
  else throw logic_error("A strange thing happened in Op::survive");
};


tuple<double, shared_ptr<Spin>, shared_ptr<Spin> >
     Op::contract(pair<shared_ptr<Index>*, shared_ptr<Spin>* >& dat, const int skip) {
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
        if (n1 == n2) throw logic_error("should not happen. Op::contract");
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


// this function make a possible permutation of indices.
pair<bool, double> Op::permute(const bool proj) {
  // if there is active daggered and no-daggered operators, you cannot do this as it changes the expression
  if (num_active_nodagger() && num_active_dagger() || (!proj && label_ == "proj"))
    return make_pair(false, 1.0);

  const vector<int> prev = perm_;
  const int size = prev.size();
  bool next = next_permutation(perm_.begin(), perm_.end());

  vector<int> map(size);
  vector<int> imap(size); // inverse map
  for (int i = 0; i != size; ++i) {
    const int ii = prev.at(i);
    for (int j = 0; j != size; ++j) {
      const int jj = perm_.at(j);
      if (ii == jj) {
        map[i] = j;
        imap[j] = i;
        break;
      }
    }
  }
  // permute indices (I don't remember why op_ is list...)
  vector<tuple<shared_ptr<Index>*, int, int> > tmp(size*2);
  auto oiter = op_.begin();
  for (int i = 0; i != size; ++i, ++oiter) {
    tmp[map[i]*2  ] = *oiter; ++oiter; 
    tmp[map[i]*2+1] = *oiter; 
  }

  // find sign.
  vector<int> act(size);
  oiter = op_.begin();
  for (int i = 0; i != size; ++i, ++oiter) {
    if ((*get<0>(*oiter))->active()) act[i] += 1; ++oiter; 
    if ((*get<0>(*oiter))->active()) act[i] += 1;
  }
  int f = 0;
  for (int i = 0; i != size; ++i) {
    const int inew = i;
    // if #active op is either zero or two, the sign won't change.
    if (act[imap[i]] != 1) continue;
    // operators originally in the left side && target is larger than mine
    for (int j = 0; j != imap[i]; ++j) {
      if (inew < map[j]) f += act[j];
    }
  }
  const double out = (f&1) ? -1.0 : 1.0;

  oiter = op_.begin();
  for (auto t = tmp.begin(); t != tmp.end(); ++t, ++oiter) *oiter = *t; 

  return make_pair(next, out);
}



// Note that spin info is not checked. 
bool Op::identical(shared_ptr<Op> o) const {
  bool out = true;
  auto j = o->op().begin();
  for (auto i = op_.begin(); i != op_.end(); ++i, ++j) {
    if (get<1>(*i) != get<1>(*j) || !(*get<0>(*i))->identical(*get<0>(*j))) {
      out = false;
      break;
    }
  }
  return out;
}


