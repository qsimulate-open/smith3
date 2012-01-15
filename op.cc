//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "op.h"

using namespace std;

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
    if (get<1>(*i)==1) ++out;
  return out;
}


int Op::num_dagger() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)==0) ++out;
  return out;
}


bool Op::contracted() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    if (get<1>(*i)!=-1) ++out;
  return out == 0;
}

void Op::print() const {
  // tensor
  cout << label_ << "(";
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    cout << (*get<0>(*i))->str() << " ";
  }
  cout << ") {";
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    if (get<1>(*i) < 0) continue;
    cout << (*get<0>(*i))->str() << (get<1>(*i)==0 ? "+" : "") << (*get<2>(*i))->str();
  }
  cout << "}";
}


void Op::refresh_indices(map<shared_ptr<Index>, int>& dict,
                         map<shared_ptr<Index>, int>& done,
                         map<shared_ptr<Spin>, int>& spin) {
  //
  // Note: the labeling is independent for those still in the operators and those
  //       already contracted. This is to make it easy to get the minus sign in the
  //       Wick theorem evaluator.
  //
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    // if this is not still contracted
    if (get<1>(*i) >= 0) {
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

    auto ster = spin.find(*get<2>(*i));
    if (get<1>(*i) >= 0) {
      if (ster == spin.end()) {
        const int c = spin.size();
        spin.insert(make_pair(*get<2>(*i), c));
        (*get<2>(*i))->set_num(c);
      }
    }
  }
}


pair<shared_ptr<Index>*, shared_ptr<Spin>* > Op::first_dagger_noactive() {
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > out;
  auto i = op_.begin();
  for (; i != op_.end(); ++i) {
    if (get<1>(*i)==0 && (*get<0>(*i))->label() != "x") { // "x" is active orbitals
      out = make_pair(get<0>(*i), get<2>(*i));
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
  else throw runtime_error("strange in survive");
};


double Op::contract(pair<shared_ptr<Index>*, shared_ptr<Spin>* >& dat, const int skip) {
  int cnt = 0;
  auto i = op_.begin();
  double fac = 0.0;
  for (; i != op_.end(); ++i) {
    if (get<1>(*i)!=1) continue;
    if (contractable((*get<0>(*i))->label(), (*dat.first)->label())) {
      if (cnt == skip) {
        const int n1 = (*dat.first)->num();
        const int n2 = (*get<0>(*i))->num();
        if (n1 == n2) throw logic_error("should not happen. Op::contract");
        fac = (abs(n2-n1) & 1) ? 1.0 : -1.0;

        *get<0>(*i) = *survive(get<0>(*i), dat.first);
        *dat.first = *get<0>(*i);
        fac *= (*dat.second == *get<2>(*i)) ? 2.0 : 1.0;
        *dat.second = *get<2>(*i);
        break;
      } else {
        ++cnt;
      }
    }
  }
  get<1>(*i) = -1;
  return fac;
}

