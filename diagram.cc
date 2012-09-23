//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: diagram.cc
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


#include "diagram.h"

using namespace std;


list<shared_ptr<Diagram> > Diagram::get_all() const {
  list<shared_ptr<Diagram> > out;

  const int max = 1 << num_general();
  for (int i = 0; i != max; ++i) {
    int j = i;
    shared_ptr<Diagram> d = copy();
    for (auto& k : d->op()) k->mutate_general(j);
    if (d->consistent_indices()) {
      out.push_back(d);
    }
  }
  return out;
}

// copying diagram with the same connectivity and so on.
shared_ptr<Diagram> Diagram::copy() const {
  // mapping of indices and spins. map<old, new>
  map<shared_ptr<Index>, shared_ptr<Index> > indexmap;
  map<shared_ptr<Spin>, shared_ptr<Spin> > spinmap;

  // creates Diagram without any info
  shared_ptr<Diagram> out(new Diagram());
  list<shared_ptr<Op> > outop;

  // loop over operators
  for (auto& it : op()) {
    // cloning...
    shared_ptr<Op> a = it->copy();

    auto j = a->op().begin();
    for (auto i = it->op().begin(); i != it->op().end(); ++i, ++j) {
      auto s = indexmap.find(*get<0>(*i));
      if (s == indexmap.end()) {
        indexmap.insert(make_pair(*get<0>(*i), *get<0>(*j)));
      } else {
        *get<0>(*j) = s->second; // one of the operators has s->second...
      }
      auto z = spinmap.find(it->rho(get<2>(*i)));
      if (z == spinmap.end()) {
        spinmap.insert(make_pair(it->rho(get<2>(*i)), a->rho(get<2>(*j))));
      } else {
        a->set_rho(get<2>(*j), z->second); // one of the operators has z->second...
      }
      get<1>(*j) = get<1>(*i);
    }
    outop.push_back(a);
  }
  out->set_op(outop);
  out->set_fac(fac_);
  if (dagger_) out->add_dagger();
  return out;
}


void Diagram::refresh_indices() {
  map<shared_ptr<Index>, int> dict;
  map<shared_ptr<Index>, int> done;
  map<shared_ptr<Spin>, int> spin;
  for (auto& i : op_)
    i->refresh_indices(dict, done, spin);
}


// this is not a const function because it refreshes the indices
void Diagram::print() {
  refresh_indices();
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto& i : op_) i->print();

  // active operators
  cout << "[";
  for (auto& i : op_) {
    if (i->num_active_nodagger() + i->num_active_dagger() != 0) {
      for (auto& j : i->op()) {
        if (get<1>(j) == -1 || get<1>(j) == 0) continue;
        cout << (*get<0>(j))->str();
      }
    }
  }
  cout << "]";
  if (dagger_) cout << " ** Daggered object added **";

  cout << endl;
  if (rdm_) rdm_->print("   ");
  cout << endl;
}


list<shared_ptr<Index> > Diagram::active_indices() const {
  list<shared_ptr<Index> > out;
  for (auto& i : op_) {
    if (i->num_active_nodagger() + i->num_active_dagger() != 0) {
      for (auto& j : i->op())
        if (get<1>(j) == 2) out.push_back(*get<0>(j));
    }
  }
  return out;
}


void Diagram::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto& i : op_) i->print();
  cout << endl;
}


bool Diagram::reduce_one_noactive(const int skip) {
  refresh_indices();

  bool found = false;
  // find the first dagger operator in list<Op>
  auto i = op_.begin();
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > data; // safe because they are held by tensors
  for (; i != op_.end(); ++i) {
    // this simultaneously eliminates one entry. op_ is therefore modified here
    data = (*i)->first_dagger_noactive();
    if (data.first) break;
  }
  if (!data.first) return false;

  // skip until it comes to skip
  int cnt = 0;
  shared_ptr<Spin> newspin, oldspin;
  for (auto j = op_.begin(); j != op_.end(); ++j) {
    // cannot contract with self
    if (i == j) continue;
    // all possible contraction pattern taken for *j (returned as a list).
    if (cnt + (*j)->num_nodagger() > skip) {
      tuple<double,shared_ptr<Spin>,shared_ptr<Spin> > tmp = (*j)->contract(data, skip-cnt);
      fac_ *= get<0>(tmp);
      newspin = get<1>(tmp);
      oldspin = get<2>(tmp);
      found = true;
      break;
    } else {
      cnt += (*j)->num_nodagger();
    }
  }
  for (auto& j : op_) {
    for (auto& k : j->rho()) {
      if (k == oldspin) k = newspin;
    } 
  }
  return found;
}


bool Diagram::valid() const {
  int out = 0;
  for (auto& i : op_) {
    if (!i->contracted()) ++out;
  }
  return out > 1;
}


bool Diagram::done() const {
  int out = 0;
  for (auto& i : op_) {
    if (!i->contracted()) ++out;
  }
  return out == 0; 
}


bool Diagram::done_noactive() const {
  int out = 0;
  for (auto& i : op_) {
    out += i->num_nodagger() + i->num_dagger();
  }
  return out == 0; 
}


int Diagram::num_dagger() const {
  int out = 0;
  for (auto& i : op_) out += i->num_dagger();
  return out;
}


int Diagram::num_general() const {
  int cnt = 0;
  for (auto& i : op_) cnt += i->num_general(); 
  return cnt;
}


bool Diagram::consistent_indices() const {
  int cnt1 = 0;
  int cnt2 = 0;
  for (auto& i : op_) {
    cnt1 += i->num_active_dagger(); 
    cnt2 += i->num_active_nodagger(); 
  }
  return cnt1 == cnt2;
}

void Diagram::active() {
  // Index should be updated.
  refresh_indices();
  list<shared_ptr<Index> >  ac = active_indices();
  if (ac.size()) {
    // Performs Wick in constructor of an Active object
    rdm_ = shared_ptr<Active>(new Active(ac));
  }
}


bool Diagram::permute(const bool proj) {
  bool found = false;
 
  // try the last one. If it returns zero, then go up.
  for (auto i = op_.rbegin(); i != op_.rend(); ++i) {
    pair<bool, double> a = (*i)->permute(proj);
    fac_ *= a.second; 
    if (!a.first) {
      continue;
    } else {
      found = true;
      refresh_indices();
      break;
    }
  }
  return found;
}




bool Diagram::identical(shared_ptr<Diagram> o) const {

  bool out = true;
  // first, they should be same size
  if (op_.size() != o->op().size()) out = false;
  // second, each indices should be the same (spin is not checked here)
  if (out) {
    for (auto i = op_.begin(), j = o->op().begin(); i != op_.end(); ++i, ++j) {
      out &= (*i)->identical(*j);
    }
  }
  // then, we check spins.
  if (out) {
    list<shared_ptr<Index> > act = active_indices();
    list<shared_ptr<Index> > oact = o->active_indices();
    map<shared_ptr<Spin>, shared_ptr<Spin> > myo;
    if (act.size() != oact.size()) {
      out = false;
    } else {
      for (auto i = act.begin(), j = oact.begin(); i != act.end(); ++i, ++j) {
        assert((*i)->identical(*j)); 
        shared_ptr<Spin> s = (*i)->spin();
        shared_ptr<Spin> os = (*j)->spin();
        auto iter = myo.find(s);
        if (myo.end() == iter) {
          // if s appears for the first time, register it
          myo.insert(make_pair(s,os));
        } else {
          if (os != iter->second) {
            out =false;
            break;
          }
        }
      }
    }
  } 

  return out;
}



