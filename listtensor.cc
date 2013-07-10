//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: listtensor.cc
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


#include "listtensor.h"

using namespace std;
using namespace smith;

ListTensor::ListTensor(shared_ptr<Diagram> d) {
  // factor
  fac_ = d->fac();
  scalar_ = d->scalar();
  braket_ = d->braket();
  // vector of tensors
  for (auto& i : d->op()) {
    // Careful, add only labeled operators! ie not excitation operators!
    if (!i->label().empty()) {
      shared_ptr<Tensor> t = make_shared<Tensor>(i);
      list_.push_back(t);
    }
  }
  if (d->rdm()) {
    shared_ptr<Tensor> t = make_shared<Tensor>(d->rdm());
    list_.push_back(t);
  }
  // dagger
  dagger_ = d->dagger();
}


void ListTensor::absorb_all_internal() {
  auto j = list_.begin();
  // first find active
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    if ((*i)->active()) j = i;
  }
  list<list<shared_ptr<Tensor>>::iterator> remove;
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    if ((*i)->all_active() && !(*i)->active() && (*i)->label() != "proj") {
      (*j)->merge(*i);
      remove.push_back(i);
    }
  }
  for (auto i = remove.begin(); i != remove.end(); ++i) list_.erase(*i);
}


void ListTensor::absorb_ket() {
  if (braket_.second) {
    assert(!braket_.first);
    // get modified rdm indices. these will be reversed in associated tensors
    list<shared_ptr<const Index>> ind;
    for (auto i = list_.begin(); i != list_.end(); ++i) {
      if ((*i)->is_gamma()) {
        ind = (*i)->index();
      }
    }

    // map indices to reverse
    list<shared_ptr<const Index>> rev_ind(ind.rbegin(),ind.rend());
    map<shared_ptr<const Index>, shared_ptr<const Index>> ind_map;
    for (auto i = ind.begin(), j = rev_ind.begin(); i != ind.end() ; ++i, ++j) ind_map[(*i)] = (*j);

    for (auto i = list_.begin(); i != list_.end(); ++i) {
      list<shared_ptr<const Index>> newind;
      if (!(*i)->is_gamma() && !ind.empty() && (*i)->label() != "proj") {
        for (auto& j : (*i)->index()) {
          if (!j->active()) {
            newind.push_back(j);
          } else {
            newind.push_back(ind_map[j]);
          }
        }
        (*i)->set_index(newind);
      }
    }
    // now braket can be reversed for this listtensor
    set_braket(make_pair(true,false));

    // reverse braket for modified rdm
    for (auto& i : list_) {
      if (i->active()) {
        list<shared_ptr<RDM>> rdms = i->active()->rdm(); 
        for (auto& j : rdms) {
          if (j->ket()) {
            j->set_bra(true); 
            j->set_ket(false); 
          } 
        } 
      }
    }
     
  }
}


bool ListTensor::has_gamma() const {
 bool found = false;
 for (auto& i : list_) {
   if (i->is_gamma()) found = true;
 }
 return found;
}


static int target_num__;
shared_ptr<Tensor> ListTensor::target() const {
  list<shared_ptr<const Index>> ind;
  for (auto t = list_.begin(); t != list_.end(); ++t) {
    list<shared_ptr<const Index>> index = (*t)->index();
    for (auto j = index.begin(); j != index.end(); ++j) {
      bool found = false;
      list<shared_ptr<const Index>>::iterator remove;
      for (auto i = ind.begin(); i != ind.end(); ++i) {
        if ((*i)->num() == (*j)->num()) {
          found = true;
          remove = i;
          break;
        }
      }
      if (found) {
        ind.erase(remove);
      } else {
        ind.push_back(*j);
      }
    }
  }
  stringstream ss;
  ss << "I" << target_num__;
  ++target_num__;
  shared_ptr<Tensor> t = make_shared<Tensor>(1.0, ss.str(), ind);
  return t;
}


shared_ptr<ListTensor> ListTensor::rest() const {
  list<shared_ptr<Tensor>> r = list_;
  r.pop_front();
  shared_ptr<ListTensor> out = make_shared<ListTensor>(fac_, scalar_, r, dagger_, braket_);
  return out;
}


void ListTensor::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << (scalar_.empty() ? "" : " * "+ scalar_) << " ";
  size_t found = false;
  for (auto& i : list_) {
    if (i->str().find("Gamma") != string::npos) found = true;
  }
  if (!found) {
    for (auto i = list_.begin(); i != list_.end(); ++i) {
      cout << (*i)->str();
      if (i == --list_.end()) cout << " <" << (braket_.first ? "I" : "0") << "|" << (braket_.second ? "I" : "0") << ">";
    }
  } else {
    for (auto& i : list_) {
      if (i->str().find("Gamma") != string::npos) {
        cout << "CI_" << i->str();
      } else {
        cout << i->str();
      }
    }
  }
  if (dagger_) cout << " ** Daggered object added **";
  cout << endl;
  for (auto& i : list_) {
    if (i->active()) i->active()->print("   ");
  }
  cout << endl;
}


