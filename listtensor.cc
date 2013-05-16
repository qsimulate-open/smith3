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
  // vector of tensors
  for (auto& i : d->op()) {
    // careful, add only labeled operators! ie not excitation operators!
    if (!i->label().empty()) {
      shared_ptr<Tensor> t(new Tensor(i));
      list_.push_back(t);
    }
  }
  if (d->rdm()) {
    shared_ptr<Tensor> t(new Tensor(d->rdm()));
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
  shared_ptr<Tensor> t(new Tensor(1.0, ss.str(), ind));
  return t;
}


shared_ptr<ListTensor> ListTensor::rest() const {
  list<shared_ptr<Tensor>> r = list_; 
  r.pop_front();     
  shared_ptr<ListTensor> out(new ListTensor(fac_, scalar_, r, dagger_));
  return out;
}


void ListTensor::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << (scalar_.empty() ? "" : " * "+ scalar_) << " ";
  for (auto& i : list_) cout << i->str();
  if (dagger_) cout << " ** Daggered object added **";
  cout << endl;
  for (auto& i : list_) {
    if (i->active()) i->active()->print("   ");
  }
  cout << endl;
}


