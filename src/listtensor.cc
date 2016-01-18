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


#include <iomanip>
#include <algorithm>
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

  // add a ci tensor if braket and if no rdm derivatives. This tensor is the overlap, cI coefficients.
  if ((d->braket().first || d->braket().second) && !d->rdm()) {
    list<shared_ptr<const Index>> in = d->target_index();
    shared_ptr<Tensor> t = make_shared<Tensor>(fac_,"dci",in);
    list_.push_back(t);
  }

  // add rdm tensors.
  // the rdm ci derivatives have an extra index
  if (d->rdm() && d->braket().first) {  // bra
    list<shared_ptr<const Index>> in = d->target_index();
    shared_ptr<Tensor> t = make_shared<Tensor>(d->rdm(), in);
    list_.push_back(t);
  } else if (d->rdm() && d->braket().second) { // ket
    list<shared_ptr<const Index>> in = d->target_index();
    shared_ptr<Tensor> t = make_shared<Tensor>(d->rdm(), in, d->rdm()->num_map());
    list_.push_back(t);
  } else if (d->rdm() && (!d->braket().first && !d->braket().second)) { // add normal rdm tensor
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
    // original reindexing from rdm derivative
    map<int, int> rdm_num_map;
    // get modified rdm indices. these will be reversed in associated tensors
    list<shared_ptr<const Index>> ind;
    for (auto i = list_.begin(); i != list_.end(); ++i) {
      if ((*i)->is_gamma()) {
        rdm_num_map = (*i)->num_map();
        list<shared_ptr<const Index>> tmp;
        tmp = (*i)->index();
        for (auto& j : tmp) {  // check to make sure not ci target index
          if (j->label() != "ci"){
            ind.push_back(j);
          }
        }
      }
    }

    for (auto i = list_.begin(); i != list_.end(); ++i) {
      list<shared_ptr<const Index>> newind;
      if (!(*i)->is_gamma() && !ind.empty() && (*i)->label() != "proj") {
        for (auto& j : (*i)->index()) {
          if (!j->active()) {
            newind.push_back(j);
          } else {
            shared_ptr<const Index> tmp = make_shared<const Index>(*j,rdm_num_map[j->num()]);
            newind.push_back(tmp);
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
          if ((*j)->label() == "ci") break;   // todo is there a better way?
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
  // make intermediate tensor
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
        if (braket().first || braket().second) cout << "dcI_" << i->str();
        else cout << i->str();
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


shared_ptr<Cost> ListTensor::calculate_cost() const {
  auto out = make_shared<Cost>();
  list<shared_ptr<const Index>> current = list_.back()->index();

  for (auto i = ++list_.rbegin(); i != list_.rend(); ++i) {
    list<shared_ptr<const Index>> sumindex, outindex;
    for (auto& a : current)
      for (auto& b : (*i)->index())
        if (a->same_num(b) && a->same_label(b))
          sumindex.push_back(a);

    current.insert(current.end(), (*i)->index().begin(), (*i)->index().end());
    for (auto& a : current) {
      bool check = false;
      for (auto& b : sumindex)
        if (a->same_num(b) && a->same_label(b))
          check = true;
      if (!check)
        outindex.push_back(a);
    }

    sumindex.insert(sumindex.end(), outindex.begin(), outindex.end());
    vector<int> cost(4);
    for (auto& a : sumindex) {
      if (a->label() == "c") cost[0] += 1;
      else if (a->label() == "x") cost[1] += 1;
      else if (a->label() == "a") cost[2] += 1;
      else if (a->label() == "ci") cost[3] += 1;
      else {
        stringstream ss; ss << "this should not happen - ListTensor::calculate_cost " << a->label() << endl;
        throw logic_error(ss.str());
      }
    }
    auto pcost = make_shared<PCost>(cost);
    out->add_pcost(*pcost);

    current = outindex;
  }

  out->sort_pcost();
  assert(list_.size()-1 == out->cost().size());
  return out;
}


void ListTensor::reorder() {
  // I need to sort list_ first
  vector<shared_ptr<Tensor>> tmp(list_.begin(), list_.end());
  sort(tmp.begin(), tmp.end());
  list<shared_ptr<Tensor>> out(tmp.begin(), tmp.end());
  list_ = out;

  shared_ptr<Cost> current;
  do {
    shared_ptr<Cost> cost = calculate_cost();
    if (!current || *cost < *current) {
      out = list_;
      current = cost;
    }
  } while (next_permutation(list_.begin(), list_.end()));
  list_ = out;
}


