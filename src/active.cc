//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: active.cc
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


#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include "active.h"

using namespace std;
using namespace smith;



Active::Active(const list<shared_ptr<const Index>>& in, pair<bool,bool> braket) : bra_(braket.first), ket_(braket.second) {
  shared_ptr<RDM> tmp;
  if (!braket.first && !braket.second) {
    tmp = make_shared<RDM00>(in, map<shared_ptr<const Index>, shared_ptr<const Index>>(), braket, 1.0);
  } else if (braket.first || braket.second) {
    // Caution, braket is passed directly so both modified rdms <I|E|0> and <0|E|I> are made here.
    tmp = make_shared<RDMI0>(in, map<shared_ptr<const Index>, shared_ptr<const Index>>(), braket, 1.0);
  } else if (braket.first && braket.second) {
    throw logic_error("Active::ctor not implemented");
  }
  // this sets list<RDM>
  reduce(tmp);

}


void Active::reduce(shared_ptr<RDM> in) {

  if (ket_) {  // in this case we should take conjugate and reindex
    list<shared_ptr<const Index>> cindex = in->conjugate();
    // map
    shared_ptr<RDM> stmp = in->copy();
    list<shared_ptr<const Index>> sindex = stmp->index();
    map<int, int> num_map;
    list<shared_ptr<const Index>> rev_ind(sindex.rbegin(), sindex.rend());
    for (auto i = sindex.begin(), j = rev_ind.begin(); i != sindex.end(); ++i, ++j) num_map[(*i)->num()] = (*j)->num();
    num_map_ = num_map;
    // reindex
    list<shared_ptr<const Index>> new_index;
    for (auto i : cindex) {
      shared_ptr<const Index> tmp = make_shared<const Index>((*i), num_map[i->num()]);
      new_index.push_back(tmp);
    }
    shared_ptr<RDM> tmp = make_shared<RDMI0>(*in, new_index);
    in = tmp;
  }


  list<int> d;
  list<pair<shared_ptr<RDM>, list<int>> > buf(1, make_pair(in,d));

  while (buf.size() != 0) {
    list<pair<shared_ptr<RDM>, list<int>> > buf2;

    for (auto& it : buf) {
      shared_ptr<RDM> tmp = it.first;
      list<int> done = it.second;
      // taking delta
      list<shared_ptr<RDM>> out = tmp->reduce_one(done);
      // this is also needed!
      out.push_back(tmp);
      for (auto& i : out) {
        if (i->reduce_done(done)) {
          rdm_.push_back(i);
        } else {
          buf2.push_back(make_pair(i,done));
        }
      }
    }
    buf = buf2;
  }

  for (auto& i : rdm_)
    i->sort();
}


void Active::print(const string& indent) const {
  for (auto& i : rdm_) i->print(indent);
}


const list<shared_ptr<const Index>> Active::index() const {
  // first find RDM object that does not have any delta functions
  bool done = false;
  auto j = rdm_.begin();
  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) {
    if (!done && (*i)->delta().size() == 0) {
      done = true;
      j = i;
#if 0
    } else if ((*i)->delta().size() == 0) {
      throw logic_error("I think this won't happen. Active::index()");
#endif
    }
  }
  return (*j)->index();
}


bool Active::operator==(const Active& o) const {
  bool out = true;
  // TODO need to add differing order allowance to make code more general
  // compare all rdms of active objects sequentially
  out &= rdm_.size() == o.rdm_.size();
  if (out) {
    for (auto i = rdm_.begin(), j = o.rdm_.begin(); i != rdm_.end(); ++i, ++j)
      out &= (**i) == (**j);
  }
  return out;
}


string Active::generate(const string indent, const string tag, const list<shared_ptr<const Index>> index, const list<shared_ptr<const Index>> merged, const string mlab, const bool use_blas) const {
  stringstream tt;

  vector<string> in_tensors;
  vector<int> req_rdm = required_rdm();
  for (auto& i : req_rdm) {
    stringstream ss; ss << "rdm" << i;
    in_tensors.push_back(ss.str());
  }
  if (!merged.empty()) {
    in_tensors.push_back(mlab);
  }

  for (auto& i : rdm_)
    tt << i->generate(indent, tag, index, merged, mlab, in_tensors, use_blas);
  return tt.str();
}


vector<int> Active::required_rdm() const {
  vector<int> out;
  if (bra_ || ket_) {
    for (auto& i : rdm_) {
      // rdm0 needs to be included as an additional tensor is now needed, <I|0>
      if (find(out.begin(), out.end(), i->rank()) == out.end())
        out.push_back(i->rank());
    }
  } else {
    for (auto& i : rdm_) {
      // rdm0 does not need to be included in header
      if (i->rank() > 0 && find(out.begin(), out.end(), i->rank()) == out.end())
        out.push_back(i->rank());
    }
  }
  sort(out.begin(), out.end());
  return out;
}


void Active::merge(shared_ptr<const Active> o, const double fac) {
  // first merge
  for (auto& i : o->rdm_) {
    shared_ptr<RDM> tmp = i->copy();
    tmp->fac() *= fac;
    rdm_.push_back(tmp);
  }
  // then simplify
  list<list<shared_ptr<RDM>>::iterator> rm;
  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) {
    auto j = i; ++j;
    for ( ; j != rdm_.end(); ++j) {
      if ((*i)->identical(*j)) {
        (*i)->fac() += (*j)->fac();
        rm.push_back(j);
      }
    }
  }
}
