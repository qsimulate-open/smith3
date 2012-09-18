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
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "active.h"

using namespace std;

void RDM::print(const string& indent) const {

  cout << indent << fixed << setw(5) << setprecision(1) << fac_ << " [";
  for (auto i = index_.begin(); i != index_.end(); ++i)
    cout << (*i)->str();
  for (auto i = delta_.begin(); i != delta_.end(); ++i)
    cout << " d(" << i->first->str(false) << i->second->str(false) << ")";
  cout << "]";

  if (done()) cout << "*";

  cout << endl;
}


shared_ptr<RDM> RDM::copy() const {
  // first clone all the indices
  list<shared_ptr<Index> > in;
  for (auto i = index_.begin(); i != index_.end(); ++i) in.push_back((*i)->clone());

  map<shared_ptr<Spin>, shared_ptr<Spin> > dict;
  auto j = in.begin();
  for (auto i = index_.begin(); i != index_.end(); ++i, ++j) {
    // get the original spin 
    shared_ptr<Spin> o = (*i)->spin();
    if (dict.find(o) == dict.end()) {
      (*j)->set_spin(o);
    } else {
      shared_ptr<Spin> s(new Spin());
      s->set_num(o->num());
      dict.insert(make_pair(o,s));
      (*j)->set_spin(s);
    }
  }

  // lastly clone all the delta functions
  list<pair<shared_ptr<Index>, shared_ptr<Index> > > d;
  for (auto i = delta_.begin(); i != delta_.end(); ++i) d.push_back(make_pair(i->first->clone(), i->second->clone())); 

  shared_ptr<RDM> out(new RDM(in, d)); 
  out->fac() = fac_;
  return out;
}

//
// An application of "Wick's theorem"
// This function is controlled by Index::num_. Not a great code, it could have been driven by pointers... lazy me.
//
list<shared_ptr<RDM> > RDM::reduce_one(list<int>& done) const {
  // first find non-daggered operator which is not aligned
  list<shared_ptr<RDM> > out;

  for (auto i = index_.begin(); i != index_.end(); ++i) {

    // looking for a non-daggered index that is not registered in done
    if ((*i)->dagger() || find(done.begin(), done.end(), (*i)->num()) != done.end())
      continue;

    // again this function is controlled by numbers... sorry... 
    const int inum = (*i)->num();

    for (auto j = i; j != index_.end(); ++j) {
      if (!(*j)->dagger() || j==i) continue;

      // if you find daggered object in the right hand side...
      shared_ptr<RDM> tmp = this->copy(); 

      // find the indices to be deleted.
      vector<list<shared_ptr<Index> >::iterator> rml;
      int cnt0 = -1;
      int cnt = 0;
      for (auto k = tmp->index().begin(); k != tmp->index().end(); ++k, ++cnt) {
        if ((*k)->same_num(*i) || (*k)->same_num(*j)) {
          cnt0 = cnt0 >= 0 ? cnt-cnt0 : cnt;
          rml.push_back(k);
        }
      }
      assert(rml.size() == 2);

      tmp->delta().push_back(make_pair(*rml[0],*rml[1]));
      // Please note that this procedure does not change the sign (you can prove it in 30sec)
      tmp->fac() *= ((cnt0-1)&1 ? -1.0 : 1.0);
      if ((*i)->same_spin(*j)) {
        tmp->fac() *= 2.0;
      } else {
        // this case we need to replace a spin 
        const shared_ptr<Spin> s0 = (*rml[0])->spin();
        const shared_ptr<Spin> s1 = (*rml[1])->spin();
        for (auto k = tmp->index().begin(); k != tmp->index().end(); ++k) {
          if ((*k)->spin() == s0) {
            (*k)->set_spin(s1);
          }
        }
      }

      // erasing indices which are push-backed in delta
      tmp->index().erase(rml[0]);
      tmp->index().erase(rml[1]);
      out.push_back(tmp);
    }
    done.push_back((*i)->num());
    break;
  }
  return out;
}


bool RDM::reduce_done(const list<int>& done) const {
  // check if there is a annihilation operator which has creation operators in his right side
  bool out = true;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    // if non-dagger and not registered in done
    if (!(*i)->dagger() && find(done.begin(), done.end(), (*i)->num()) == done.end()) {
      for (auto j = i; j != index_.end(); ++j) {
        if ((*j)->dagger()) out = false;
      } 
      break;
    }
  }
  return out;
}


void RDM::sort() {

  // of course this is not the fastest code, but I am fine. 

  // first we align indices so that
  // 0+ 0 1+ 1 2+ 2...
  // actually this might be better for actually implementation. 
  vector<shared_ptr<Spin> > done_spin;
  while (!done()) {

    list<shared_ptr<Index> > buf;
    auto i = index_.begin();
    // continue to spin which is not processed
    for (; i != index_.end(); ++i) {
      if (find(done_spin.begin(), done_spin.end(), (*i)->spin()) != done_spin.end()) {
        buf.push_back(*i); 
      } else {
        break;
      }
    }
    {
      shared_ptr<Spin> cs = (*i)->spin();
      const bool dagger = (*i)->dagger();
      auto j = i;
      int cnt = 0;
      bool found = false;

      if (dagger) {
        // if dagger, move it to right before the nondagger of the same spin
        for (++j; j != index_.end(); ++j) {
          if ((*j)->spin() == cs) {
            buf.push_back(*i); 
            buf.push_back(*j); 
            assert(!(*j)->dagger()); 
            found = true;
          } else {
            buf.push_back(*j); 
            if (!found) ++cnt;
          }
        }
      } else {
        // if nodagger, move it to right after the dagger of the same spin
        for (++j; j != index_.end(); ++j) {
          if ((*j)->spin() == cs) {
            buf.push_back(*j);
            assert((*j)->dagger()); 
            buf.push_back(*i); 
            found = true;
          } else {
            buf.push_back(*j); 
            if (!found) ++cnt;
          }
        }
      }
      fac_ *= cnt&1 ? 1 : -1;
      done_spin.push_back(cs);
    }
    if (index_.size() != buf.size()) {
      for (auto z = buf.begin(); z != buf.end(); ++z) (*z)->print();
      throw logic_error("RDM::sort()");
    }
    index_ = buf;
  }
}


bool RDM::done() const {
  // if operators are aligned as a0+ a0 a1+ a1...
  bool out = true;
  assert((index_.size()&1) == 0); // for sure..

  int cnt = 0;
  shared_ptr<Spin> prev;
  for (auto i = index_.begin(); i != index_.end(); ++i, ++cnt) {
    // even number, then (*i) should be daggered.
    if ((cnt & 1) == 0) {
      if (!(*i)->dagger()) {
        out = false;
        break;
      }
      prev = (*i)->spin();
    } else {
      if ((*i)->dagger() || (*i)->spin() != prev) {
        out = false;
        break;
      }
    }
  }
  return out;
}


Active::Active(const list<shared_ptr<Index> >& in) {

  list<pair<shared_ptr<Index>, shared_ptr<Index> > > t;
  shared_ptr<RDM> tmp(new RDM(in, t, 1.0));

  // this sets list<RDM>
  reduce(tmp);

}


void Active::reduce(shared_ptr<RDM> in) {
  list<int> d;
  list<pair<shared_ptr<RDM>, list<int> > > buf(1, make_pair(in,d));

  while (buf.size() != 0) {
    list<pair<shared_ptr<RDM>, list<int> > > buf2;

    for (auto iter = buf.begin(); iter != buf.end(); ++iter) { 
      shared_ptr<RDM> tmp = iter->first;
      list<int> done = iter->second; 
      
      // taking delta
      list<shared_ptr<RDM> > out = tmp->reduce_one(done);
      // this is also needed!
      out.push_back(tmp);
      for (auto i = out.begin(); i != out.end(); ++i) {
        if ((*i)->reduce_done(done)) {
          rdm_.push_back(*i);
        } else {
          buf2.push_back(make_pair(*i,done));
        }
      }
    }
    buf = buf2;
  }

  for (auto i = rdm_.begin(); i != rdm_.end(); ++i)
    (*i)->sort();
}


void Active::print(const string& indent) const {
  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) (*i)->print(indent);
}


const list<shared_ptr<Index> > Active::index() const {
  // first find RDM object that does not have any delta functions
  bool done = false;
  auto j = rdm_.begin();
  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) {
    if (!done && (*i)->delta().size() == 0) {
      done = true;
      j = i; 
    } else if ((*i)->delta().size() == 0) {
      throw logic_error("I think this won't happen. Active::index()");
    }
  }
  return (*j)->index();
}

string Active::generate_get_block(const string cindent, const string lab, const bool move, const string label) const {

  string lbl = label;
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  bool trans = false;
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
    trans = true;
  }

  stringstream tt;
  tt << cindent << "std::vector<size_t> " << lab << "hash = vec(";
  if (!trans) {
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
      tt << (*iter)->str_gen() << "->key()";
    }
  } else {
    assert(!(index_.size() & 1));
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
      string i0 = (*iter)->str_gen() + "->key()";
      ++iter;
      tt << (*iter)->str_gen() << "->key()" << ", " << i0;
    }
  }
  tt << ");" << endl;
  {
    tt << cindent << "std::unique_ptr<double[]> " << lab << "data = "
                  << lbl << "->" << (move ? "move" : "get") << "_block(" << lab << "hash);" << endl;
  }
  return tt.str();
}

