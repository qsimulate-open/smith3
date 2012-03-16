//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tensor.cc
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


#include "tensor.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

Tensor::Tensor(const shared_ptr<Op> op) : factor_(1.0) {
  // label
  label_ = op->label();
  // op
  for (auto i = op->op().begin(); i != op->op().end(); ++i) {
    shared_ptr<Index> in = *get<0>(*i);
    index_.push_back(in);
  }

} 


Tensor::Tensor(const shared_ptr<Active> activ) : factor_(1.0) {
  // label
  label_ = "Gamma";
  // op
  index_ = activ->index();
  active_ = activ;
} 


std::string Tensor::str() const {
  stringstream ss;
  if (factor_ != 1.0) ss << " " << fixed << setw(4) << setprecision(2) << factor_ << " "; 
  ss << label_ << "(";
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    // we don't need the spin part here
    if (i != index_.begin()) ss << ", ";
    ss << (*i)->str(false);
  }
  ss << ")";

  if (merged_) ss << " << " << merged_->str();
  return ss.str();
}


// returns if all the indices are of active orbitals
bool Tensor::all_active() const {
  bool out = true;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    out &= (*i)->active(); 
  }
  return out;
}


// adds all-active tensor to Active_;
void Tensor::merge(shared_ptr<Tensor> a) {
  assert(active_);  
  merged_ = a; 
  list<list<shared_ptr<Index> >::iterator> remove;
  // remove Index that belongs to a
  for (auto i = a->index().begin(); i != a->index().end(); ++i) {
    const int n = (*i)->num();
    for (auto j = index_.begin(); j != index_.end(); ++j) {
      // careful, this code is driven by numbers
      if ((*j)->num() == n) {
        remove.push_back(j);
        break;
      }
    }
  }
  for (auto i = remove.begin(); i != remove.end(); ++i) {
    index_.erase(*i);
  }
}


bool Tensor::operator==(const Tensor& o) const {
  bool out = true;
  out &= label_ == o.label();
  out &= factor_ == o.factor();
  // TODO for the time being, active is not assumed to be factorized
  if (active() || o.active()) out = false;
  // index
  out &= index_.size() == o.index().size();
  if (out) {
    for (auto i = index_.begin(), j = o.index().begin(); i != index_.end(); ++i, ++j) {
      out &= (*i)->identical(*j);
    } 
  }
  return out;
}


string Tensor::constructor_str(std::string indent) const {
  stringstream ss;
  ss << indent << "std::vector<IndexRange> " << label_ << "_index";
  if (index_.empty()) {
    ss << ";" << endl;
  } else {
    ss << " = vec(";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << (i != index_.rbegin() ? ", this->" : "this->") << (*i)->generate();
    ss << ");" << endl;
  }
  ss << indent << "std::shared_ptr<Tensor<T> > " << label_ << "(new Tensor<T>(" << label_ << "_index, false));";
  return ss.str();
}

string Tensor::generate_get_block(const string cindent, const string lab, const bool move) const {
  stringstream tt;
  tt << cindent << "std::vector<size_t> " << lab << "hash = vec(";
  for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
    if (iter != index_.rbegin()) tt << ", ";
    tt << (*iter)->str_gen() << "->key()";
  }   
  tt << ");" << endl; 
  {   
    tt << cindent << "std::unique_ptr<double[]> " << lab << "data = "
                  << (label_ == "proj" ? "r" : label_) << "->" << (move ? "move" : "get") << "_block(" << lab << "hash);" << endl;
  }   
  return tt.str();
}


// TODO replace by a standard function (since I am aboard, I cannot google..)
static double abs__(const double& a) { return a > 0 ? a : -a; };
static string prefac__(const double& factor_) {
  stringstream ss;
  if (abs__(factor_-1.0) < 1.0e-10) { ss << "1,1";
  } else if (abs__(factor_+1.0) < 1.0e-10) { ss << "-1,1";
  } else if (abs__(factor_-2.0) < 1.0e-10) { ss << "2,1";
  } else if (abs__(factor_+2.0) < 1.0e-10) { ss << "-2,1";
  } else if (abs__(factor_-4.0) < 1.0e-10) { ss << "4,1";
  } else if (abs__(factor_+4.0) < 1.0e-10) { ss << "-4,1";
  } else if (abs__(factor_-8.0) < 1.0e-10) { ss << "8,1";
  } else if (abs__(factor_+8.0) < 1.0e-10) { ss << "-8,1";
  } else if (abs__(factor_-0.5) < 1.0e-10) { ss << "1,2";
  } else if (abs__(factor_+0.5) < 1.0e-10) { ss << "-1,2";
  } else if (abs__(factor_-0.25) < 1.0e-10) { ss << "1,4";
  } else if (abs__(factor_+0.25) < 1.0e-10) { ss << "-1,4";
  } else {
    ss << "this case is not yet considered " << factor_ << " in Tensor::generate_sort_indices()";
    throw runtime_error(ss.str());
  }
  return ss.str();
}


string Tensor::generate_scratch_area(const string cindent, const string lab, const bool zero) const {
  stringstream ss;
  ss << cindent << "std::unique_ptr<double[]> " << lab << "data_sorted(new double["
                << (label_ == "proj" ? "r" : label_) << "->get_size(" << lab << "hash)]);" << endl;
  if (zero) {
    ss << cindent << "std::fill(" << lab << "data_sorted.get(), " << lab << "data_sorted.get()+"
                  << (label_ == "proj" ? "r" : label_) << "->get_size(" << lab << "hash), 0.0);" << endl;
  }
  return ss.str();
}

string Tensor::generate_sort_indices(const string cindent, const string lab, const list<shared_ptr<Index> >& loop, const bool op) const {
  stringstream ss;
  if (!op) {
    ss << generate_scratch_area(cindent, lab);
  }
  vector<int> map(index_.size());
  ss << cindent << "sort_indices<";
  // determine mapping
  // first loop indices. order as in loop
  vector<int> done;
  for (auto i = loop.rbegin(); i != loop.rend(); ++i) {
    // count
    int cnt = 0;
    for (auto j = index_.rbegin(); j != index_.rend(); ++j, ++cnt) {
      if ((*i)->identical(*j)) break;
    }
    if (cnt == index_.size()) throw logic_error("should not happen.. Tensor::generate_sort_indices");
    ss << cnt << ",";
    done.push_back(cnt);
  }
  // then fill out others
  for (int i = 0; i != index_.size(); ++i) {
    if (find(done.begin(), done.end(), i) == done.end())
      ss << i << ",";
  }

  string target_label = op ? "odata" : lab + "data_sorted";

  ss << "1,1," << prefac__(factor_);
  ss << ">(" << lab << "data, " << target_label; 
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) ss << ", " << (*i)->str_gen() << "->size()";
  ss << ");" << endl;
  return ss.str();
}


string Tensor::generate_sort_indices_target(const string cindent, const string lab, const list<shared_ptr<Index> >& loop,
                                            const shared_ptr<Tensor> a, const shared_ptr<Tensor> b) const {
  stringstream ss;
  vector<int> map(index_.size());
  ss << cindent << "sort_indices<";

  // determine mapping
  // first obtain the ordering of indices from dgemm
  list<shared_ptr<Index> > source;
  {
    list<shared_ptr<Index> > aind = a->index();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto j = loop.begin(); j != loop.end(); ++j)
        if ((*i)->identical(*j)) found = true;
      if (!found) source.push_back(*i);
    }
    aind = b->index();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto j = loop.begin(); j != loop.end(); ++j)
        if ((*i)->identical(*j)) found = true;
      if (!found) source.push_back(*i);
    }
  }
  
  // note that source is already reversed.
  for (auto i = source.begin(); i != source.end(); ++i) {
    // count
    int cnt = 0;
    // but index_ is not; that's why the reverse_iterator is used here
    for (auto j = index_.rbegin(); j != index_.rend(); ++j, ++cnt) {
      if ((*i)->identical(*j)) break;
    }
    if (cnt == index_.size()) throw logic_error("should not happen.. Tensor::generate_sort_indices_target");
    ss << cnt << ",";
  }

  ss << "1,1," << prefac__(factor_);
  ss << ">(" << lab << "data_sorted, " << lab << "data"; 
  for (auto i = source.begin(); i != source.end(); ++i) ss << ", " << (*i)->str_gen() << "->size()";
  ss << ");" << endl;
  return ss.str();
}


pair<string, string> Tensor::generate_dim(const list<shared_ptr<Index> >& di) const {
  vector<string> s, t;
  // first indices which are not shared
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    bool shared = false;
    for (auto j = di.begin(); j != di.end(); ++j) {
      if ((*i)->identical(*j)) {
        shared = true;
        break;
      }
    }
    if (shared) {
      t.push_back((*i)->str_gen() + "->size()");
    } else {
      s.push_back((*i)->str_gen() + "->size()");
    }
  }

  stringstream ss, tt;
  for (auto i = s.begin(); i != s.end(); ++i) {
    if (i != s.begin()) ss << "*";
    ss << *i;
  }
  for (auto i = t.begin(); i != t.end(); ++i) {
    if (i != t.begin()) tt << "*";
    tt << *i;
  }
  return make_pair(ss.str(), tt.str());
}
