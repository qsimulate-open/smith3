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
