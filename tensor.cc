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
#include "constants.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace smith;

Tensor::Tensor(const shared_ptr<Op> op) : factor_(1.0) {
  // label
  label_ = op->label();
  // op
  for (auto& i : op->op()) {
    shared_ptr<Index> in = *get<0>(i);
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


string Tensor::str() const {
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
  for (auto& i : a->index()) {
    const int n = i->num();
    for (auto j = index_.begin(); j != index_.end(); ++j) {
      // careful, this code is driven by numbers
      if ((*j)->num() == n) {
        remove.push_back(j);
        break;
      }
    }
  }
  for (auto& i : remove) index_.erase(i);
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


string Tensor::constructor_str(string indent) const {
  stringstream ss;
  ss << indent << "std::vector<IndexRange> " << label_ << "_index";
  if (index_.empty()) {
    ss << ";" << endl;
  } else {
    // mkm vec() not working in bagel
    //ss << " = vec(";
    ss << " = {";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << (i != index_.rbegin() ? ", this->" : "this->") << (*i)->generate();
    //ss << ");" << endl;
    ss << "};" << endl;
  }
  ss << indent << "std::shared_ptr<Tensor<T> > " << label_ << "(new Tensor<T>(" << label_ << "_index, false));";
  return ss.str();
}

string Tensor::generate_get_block(const string cindent, const string lab, const bool move) const {
  string lbl = label_;
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  bool trans = false;
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
    trans = true;
  }

  stringstream tt;
  tt << cindent << "std::vector<size_t> " << lab << "hash = {";
  if (!trans) {
    tt << list_keys(index_);
  } else {
    assert(!(index_.size() & 1));
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
      string i0 = (*iter)->str_gen() + ".key()";
      ++iter;
      tt << (*iter)->str_gen() << ".key()" << ", " << i0;
    }
  }
  tt << "};" << endl;
  {
    tt << cindent << "std::unique_ptr<double[]> " << lab << "data = "
                  << lbl << "->" << (move ? "move" : "get") << "_block(" << lab << "hash);" << endl;
  }
  return tt.str();
}


string Tensor::generate_scratch_area(const string cindent, const string lab, const bool zero) const {
  string lbl = label_;
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
  }

  stringstream ss;
  ss << cindent << "std::unique_ptr<double[]> " << lab << "data_sorted(new double["
                << lbl << "->get_size(" << lab << "hash)]);" << endl;
  if (zero) {
    ss << cindent << "std::fill(" << lab << "data_sorted.get(), " << lab << "data_sorted.get()+"
                  << lbl << "->get_size(" << lab << "hash), 0.0);" << endl;
  }
  return ss.str();
}

string Tensor::generate_sort_indices(const string cindent, const string lab, const list<shared_ptr<Index> >& loop, const bool op) const {
  stringstream ss;
  if (!op) ss << generate_scratch_area(cindent, lab);

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
    done.push_back(cnt);
  }
  // then fill out others
  for (int i = 0; i != index_.size(); ++i) {
    if (find(done.begin(), done.end(), i) == done.end())
      done.push_back(i);
  }
  // if trans, transpose here!
  size_t found = label_.find("dagger");
  const bool trans = found != string::npos;
  if (trans && index_.size() & 1) throw logic_error("transposition not possible with 3-index objects");
  if (trans) {
    vector<int> tmp;
    for (int i = 0; i != done.size(); i += 2) {
      tmp.push_back(done[i+1]);
      tmp.push_back(done[i]);
    }
    done = tmp;
  }

  // then write them out.
  for (auto i = done.begin(); i != done.end(); ++i)
    ss << *i << ",";

  string target_label = op ? "odata" : lab + "data_sorted";

  ss << (op ? 1 : 0) << ",1," << prefac__(factor_);
  ss << ">(" << lab << "data, " << target_label;
  if (!trans) {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << ", " << (*i)->str_gen() << ".size()";
  } else {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
      string tmp = ", " + (*i)->str_gen() + ".size()";
      ++i;
      ss << ", " << (*i)->str_gen() << ".size()" << tmp;
    }
  }
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
      for (auto& j : loop)
        if ((*i)->identical(j)) found = true;
      if (!found) source.push_back(*i);
    }
    aind = b->index();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto& j : loop)
        if ((*i)->identical(j)) found = true;
      if (!found) source.push_back(*i);
    }
  }

  for (auto j = index_.rbegin(); j != index_.rend(); ++j) {
    // count
    int cnt = 0;
    // note that source is already reversed!
    for (auto i = source.begin(); i != source.end(); ++i, ++cnt) {
      if ((*i)->identical(*j)) break;
    }
    if (cnt == index_.size()) throw logic_error("should not happen.. Tensor::generate_sort_indices_target");
    ss << cnt << ",";
  }

  ss << "1,1," << prefac__(factor_);
  ss << ">(" << lab << "data_sorted, " << lab << "data";
  for (auto i = source.begin(); i != source.end(); ++i) ss << ", " << (*i)->str_gen() << ".size()";
  ss << ");" << endl;
  return ss.str();
}


pair<string, string> Tensor::generate_dim(const list<shared_ptr<Index> >& di) const {
  vector<string> s, t;
  // first indices which are not shared
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    bool shared = false;
    for (auto& j : di) {
      if ((*i)->identical(j)) {
        shared = true;
        break;
      }
    }
    if (shared) {
      t.push_back((*i)->str_gen() + ".size()");
    } else {
      s.push_back((*i)->str_gen() + ".size()");
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


string Tensor::generate_active(string indent, const string tag) const {
  assert(label() == "Gamma");
  stringstream tt;
  if (!merged_) {
    tt << active()->generate(indent, tag, index());
  } else {
    // add merge loops
    tt << indent <<"// associated with merged" << endl;

    vector<string> mclose;
    list<shared_ptr<Index> >& merged = merged_->index();
    // the m->generate() makes active_ this label must be identical on bagel side (src/smith/spinfreebase.h)
    for (auto& m : merged) {
      tt << indent << "for (auto& " << m->str_gen() << " : "  <<  m->generate()  << ") {" << endl;
      mclose.push_back(indent + "}");
      indent += "  ";
    }
    // make get block for merge obj
    tt << indent << "std::vector<size_t> fhash = {" << list_keys(merged) << "};" << endl;
    tt << indent << "std::unique_ptr<double[]> fdata = " << merged_->label() << "->get_block(fhash);" << endl;

    tt << active()->generate_merged(indent, tag, index(), merged_->index(), merged_->label());

    // close merge for loops
    for (auto iter = mclose.rbegin(); iter != mclose.rend(); ++iter)
      tt << *iter << endl;

  }
  return tt.str();
}


string Tensor::generate_loop(string& indent, vector<string>& close) const {
  stringstream tt;
  for (auto iter = index_.begin(); iter != index_.end(); ++iter, indent += "  ") {
    string cindex = (*iter)->str_gen();
    tt << indent << "for (auto& " << cindex << " : " << (*iter)->generate() << ") {" << endl;
    close.push_back(indent + "}");
  }
  return tt.str();
}


string Tensor::generate_gamma(string& indent, vector<string>& close, string tag, const bool move) const {
  assert(label() == "Gamma");
  stringstream tt;
  tt << generate_loop(indent, close);
  // generate gamma get block, true does a move_block
  tt << generate_get_block(indent, tag, move);
  // now generate codes for rdm
  tt << generate_active(indent, tag);
  return tt.str();
}
