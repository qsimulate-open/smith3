//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tree.cc
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


#include "tree.h"
#include "constants.h"

using namespace std;
using namespace smith;


Tree::Tree(shared_ptr<Equation> eq, string lab) : parent_(NULL), tree_name_(eq->name()), num_(-1), label_(lab), root_targets_(eq->targets()) {
  // First make ListTensor for all the diagrams
  list<shared_ptr<Diagram>> d = eq->diagram();

  const bool rt_targets = eq->targets();

  for (auto& i : d) {
    shared_ptr<ListTensor> tmp = make_shared<ListTensor>(i);
    // All internal tensor should be included in the active part
    tmp->absorb_all_internal();

    // rearrange brakets and reindex associated tensors, ok if not complex
    tmp->absorb_ket();

    shared_ptr<Tensor> first = tmp->front();
    shared_ptr<ListTensor> rest = tmp->rest();

    // reorder to minimize the cost
    rest->reorder();

    // convert to tree and then bc
    list<shared_ptr<Tree>> lt { make_shared<Tree>(rest, lab, rt_targets) };
    shared_ptr<BinaryContraction> b = make_shared<BinaryContraction>(lt, first, i->target_index());
    bc_.push_back(b);
  }

  factorize();
  move_up_operator();
  set_parent_sub();
  set_target_rec();
}


Tree::Tree(const shared_ptr<ListTensor> l, string lab, const bool t) : num_(-1), label_(lab), root_targets_(t) {
  target_ = l->target();
  dagger_ = l->dagger();
  if (l->length() > 1) {
    shared_ptr<BinaryContraction> bc = make_shared<BinaryContraction>(target_, l, lab, root_targets_);
    bc_.push_back(bc);
  } else {
    shared_ptr<Tensor> t = l->front();
    t->set_factor(l->fac());
    t->set_scalar(l->scalar());
    op_.push_back(t);
  }
}


BinaryContraction::BinaryContraction(shared_ptr<Tensor> o, shared_ptr<ListTensor> l, string lab, const bool rt) : label_(lab) {
  // o is a target tensor NOT included in l
  target_ = o;
  tensor_ = l->front();
  shared_ptr<ListTensor> rest = l->rest();
  subtree_.push_back(make_shared<Tree>(rest, label_, rt));
}


void Tree::factorize() {
  list<list<shared_ptr<BinaryContraction>>::iterator> done;
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    auto j = i; ++j;
    for ( ; j != bc_.end(); ++j) {
      // bc is factorized according to excitation target indices
      if (*(*i)->tensor() == *(*j)->tensor() && (*i)->dagger() == (*j)->dagger() && (*i)->target_index_str() == (*j)->target_index_str()) {
        done.push_back(j);
        (*i)->subtree().insert((*i)->subtree().end(), (*j)->subtree().begin(), (*j)->subtree().end());
      }
    }
  }
  for (auto& i : done) bc_.erase(i);
  for (auto& i : bc_) i->factorize();
}


void BinaryContraction::factorize() {
  list<list<shared_ptr<Tree>>::iterator> done;
  // and then merge subtrees
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    auto j = i; ++j;
    for ( ; j != subtree_.end(); ++j) {
      if (find(done.begin(), done.end(), j) != done.end()) continue;
      if ((*i)->merge(*j)) done.push_back(j);
    }
  }
  for (auto i = done.begin(); i != done.end(); ++i) subtree_.erase(*i);
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) (*i)->factorize();
}


void Tree::move_up_operator() {
  for (auto& b : bc_)
    b->move_up_operator();
}


void BinaryContraction::move_up_operator() {
  if (subtree_.size() == 1 && subtree_.front()->can_move_up()) {
    source_ = subtree_.front()->op().front();
    subtree_.clear();
  }

  for (auto& i : subtree_)
    i->move_up_operator();
}


bool Tree::merge(shared_ptr<Tree> o) {
  bool out = false;
  if (o->bc_.size() > 0) {
    if (bc_.size() == 0 || *bc_.front()->tensor() == *o->bc_.front()->tensor()) {
      out = true;
      bc_.insert(bc_.end(), o->bc_.begin(), o->bc_.end());
      for (auto& i : bc_) i->set_target(target_);
    }
  } else if (o->op_.size() > 0) {
    out = true;
    op_.insert(op_.end(), o->op_.begin(), o->op_.end());
  }
  return out;
}


void BinaryContraction::set_target_rec() {
  if (!subtree_.empty()) {
    auto i = subtree_.begin();
    shared_ptr<Tensor> first = (*i)->target();
    for (++i ; i != subtree_.end(); ++i) {
      (*i)->set_target(first);
    }
    for (i = subtree_.begin() ; i != subtree_.end(); ++i)
      (*i)->set_target_rec();
  }
}


void Tree::set_target_rec() {
  if (!bc_.empty()) {
    for (auto& i : bc_) i->set_target_rec();
  }
}


void BinaryContraction::set_parent_subtree() {
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    (*i)->set_parent(this);
    (*i)->set_parent_sub();
  }
}


void Tree::set_parent_sub() {
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    (*i)->set_parent(this);
    (*i)->set_parent_subtree();
  }
}


int BinaryContraction::depth() const { return parent_->depth(); }

int Tree::depth() const { return parent_ ? parent_->parent()->depth()+1 : 0; }


bool BinaryContraction::dagger() const {
  return subtree_.front()->dagger();
}


shared_ptr<Tensor> BinaryContraction::next_target() {
  return !subtree_.empty() ? subtree().front()->target() : source_;
}


void Tree::sort_gamma(list<shared_ptr<Tensor>> o) {
  gamma_ = o;
  list<shared_ptr<Tensor>> g = gather_gamma();
  for (auto& i : g) find_gamma(i);
}


void Tree::find_gamma(shared_ptr<Tensor> o) {
  bool found = false;
  for (auto& i : gamma_) {
    if ((*i) == (*o)) {
      found = true;
      o->set_alias(i);
      break;
    }
  }
  if (!found) gamma_.push_back(o);
}


// return gamma below this
list<shared_ptr<Tensor>> Tree::gather_gamma() const {
  list<shared_ptr<Tensor>> out;
  for (auto& i : bc_) {
    for (auto& j : i->subtree()) {
      // recursive call
      list<shared_ptr<Tensor>> tmp = j->gather_gamma();
      out.insert(out.end(), tmp.begin(), tmp.end());
    }
    if (i->source() && i->source()->label().find("Gamma") != string::npos)
      out.push_back(i->source());
  }
  for (auto& i : op_)
    if (i->label().find("Gamma") != string::npos) out.push_back(i);
  for (auto& b : bc_)
    if (b->tensor()->label().find("Gamma") != string::npos) out.push_back(b->tensor());
  return out;
}


string BinaryContraction::target_index_str() const {
  stringstream zz;
  if (!target_index_.empty()) {
    zz << "(";
    for (auto i = target_index_.begin(); i != target_index_.end(); ++i) {
      if (i != target_index_.begin()) zz << ", ";
      zz << (*i)->str(false);
    }
    zz << ")";
  }
  return zz.str();
}


void BinaryContraction::print() const {
  string indent = "";
  for (int i = 0; i != depth(); ++i) indent += "  ";
  cout << indent << (target_ ? (target_->str() + " = ") : "") << tensor_->str() << (depth() == 0 ? target_index_str() : "") << " * "
                 << (subtree_.empty() ? source_->str() : subtree_.front()->target()->str()) << endl;
  for (auto& i : subtree_) i->print();
}


void Tree::print() const {
  string indent = "";
  for (int i = 0; i != depth(); ++i) indent += "  ";
  if (target_) {
    for (auto i = op_.begin(); i != op_.end(); ++i)
      cout << indent << target_->str() << " += " << (*i)->str() << (dagger_ ? " *" : "") << endl;
  }
  for (auto i = bc_.begin(); i != bc_.end(); ++i)
    (*i)->print();
}


vector<shared_ptr<Tensor>> BinaryContraction::tensors_vec() {
  vector<shared_ptr<Tensor>> out;
  if (target_) out.push_back(target_);
  out.push_back(tensor_);
  if (!subtree_.empty())
    out.push_back(subtree_.front()->target());
  else if (source_)
    out.push_back(source_);

  return out;
}


bool BinaryContraction::diagonal_only() const {
  return tensor_->label().find("Gamma") == string::npos
      && all_of(subtree_.begin(), subtree_.end(), [](shared_ptr<Tree> i){ return i->diagonal_only(); })
      && (!source_ || source_->label().find("Gamma") == string::npos)
      && nogamma_upstream();
}


bool BinaryContraction::nogamma_upstream() const {
  return tensor_->label().find("Gamma") == std::string::npos
      && parent_->nogamma_upstream();
}


vector<int> BinaryContraction::required_rdm(vector<int> orig) const {
  vector<int> out = orig;
  for (auto& i : subtree_)
    out = i->required_rdm(out);

  sort(out.begin(), out.end());
  return out;
}

vector<int> Tree::required_rdm(vector<int> orig) const {
  vector<int> out = orig;
  for (auto& i : bc_)
    out = i->required_rdm(out);

  for (auto& i : op_) {
    if (i->label().find("Gamma") == string::npos) continue;
    vector<int> rdmn = i->active()->required_rdm();
    for (auto& j: rdmn)
      if (find(out.begin(), out.end(), j) == out.end()) out.push_back(j);
  }
  sort(out.begin(), out.end());
  return out;
}



