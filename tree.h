//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tree.h
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


#ifndef __TREE_H
#define __TREE_H

#include <memory>
#include "equation.h"
#include "listtensor.h"

class Tree;

class BinaryContraction {
  protected:
    // target_ = tensor_ * subtrees
    std::shared_ptr<Tensor> target_;
    std::shared_ptr<Tensor> tensor_;
    std::list<std::shared_ptr<Tree> > subtree_;

    // Tree that has this
    // TODO in principle, I can write a very general iterator by using this.
    // could not used shared_ptr since it creates cyclic reference 
    Tree* parent_;

  public:
    BinaryContraction(std::shared_ptr<Tensor> o, std::shared_ptr<ListTensor> l);
    BinaryContraction(std::list<std::shared_ptr<Tree> > o, std::shared_ptr<Tensor> t) : subtree_(o), tensor_(t) {};
    ~BinaryContraction() {};

    std::list<std::shared_ptr<Tree> >& subtree() { return subtree_; };
    std::shared_ptr<Tensor> tensor() { return tensor_; };
    std::shared_ptr<Tensor> target() { return target_; };
    void print() const;

    void factorize();
    void set_parent(Tree* o) { parent_ = o; };
    void set_parent_subtree();
    void set_target_rec();
    void set_target(std::shared_ptr<Tensor> o) { target_ = o; };

    bool dagger() const;
    Tree* parent() { return parent_; };
    const Tree* parent() const { return parent_; };

    int depth() const;

    std::vector<std::shared_ptr<Tensor> > tensors_str();
    std::pair<std::string, std::string> generate_task_list() const;
    std::string generate_depend() const;

};


class Tree : public std::enable_shared_from_this<Tree> {
  protected:
    // recursive structure. BinaryContraction contains Tree.
    std::list<std::shared_ptr<BinaryContraction> > bc_;
    // note that target_ can be NULL (at the very beginning) 
    std::shared_ptr<Tensor> target_;

    std::list<std::shared_ptr<Tensor> > op_;

    bool dagger_;

    // TODO raw pointer should be replaced by some better mean
    BinaryContraction* parent_;

    // std::string
    std::string tree_name_;

    // when we do generate, a counter is used to generate a list of tasks.
    mutable int num_;

  public:
    Tree(const std::shared_ptr<Equation> eq);
    Tree(const std::shared_ptr<ListTensor> l);
    ~Tree() {};

    bool done() const;
    void print() const;
    int num() const {
      if (num_ < 0) throw std::logic_error("it seems that the logic is broken - Tree::num_ is not initialized.");
      return num_;
    };

    std::shared_ptr<Tensor> target() const { return target_; };

    std::string generate() const;

    std::string tree_name() const { return tree_name_; };

    void set_parent(BinaryContraction* o) { parent_ = o; };
    void set_parent_sub();
    bool dagger() const { return dagger_; };

    // factorize (i.e., perform factorize in BinaryContraction)
    void factorize();

    // tidy up target tensors. (for some reasons I could not do it on the fly... lazy)
    void set_target_rec();
    void set_target(std::shared_ptr<Tensor> o) {
      target_ = o;
      for (auto i = bc_.begin(); i != bc_.end(); ++i) (*i)->set_target(o);
    };

    bool merge(std::shared_ptr<Tree> o);

    int depth() const;

    // code generators!
    std::pair<std::string,std::string> generate_task_list() const;

};


#endif