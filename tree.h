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

#include <string>
#include <memory>
#include "equation.h"
#include "listtensor.h"

namespace smith {

class Tree;

/// Class framework for tensor multiplication and factorization, contains tree.
class BinaryContraction {
  protected:
    /// An intermediate, target_ = tensor_ * subtrees
    std::shared_ptr<Tensor> target_;
    /// Tensor to be defined (supplied) by BAGEL.
    std::shared_ptr<Tensor> tensor_;
    /// A list of trees
    std::list<std::shared_ptr<Tree>> subtree_;

    // Tree that has this
    // TODO in principle, I can write a very general iterator by using this.
    // could not used shared_ptr since it creates cyclic reference
    /// Points up in graph.
    Tree* parent_;

  public:
    /// Construct binary contraction from tensor and listtensor pointers.
    BinaryContraction(std::shared_ptr<Tensor> o, std::shared_ptr<ListTensor> l);
    /// Construct binary contraction from  subtree and tensor.
    BinaryContraction(std::list<std::shared_ptr<Tree>> o, std::shared_ptr<Tensor> t) : tensor_(t), subtree_(o) {};
    ~BinaryContraction() {};

    /// Return list of trees below.
    std::list<std::shared_ptr<Tree>>& subtree() { return subtree_; };
    /// Return current tensor. For example f1, h1, v2, t2, etc.
    std::shared_ptr<Tensor> tensor() { return tensor_; };
    /// Return target tensor. An intermediate tensor, for example I0.
    std::shared_ptr<Tensor> target() { return target_; };
    /// Retrieve next target--next intermediate below, for example I1. This is the front of subtree of target.
    std::shared_ptr<Tensor> next_target();
    /// Print binary contraction.
    void print() const;

    /// Do factorization and then merge subtrees.
    void factorize();
    /// Set parent to this tree pointer.
    void set_parent(Tree* o) { parent_ = o; };
    /// Sets parent and does set_parent_sub for subtree (bc and tree alternate in graph).
    void set_parent_subtree();
    /// Does set_target_rec for subtrees.
    void set_target_rec();
    /// Set target to this tensor pointer.
    void set_target(std::shared_ptr<Tensor> o) { target_ = o; };

    /// Returns a list of inner loop indices, i.e., those which are the same between tensor_ and target_.
    std::list<std::shared_ptr<Index>> loop_indices();
    /// Returns a list of target indices.
    std::list<std::shared_ptr<Index>> target_indices();

    /// If transpose.
    bool dagger() const;
    /// Returns node above.
    Tree* parent() { return parent_; };
    /// Returns const node above.
    const Tree* parent() const { return parent_; };

    /// Returns the ranks of RDMs in subtree_.
    std::vector<int> required_rdm(std::vector<int> input) const;

    /// Returns depth in graph.
    int depth() const;

    /// Returns vector of tensor with target tensor. 
    std::vector<std::shared_ptr<Tensor>> tensors_str();
    /// Calls generate_task_list for subtree.
    std::pair<std::string, std::string> generate_task_list(const bool enlist = false) const;
    /// TODO this should probably be removed, double check.
    std::string generate_depend() const;

};

/// Used in conjunction with BinaryContraction for overall graph structure. Maps derived equation into tensor contraction tasks.  Task and dependency generation also initiated here.
class Tree : public std::enable_shared_from_this<Tree> {
  protected:
    /// Recursive structure. BinaryContraction contains Tree. A pair of tensors.
    std::list<std::shared_ptr<BinaryContraction>> bc_;
    /// Pointer to target_ tensor. Note that target_ can be NULL (at the very beginning)
    std::shared_ptr<Tensor> target_;
  
    /// Collection of operator tensors.
    std::vector<std::shared_ptr<Tensor>> op_;

    /// If transpose.
    bool dagger_;

    /// Points above. TODO raw pointer should be replaced by some better mean
    BinaryContraction* parent_;

    /// std::string tree name used in generated code.
    std::string tree_name_;

    /// This is only used in the top level for unique Gamma tensors.
    std::list<std::shared_ptr<Tensor>> gamma_;

    /// Add dependency tasks.
    std::string add_depend(const std::shared_ptr<const Tensor> o) const;

    /// When we generate, a counter is used to generate a list of tasks.
    mutable int num_;

    /// Tree label, used for graph-specific code generation. 
    std::string label_;


  public:
    /// Construct tree of equation pointers and set tree label.
    Tree(const std::shared_ptr<Equation> eq, std::string lab="");
    /// Construct tree of listtensor pointers.
    Tree(const std::shared_ptr<ListTensor> l);
    ~Tree() {};

    /// Return label of tree.
    std::string label() { return label_; };
  

    /// TODO check, this seems to be not needed.
    bool done() const;

    /// Prints tree which is the list of operator tensors (op_) if target_ otherwise the binary contraction list (bc_).
    void print() const;
    /// Returns num_. Should be greater than zero, otherwise throws an error.
    int num() const {
      if (num_ < 0) throw std::logic_error("it seems that the logic is broken - Tree::num_ is not initialized.");
      return num_;
    };

    /// Return pointer to target tensor.
    std::shared_ptr<Tensor> target() const { return target_; };
   
    /// TODO check, this seems to be not needed.
    std::string generate() const;

    /// Returns tree name, for generated code..
    std::string tree_name() const { return tree_name_; };

    void set_parent(BinaryContraction* o) { parent_ = o; };
    /// Sets parent and does set_parent_subtree for bc_ (bc_ and tree alternate in graph).
    void set_parent_sub();
    /// If transpose.
    bool dagger() const { return dagger_; };

    /// Factorize (i.e., perform factorize in BinaryContraction).
    void factorize();

    /// Tidy up target tensors. Could be done on the fly. sets for all in bc_ list.
    void set_target_rec();
    /// Set target to this tensor pointer.
    void set_target(std::shared_ptr<Tensor> o) {
      target_ = o;
      for (auto i = bc_.begin(); i != bc_.end(); ++i) (*i)->set_target(o);
    };

    /// Combine trees if tensors are equal, or if have operator tensors.
    bool merge(std::shared_ptr<Tree> o);
    /// Function runs gather_gamma and find_gamma, called from top level (main.cc).
    void sort_gamma(std::list<std::shared_ptr<Tensor>> o = std::list<std::shared_ptr<Tensor>>());
    /// Updates gamma_ with a list of unique Gamma tensors.
    void find_gamma(std::shared_ptr<Tensor> o);

    /// Recursive function to collect all Gamma tensors in graph.
    std::list<std::shared_ptr<Tensor>> gather_gamma() const;
    /// Returns gamma_, list of unique Gamma tensors.
    std::list<std::shared_ptr<Tensor>> gamma() const { return gamma_; };

    /// Returns depth, 0 is top of graph.
    int depth() const;

    /// This function returns the rank of required RDMs here + inp. Goes through bc_ and op_ tensor lists.
    std::vector<int> required_rdm(std::vector<int> inp) const;

    // code generators!
    /// Generate task and task list files.
    std::pair<std::string,std::string> generate_task_list(const bool enlist = false,
        const std::list<std::shared_ptr<Tree>> energy = std::list<std::shared_ptr<Tree>>()) const;

    /// Generate task header.
    std::string generate_compute_header(const int, const std::list<std::shared_ptr<Index>> ti, const std::vector<std::shared_ptr<Tensor>>, const bool, const bool = false) const;
    /// Generate task header.
    std::string generate_compute_footer(const int, const std::list<std::shared_ptr<Index>> ti, const std::vector<std::shared_ptr<Tensor>>, const bool) const;
    /// Generate task for operator task (ie not a binary contraction task).
    std::string generate_compute_operators(const std::string, const std::shared_ptr<Tensor>, const std::vector<std::shared_ptr<Tensor>>,
                                           const bool dagger = false) const;
    /// Generate a task. Here ip is the tag of parent, ic is the tag of this.
    std::string generate_task(const std::string, const int ip, const int ic, const std::vector<std::string>, const bool enlist, const std::string scalar = "") const;
    std::string generate_task(const std::string, const int, const std::vector<std::shared_ptr<Tensor>>, const bool enlist) const;

    /// TODO this should probably be removed, double check.
    std::shared_ptr<Index> generate_rdms() const;

};

}


#endif
