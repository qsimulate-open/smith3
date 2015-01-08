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

#include "equation.h"
#include "listtensor.h"

namespace smith {

class Tree;

/// Class framework for tensor multiplication and factorization, contains tree.
class BinaryContraction {
  protected:
    /// Intermediate on LHS of equation:  target_ = tensor_ * subtrees
    std::shared_ptr<Tensor> target_;
    /// Tensor. Can be defined (supplied) by BAGEL.
    std::shared_ptr<Tensor> tensor_;
    /// A list of trees
    std::list<std::shared_ptr<Tree>> subtree_;

    /// label for code generation specifics, corresponds to derived tree name.
    std::string label_;

    // Tree that has this
    // TODO in principle, I can write a very general iterator by using this.
    // could not used shared_ptr since it creates cyclic reference
    /// Points up in graph.
    Tree* parent_;

    /// Target indices, could be from excitation operator target indices, or ci derivative target index.
    std::list<std::shared_ptr<const Index>> target_index_;

  public:
    /// Construct binary contraction from subtree and tensor if diagram has excitation operator target indices, index list will not be empty.
    BinaryContraction(std::list<std::shared_ptr<Tree>> o, std::shared_ptr<Tensor> t, std::list<std::shared_ptr<const Index>> ti) : tensor_(t), subtree_(o), target_index_(ti) { }
    /// Construct binary contraction from tensor and listtensor pointers.
    BinaryContraction(std::shared_ptr<Tensor> o, std::shared_ptr<ListTensor> l, std::string lab, const bool rt);
    ~BinaryContraction() { }


    /// Return list of trees below.
    std::list<std::shared_ptr<Tree>>& subtree() { return subtree_; }
    /// Return current tensor. For example f1, h1, v2, t2, etc.
    std::shared_ptr<Tensor> tensor() { return tensor_; }
    /// Return target tensor. An intermediate tensor, for example I0.
    std::shared_ptr<Tensor> target() { return target_; }
    /// Retrieve next target--next intermediate below, for example I1. This is the front of subtree of target.
    std::shared_ptr<Tensor> next_target();
    /// Returns vector of tensor with target tensor.
    std::vector<std::shared_ptr<Tensor>> tensors_vec();

    /// Print binary contraction.
    void print() const;
    /// Returns the target indices.
    std::string target_index_str() const;
    /// Returns bc label.
    std::string label() { return label_; }

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
    std::list<std::shared_ptr<const Index>> loop_indices();
    /// Returns a list of target indices..these are to be stored (via put_block).
    std::list<std::shared_ptr<const Index>> target_indices();
    /// Return excitation target indices.
    std::list<std::shared_ptr<const Index>> target_index() { return target_index_; }
    /// True if bc is associated with excitation target indices.
    bool has_target_index() { return !target_index_.empty(); }

    /// If transpose.
    bool dagger() const;

    /// Returns node above.
    Tree* parent() { return parent_; }
    /// Returns const node above.
    const Tree* parent() const { return parent_; }

    /// Returns the ranks of RDMs in subtree_.
    std::vector<int> required_rdm(std::vector<int> input) const;

    /// Returns depth in graph.
    int depth() const;

    /// Calls generate_task_list for subtree.
    std::tuple<OutStream, int, int, std::vector<std::shared_ptr<Tensor>>>
        generate_task_list(int tcnt, int t0, const std::list<std::shared_ptr<Tensor>> gamma, std::vector<std::shared_ptr<Tensor>> itensors) const;

};

/// Base class for specific trees (derived classes). Used in conjunction with BinaryContraction for overall graph structure. Maps derived explicit equation into tensor contraction tasks.  Task and dependency generation also initiated here.
class Tree {
  protected:
    /// Recursive structure. BinaryContraction contains Tree. A pair of tensors.
    std::list<std::shared_ptr<BinaryContraction>> bc_;
    /// Pointer to target_ tensor. Note that target_ can be NULL (at the very beginning).
    std::shared_ptr<Tensor> target_;

    /// Collection of operator tensors.
    std::vector<std::shared_ptr<Tensor>> op_;

    /// If transpose.
    bool dagger_;

    /// Points above. TODO raw pointer should be replaced by some better means.
    BinaryContraction* parent_;

    /// This is only used for finding gamma with in the tree.
    std::list<std::shared_ptr<Tensor>> gamma_;

    /// For code generation.
    std::string tree_name_;

    /// Add dependency tasks.
    std::string add_depend(const std::shared_ptr<const Tensor> o, const std::list<std::shared_ptr<Tensor>> gamma) const;

    /// When we generate, a counter is used to generate a list of tasks.
    mutable int num_;

    /// Tree label, used for graph-specific code generation.
    std::string label_;

    /// If top of tree has target indices.
    const bool root_targets_;


  public:
    /// Construct tree from equation and set tree label. Tree construction starts here.
    Tree(const std::shared_ptr<Equation> eq, std::string lab = "");
    /// Construct tree of listtensors.
    Tree(const std::shared_ptr<ListTensor> l, std::string lab, const bool t);
    virtual ~Tree() { }


    /// Return label of tree.
    virtual std::string label() const = 0;

    /// Returns depth, 0 is top of graph.
    int depth() const;

    /// Prints tree which consists of tensor binary contrations between tensors and tensor additions. If excitation targets are present,  these printed without tensor label at top of tree.
    void print() const;

    /// Returns num_. Should be greater than zero, otherwise throws an error.
    int num() const {
      if (num_ < 0) throw std::logic_error("it seems that the logic is broken - Tree::num_ is not initialized.");
      return num_;
    }

    /// Return pointer to target tensor.
    std::shared_ptr<Tensor> target() const { return target_; }
    /// Returns true if top of tree has excitation targets.
    const bool root_targets() const { return root_targets_; }
    /// Return operator tensors.
    std::vector<std::shared_ptr<Tensor>> op() const { return op_; }
    /// Return bc.
    std::list<std::shared_ptr<BinaryContraction>> bc() const { return bc_; }
    /// Return tree name for code generation.
    std::string tree_name() const { return tree_name_; }
    /// If transpose.
    bool dagger() const { return dagger_; }

    /// Factorize (i.e., perform factorize in BinaryContraction).
    void factorize();
    /// Combine trees if tensors are equal, or if have operator tensors.
    bool merge(std::shared_ptr<Tree> o);

    void set_parent(BinaryContraction* o) { parent_ = o; }
    /// Sets parent and does set_parent_subtree for bc_ (bc_ and tree alternate in graph).
    void set_parent_sub();

    /// Tidy up target tensors. Could be done on the fly. sets for all in bc_ list.
    void set_target_rec();
    /// Set target to this tensor pointer.
    void set_target(std::shared_ptr<Tensor> o) {
      target_ = o;
      for (auto i = bc_.begin(); i != bc_.end(); ++i) (*i)->set_target(o);
    };

    /// Function runs gather_gamma and find_gamma, called from top level (main.cc).
    void sort_gamma(std::list<std::shared_ptr<Tensor>> o = std::list<std::shared_ptr<Tensor>>());
    /// Updates gamma_ with a list of unique Gamma tensors.
    void find_gamma(std::shared_ptr<Tensor> o);
    /// Recursive function to collect all Gamma tensors in graph.
    std::list<std::shared_ptr<Tensor>> gather_gamma() const;
    /// Returns gamma_, list of unique Gamma tensors.
    std::list<std::shared_ptr<Tensor>> gamma() const { return gamma_; }

    /// This function returns the rank of required RDMs here + inp. Goes through bc_ and op_ tensor lists.
    std::vector<int> required_rdm(std::vector<int> inp) const;


    // code generators!
    /// Generate task and task list files.
    std::tuple<OutStream, int, int, std::vector<std::shared_ptr<Tensor>>>
        generate_task_list(int tcnt, int t0, const std::list<std::shared_ptr<Tensor>> gamma, std::vector<std::shared_ptr<Tensor>> itensors) const;
    /// Generate code by stepping through op and bc.
    std::tuple<OutStream, int, int, std::vector<std::shared_ptr<Tensor>>>
        generate_steps(const std::string indent, int tcnt, int t0, const std::list<std::shared_ptr<Tensor>> gamma, std::vector<std::shared_ptr<Tensor>> itensors) const;
    /// Generate task in dependency file with ic as task number. Caution also have a virtual generate_task.
    OutStream generate_task(const int ic, const std::vector<std::shared_ptr<Tensor>>, const std::list<std::shared_ptr<Tensor>> g, const int i0 = 0) const;

    /// Generate task for operator task (ie not a binary contraction task). Dagger arguement refers to front subtree used at top level.
    OutStream generate_compute_operators(const std::shared_ptr<Tensor>, const std::vector<std::shared_ptr<Tensor>>, const bool dagger = false) const;

    // Tree specific code generation moved to derived classes.
    /// Needed for zero level target tensors. Generates a Task '0' ie task to initialize top (zero depth) target tensor also sets up dependency queue.
    virtual OutStream create_target(const int i) const = 0;
    /// Create new tensor based on derived tree.
    virtual std::shared_ptr<Tensor> create_tensor(std::list<std::shared_ptr<const Index>>) const = 0;

    /// Generate a task. Here ip is the tag of parent, ic is the tag of this.
    virtual OutStream generate_task(const int ip, const int ic, const std::vector<std::string>, const std::string scalar = "", const int i0 = 0, bool der = false) const = 0;
    /// Generate task header.
    virtual OutStream generate_compute_header(const int, const std::list<std::shared_ptr<const Index>> ti, const std::vector<std::shared_ptr<Tensor>>, const bool = false) const = 0;
    /// Generate task footer.
    virtual OutStream generate_compute_footer(const int, const std::list<std::shared_ptr<const Index>> ti, const std::vector<std::shared_ptr<Tensor>>) const = 0;
    /// Generate Binary contraction code.
    virtual OutStream generate_bc(const std::shared_ptr<BinaryContraction>) const = 0;

};

}


#endif
