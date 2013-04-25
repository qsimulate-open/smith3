//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tensor.h
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


#ifndef __TENSOR_H
#define __TENSOR_H

// In order to treat Active and Op on the same footing.
// To make it easy, this is again driven by const Index::num_.
// Perhaps it's better if it is done by pointer itself, but for the time being it's fine...
//
// const Index::spin_ is not set.

// this object should generate some codes.

#include <memory>
#include <string>
#include <list>
#include "operator.h"
#include "op.h"
#include "ex.h"
#include "active.h"

namespace smith {

/// A class for Tensors. May be active_ (contain all active indices), or be merged (contain additional tensor), or have have alias (equivalent tensor). 
class Tensor {
  protected:
    /// Tensor prefactor.
    double factor_;
    /// Scalar to be supplied later by BAGEL.
    std::string scalar_;
    /// Label of this tensor.
    std::string label_;
    /// List of indices.
    std::list<std::shared_ptr<const Index>> index_;

    /// If this tensor is active, it has an internal structure.
    std::shared_ptr<Active> active_;
    /// If merged, tensor should be multiplied by additional tensor.
    std::shared_ptr<Tensor> merged_;

    /// Alias tensor if any.
    std::shared_ptr<Tensor> alias_;
    /// For counting Gamma tensors, used when adding all active tensor, see merge().
    mutable int num_;


  public:
    /// Constructor for tensor with scalar.
    Tensor(const double& d, const std::string s, const std::string& l, const std::list<std::shared_ptr<const Index>>& i)
      : factor_(d), scalar_(s), label_(l), index_(i) { };
    /// Constructor for tensor without scalar.
    Tensor(const double& d, const std::string& l, const std::list<std::shared_ptr<const Index>>& i)
      : factor_(d), label_(l), index_(i) { };
    /// Constructor for const operator tensor, creates index list and checks for target indices. Called from listtensor after labels are checked in listtensor constructor. 
    Tensor(const std::shared_ptr<Operator> op);
    /// Constructor for const active tensor.
    Tensor(const std::shared_ptr<Active> active);
    /// Construct tensor with no arguements.
    Tensor() {};
    ~Tensor() {};

    /// Returns list of index pointers for tensor.
    std::list<std::shared_ptr<const Index>>& index() { return index_; };
    /// Returns const list of index pointers for tensor.
    const std::list<std::shared_ptr<const Index>>& index() const { return index_; };


    /// Returns const Tensor pointer.
    const std::shared_ptr<const Tensor> merged() const { return merged_; };
    /// Returns tensor rank, cannot be called by DF tensors so far. 
    int rank() const {
      if (index_.size() & 1) throw std::logic_error("Tensor::rank() cannot be called by DF tensors so far.");
      return index_.size() >> 1;
    };
    
    /// Returns string with tensor prefactors, label, indices and those of merged and alias tensors. 
    std::string str() const;
    void print(std::string indent = "") const { std::cout << indent << str() << std::endl; };
    /// Set prefactor for tensor.
    void set_factor(const double a) { factor_ = a; };
    /// Set name of scalar. Actual value is defined later on BAGEL side, eg e0.
    void set_scalar(const std::string s) { scalar_ = s; };

    /// Returns tensor prefactor.
    double factor() const { return factor_; };
    /// Returns scalar name.
    std::string scalar() const { return scalar_; };
    /// Returns tensor name.
    std::string label() const { return alias_ ? alias_->label() : label_; };
    /// Returns active tensor pointer.
    std::shared_ptr<Active> active() { return active_; };
    /// Returns const active tensor pointer.
    const std::shared_ptr<Active> active() const { return active_; };

    /// Returns true if all the indices are of active orbitals.
    bool all_active() const;

    /// Used for factorization of trees.
    bool operator==(const Tensor& o) const;

    /// Adds all-active tensor to Active_.
    void merge(std::shared_ptr<Tensor> o);
    /// Sets alias used for equivalent Gamma tensor. Used in Tree::find_gamma(). The alias is given to tensor o.
    void set_alias(std::shared_ptr<Tensor> o) { alias_ = o; };


    /// Generates string for constructor for tensors in Method.h file
    std::string constructor_str(std::string indent) const;
    /// Generates code for get_block - source block to be added later to target (move) block.
    std::string generate_get_block(const std::string, const std::string, const std::string, const bool move = false, const bool noscale = false) const;
    /// Generate code for unique_ptr scratch arrays.
    std::string generate_scratch_area(const std::string, const std::string, const std::string tensor_lab, const bool zero = false) const;
    /// Generate code for sort_indices. Based on operations needed to sort input tensor to output tensor. 
    std::string generate_sort_indices(const std::string, const std::string, const std::string, const std::list<std::shared_ptr<const Index>>&, const bool op = false) const;
    /// Generate code for final sort_indices back to target indices (those not summed over).
    std::string generate_sort_indices_target(const std::string, const std::string, const std::list<std::shared_ptr<const Index>>&,
                                             const std::shared_ptr<Tensor>, const std::shared_ptr<Tensor>) const;
    /// Obtain dimensions for code for tensor multiplication in dgemm.
    std::pair<std::string, std::string> generate_dim(const std::list<std::shared_ptr<const Index>>&) const;
    /// Generates code for RDMs. 
    std::string generate_active(const std::string indent, const std::string tag, const int ninptensors, const bool) const;
    /// Generate for loops.
    std::string generate_loop(std::string&, std::vector<std::string>&) const;
    /// Generate code for Gamma (overlap) task.
    std::string generate_gamma(const int, const bool) const;
    /// Returns Gamma (overlap) number.
    int num() const { assert(label_.find("Gamma") != std::string::npos); return num_; }; 
    /// Set Gamma (overlap) number.
    void set_num(const int n) const { assert(label_.find("Gamma") != std::string::npos); num_ = n; };

};

}

#endif

