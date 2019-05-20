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
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __TENSOR_H
#define __TENSOR_H

// In order to treat Active and Op on the same footing.
// To make it easy, this is again driven by const Index::num_.
// Perhaps it's better if it is done by pointer itself, but for the time being it's fine...
//
// const Index::spin_ is not set.

// this object should generate some codes.

#include <algorithm>
#include "active.h"
#include "output.h"

namespace smith {

/// A class for Tensors. May be active_ (contain all active indices), or be merged (contain additional tensor), or have alias (equivalent tensor).
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

    /// If deriviative tensor has additional index
    std::list<std::shared_ptr<const Index>> der_;

    /// For counting Gamma tensors, used when adding all active tensor, see merge().
    mutable int num_;

    // For tensor reindexing in case of ket.
    std::map<int, int> num_map_;

  public:
    /// Constructor for intermediate tensors, and also ci tensor. todo what about scalar--needed for intermediates or ci tensor? check!
    Tensor(const double& d, const std::string& l, const std::list<std::shared_ptr<const Index>>& i)
      : factor_(d), label_(l), index_(i) { }
    /// Constructor for const operator tensor, creates index list and checks for target indices. Called from listtensor after labels are checked in listtensor constructor.
    Tensor(const std::shared_ptr<Operator> op);
    /// Constructor for const active tensor.
    Tensor(const std::shared_ptr<Active> active);
    /// Constructor for const active tensor and index list. For rdm ci derivatives have target index and kets.
    Tensor(const std::shared_ptr<Active> active, const std::list<std::shared_ptr<const Index>>& i, std::map<int, int> m = std::map<int, int>());
    ~Tensor() { }


    /// Returns list of index pointers for tensor.
    std::list<std::shared_ptr<const Index>>& index() { return index_; }
    /// Returns const list of index pointers for tensor.
    const std::list<std::shared_ptr<const Index>>& index() const { return index_; }

    /// Returns const Tensor pointer.
    const std::shared_ptr<const Tensor> merged() const { return merged_; }

    /// Returns tensor rank, cannot be called by DF tensors so far.
    int rank() const {
      if (index_.size() & 1) throw std::logic_error("Tensor::rank() cannot be called by DF tensors so far.");
      return index_.size() >> 1;
    }

    /// prints string with tensor prefactors, label, indices and those of merged and alias tensors.
    std::string str() const;
    void print(std::string indent = "") const { std::cout << indent << str() << std::endl; }
    /// Set prefactor for tensor.
    void set_factor(const double a) { factor_ = a; }
    /// Set name of scalar. Actual value is defined later on BAGEL side, eg e0.
    void set_scalar(const std::string s) { scalar_ = s; }
    /// Used to reindex tensor in absorb_ket().
    void set_index(std::list<std::shared_ptr<const Index>> i) { index_ = i; }

    /// Returns tensor prefactor.
    double factor() const { return factor_; }
    /// Returns scalar name.
    std::string scalar() const { return scalar_; }
    /// Returns tensor name.
    std::string label() const { return alias_ ? alias_->label() : label_; }
    /// Returns active tensor pointer.
    std::shared_ptr<Active> active() { return active_; }
    /// Returns const active tensor pointer.
    const std::shared_ptr<Active> active() const { return active_; }

    /// Returns true if all the indices are of active orbitals.
    bool all_active() const { return std::all_of(index_.begin(), index_.end(), [](std::shared_ptr<const Index> i){ return i->active(); }); }

    /// Used for factorization of trees.
    bool operator==(const Tensor& o) const;

    /// Adds all-active tensor to Active_.
    void merge(std::shared_ptr<Tensor> o);
    /// Sets alias used for equivalent Gamma tensor. Used in Tree::find_gamma(). The alias is given to tensor o.
    void set_alias(std::shared_ptr<Tensor> o) { alias_ = o; }
    /// if tensor is a repeat.
    bool has_alias() const { return !!alias_; }
    /// Checks if tensor is gamma.
    bool is_gamma() const { return label_.find("Gamma") != std::string::npos; }
    /// if deriviative tensor
    bool der() { return !der_.empty(); }

    /// Return number map for tensor. Originally generated when rdm is reindexed in active.
    std::map<int, int> num_map() const { return num_map_; }

    /// Generates string for constructor for tensors in Method.cc file
    std::string constructor_str_ci(const bool diagonal = false) const;
    /// Generates string for constructor for tensors in Method.cc file
    std::string constructor_str(const bool diagonal = false) const;
    /// Generates code for get_block - source block to be added later to target (move) block.
    std::string generate_get_block(const std::string, const std::string, const std::string, const bool move = false, const bool noscale = false, int number = -2, bool merged = false, const std::list<std::shared_ptr<const Index>>& mergedlist = (std::list<std::shared_ptr<const Index>>())) const;
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
    std::string generate_active_sources(const std::string indent, const std::string tag, const int ninptensors, const bool, const std::shared_ptr<Tensor>) const;
    /// Generate for loops.
    std::string generate_loop(std::string&, std::vector<std::string>&) const;
    /// Generate code for Gamma task.
    OutStream generate_gamma_sources(const int, const bool use_blas, const bool der, const std::shared_ptr<Tensor> source, const std::list<std::shared_ptr<const Index>> di) const;
    OutStream generate_gamma_header_sources(const int, const bool use_blas, const bool der, const int nindex) const;
    OutStream generate_gamma_body_sources(const int, const bool, const bool, const int, const int, const std::list<std::shared_ptr<const Index>>&, const std::shared_ptr<Tensor> source, const std::list<std::shared_ptr<const Index>> di) const;
    OutStream generate_gamma_footer_sources(const int, const bool, const bool, const int, const int, const std::list<std::shared_ptr<const Index>>&) const;

    OutStream generate_gamma(const int, const bool use_blas, const bool der) const;
    OutStream generate_gamma_header(const int, const bool use_blas, const bool der, const int nindex, const int ninptensors) const;
    OutStream generate_gamma_body(const int, const bool, const bool, const int, const int, const std::list<std::shared_ptr<const Index>>&) const;
    OutStream generate_gamma_footer(const int, const bool, const bool, const int, const int, const std::list<std::shared_ptr<const Index>>&) const;
    /// Returns Gamma number.
    int num() const { assert(is_gamma()); return num_; }
    /// Set Gamma number.
    void set_num(const int n) const { assert(is_gamma()); num_ = n; }

    /// Static function for determining canonical ordering of Tensor
    static bool comp(std::shared_ptr<const Tensor> a, std::shared_ptr<const Tensor> b) {
      bool out;
      const std::string alabel = a->label();
      const std::string blabel = b->label();
           if (alabel.find("Gamma") != std::string::npos) out = true;
      else if (blabel.find("Gamma") != std::string::npos) out = false;
      else if (alabel == "h1") out = true;
      else if (blabel == "h1") out = false;
      else if (alabel == "f1") out = true;
      else if (blabel == "f1") out = false;
      else if (alabel == "v2") out = true;
      else if (blabel == "v2") out = false;
      else if (alabel == "t2dagger") out = true;
      else if (blabel == "t2dagger") out = false;
      else if (alabel == "t2") out = true;
      else if (blabel == "t2") out = false;
      else if (alabel == "l2dagger") out = true;
      else if (blabel == "l2dagger") out = false;
      else if (alabel == "l2") out = true;
      else if (blabel == "l2") out = false;
      else {
         std::cout << alabel << " " << blabel << std::endl;
         throw std::logic_error("I have not thought about this yet - tensor.h");
      }
      return out;
    }

};

}

#endif

