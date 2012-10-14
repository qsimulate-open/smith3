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
// To make it easy, this is again driven by Index::num_.
// Perhaps it's better if it is done by pointer itself, but for the time being it's fine...
//
// Index::spin_ is not set.

// this object should generate some codes.

#include <memory>
#include <string>
#include <list>
#include "op.h"
#include "active.h"

namespace smith {

class Tensor {
  protected:
    // factor
    double factor_;
    // scalar to be supplied later by bagel
    std::string scalar_;
    // label of this tensor.
    std::string label_;
    // a list of indices
    std::list<std::shared_ptr<Index> > index_;

    // if this tensor is active, it has an internal structure
    std::shared_ptr<Active> active_;
    std::shared_ptr<Tensor> merged_;

  public:
    Tensor(const double& d, const std::string s, const std::string& l, const std::list<std::shared_ptr<Index> >& i)
      : factor_(d), scalar_(s), label_(l), index_(i) { };
    Tensor(const double& d, const std::string& l, const std::list<std::shared_ptr<Index> >& i)
      : factor_(d), label_(l), index_(i) { };
    Tensor(const std::shared_ptr<Op> op);
    Tensor(const std::shared_ptr<Active> active);
    Tensor() {};
    ~Tensor() {};

    std::list<std::shared_ptr<Index> >& index() { return index_; };
    const std::list<std::shared_ptr<Index> >& index() const { return index_; };

    const std::shared_ptr<const Tensor> merged() const { return merged_; };

    int rank() const {
      if (index_.size() & 1) throw std::logic_error("Tensor::rank() cannot be called by DF tensors so far.");
      return index_.size() >> 1;
    };

    std::string str() const;
    void print(std::string indent = "") const { std::cout << indent << str() << std::endl; };
    void set_factor(const double a) { factor_ = a; };
    void set_scalar(const std::string s) { scalar_ = s; };

    double factor() const { return factor_; };
    std::string scalar() const {return scalar_; };
    std::string label() const { return label_; };
    std::shared_ptr<Active> active() { return active_; };
    const std::shared_ptr<Active> active() const { return active_; };

    bool all_active() const;

    // used for factorization of trees
    bool operator==(const Tensor& o) const;

    void merge(std::shared_ptr<Tensor> o);

    std::string constructor_str(std::string indent) const;

    std::string generate_get_block(const std::string, const std::string, const bool move = false, const bool noscale = false) const;
    std::string generate_scratch_area(const std::string, const std::string, const bool zero = false) const;
    std::string generate_sort_indices(const std::string, const std::string, const std::list<std::shared_ptr<Index> >&, const bool op = false) const;
    std::string generate_sort_indices_target(const std::string, const std::string, const std::list<std::shared_ptr<Index> >&,
                                             const std::shared_ptr<Tensor>, const std::shared_ptr<Tensor>) const;

    std::pair<std::string, std::string> generate_dim(const std::list<std::shared_ptr<Index> >&) const;

    std::string generate_active(const std::string indent, const std::string tag) const;
    std::string generate_loop(std::string&, std::vector<std::string>&) const;
    std::string generate_gamma(std::string& indent, std::vector<std::string>& close, std::string tag, const bool) const;

};

}

#endif

