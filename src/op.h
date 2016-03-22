//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: op.h
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


#ifndef __TWOOP_H
#define __TWOOP_H

#include "operator.h"

namespace smith {

/// Derived class for spin-summed operators which produce tensors.
/// for a historical reason, this is a derived class of Operator, although it does not have to be so.
class Op : public Operator {
  protected:
    /// Related to tensor info.
    std::string label_;

  public:
    /// Create two-body tensor operator.
    Op(const std::string lab, const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td,
       const bool alpha1 = false, const bool alpha2 = false);
    /// Create one-body tensor operator. When alpha is true, the operators are alpha only.
    Op(const std::string lab, const std::string& ta, const std::string& tb, const bool alpha = false);
    /// Create one-body tensor with spin information. No operator is created.
    Op(const std::string lab, std::shared_ptr<Index> ta, std::shared_ptr<Index> tb, std::shared_ptr<Spin> ts);
    /// Create operator with label.
    Op(const std::string lab = "") : label_(lab) { }

    /// Create operator without label.
    Op(const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td, const bool alpha1, const bool alpha2)
     : Op("", ta, tb, tc, td, alpha1, alpha2) { }
    /// Create operator without label. When alpha is true, the operators are alpha only.
    Op(const std::string& ta, const std::string& tb, const bool alpha = false) : Op("", ta, tb, alpha) { }

    virtual ~Op() { }

    bool is_ex() const { return label_.empty(); }

    /// Returns operator name.
    std::string label() const override { return label_; }

    /// Print out operator.
    void print() const override;

    /// Makes a possible permutation of indices. Cannot permute if there are active daggered and no-daggered operators or if label is proj.
    std::pair<bool, double> permute(const bool proj) override;

    /// Checks label, and first two operator tuple fields (index and operator contraction info). **NOTE** that spin info (third op field) is not checked.
    bool identical(std::shared_ptr<Operator> o) const override;

    /// Creates a new Operator pointer.
    std::shared_ptr<Operator> copy() const override;


};

}

#endif
