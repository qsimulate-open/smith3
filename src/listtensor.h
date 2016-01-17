//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: listtensor.h
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


#ifndef __LISTTENSOR_H
#define __LISTTENSOR_H

#include "cost.h"
#include "tensor.h"
#include "diagram.h"

namespace smith {

/// Class for a list of tensors.
class ListTensor {
  protected:
    /// Prefactor for listtensor.
    double fac_;
    /// Scalar for listtensor, to be defined in BAGEL.
    std::string scalar_;
    /// A list of tensors in a diagram.
    std::list<std::shared_ptr<Tensor>> list_;
    /// If dagger (transpose).
    bool dagger_;
    /// Braket information.
    std::pair<bool, bool> braket_;


  public:
    /// Constructs a list of tensors in a diagram by constructing tensors from the operators in diagram, IF they have labels.
    ListTensor(std::shared_ptr<Diagram> d);
    /// Construct listtensor using prefactor, scalar, list of tensors and dagger information.
    ListTensor(double f, std::string sc, std::list<std::shared_ptr<Tensor>> ve, bool d, std::pair<bool, bool> bk)
      : fac_(f), scalar_(sc), list_(ve), dagger_(d), braket_(bk) { }
    ~ListTensor() { }


    /// Prints prefactor, if available: scalar, dagger. Finally prints out each tensor in list.
    void print() const;

    /// Combines tensors and removes one from list. To do this, finds active tensor then merges other tensor if other tensor is all_active (has all active indices) but not if active and if not proj. Eg f1 tensor can be absorbed if all active.
    void absorb_all_internal();
    /// Careful, only valid if wave function is not complex. This will reverse braket for gamma and reindex tensors in case of ket, allowing gamma tensors from bra case to be reused.
    void absorb_ket();

    /// check if listtensor has rdm(s).
    bool has_gamma() const;

    /// Returns list_
    const std::list<std::shared_ptr<Tensor>>& tensors() const { return list_; }
    /// Returns the size of the list of tensors.
    int length() const { return list_.size(); }
    /// Returns first tensor in the list of tensors.
    std::shared_ptr<Tensor> front() const { return list_.front(); }
    /// Returns the list of tensors (listtensor) minus the front tensor.
    std::shared_ptr<ListTensor> rest() const ;
    /// Creates! and returns a target tensor from the list of tensors. The intermediate tensors are made here. Called from tree ctor.
    std::shared_ptr<Tensor> target() const;

    /// Returns the prefactor for listtensor.
    double fac() const { return fac_; }
    /// Returns scalar for listtensor.
    std::string scalar() const { return scalar_; }
    /// Returns braket for listtensor.
    std::pair<bool,bool> braket() const { return braket_; }
    /// set braket information for listtensor, used in absorb_ket.
    void set_braket(std::pair<bool,bool> bk) { braket_ = bk; }
    /// Returns dagger (if transpose) for listtensor.
    bool dagger() const { return dagger_; }

    /// evaluate the cost of computing this diagram as in the current order
    std::shared_ptr<Cost> calculate_cost() const;
    /// reorder the tensors so that the cost is minimal
    void reorder();
};

}

#endif
