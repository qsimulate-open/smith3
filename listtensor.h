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


#ifndef __LISTTENSOR_H
#define __LISTTENSOR_H

#include "tensor.h"
#include "diagram.h"
#include <memory>
#include <list>

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


  public:
    /// Constructs a list of tensors in a diagram by constructing tensors from the operators in diagram, IF they have labels.
    ListTensor(std::shared_ptr<Diagram> d);
    /// Construct listtensor using prefactor, scalar, list of tensors and dagger information.
    ListTensor(double f, std::string sc, std::list<std::shared_ptr<Tensor>> ve, bool d)
      : fac_(f), scalar_(sc), list_(ve), dagger_(d) {};
    ~ListTensor() {};
    
    /// Prints prefactor, if available: scalar, dagger. Finally prints out each tensor in list.
    void print() const;
    /// Combines tensors and removes one from list. To do this, finds active tensor then merges other tensor if other tensor is all_active (has all active indices) but not if active and if not proj. 
    void absorb_all_internal();

    /// Returns the size of the list of tensors.
    int length() const { return list_.size(); };
    /// Returns first tensor in the list of tensors.
    std::shared_ptr<Tensor> front() const { return list_.front(); };
    /// Returns the list of tensors (listtensor) minus the front tensor.
    std::shared_ptr<ListTensor> rest() const ;
    /// Creates! and returns a target tensor from the list of tensors.
    std::shared_ptr<Tensor> target() const;

    /// Returns the prefactor for listtensor.
    double fac() const { return fac_; };
    /// Returns scalar for listtensor.
    std::string scalar() const { return scalar_; };
    /// Returns dagger (if transpose) for listtensor.
    bool dagger() const { return dagger_; };
};

}

#endif
