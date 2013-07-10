//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: active.h
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


//
// implements the active part
//


#ifndef __ACTIVE_H
#define __ACTIVE_H

#include <string>
#include <list>
#include <memory>
#include "operator.h"
#include "op.h"
#include "ex.h"
#include "rdm.h"
#include "rdm00.h"
#include "rdmI0.h"

namespace smith {

/// A class for active tensors.
class Active {
  protected:
    /// List of RDMs.
    std::list<std::shared_ptr<RDM>> rdm_;
    /// This function calls RDM::reduce_one and RDM::reduce_done functions and does sort to apply Wick's theorem to this RDM. Use anticommutator property to rearrange indices
    void reduce(std::shared_ptr<RDM> in);

    /// TODO double check if needed.
    mutable int count__;

  public:
    /// Make active object from const list index and braket.
    Active(const std::list<std::shared_ptr<const Index>>& in, std::pair<bool, bool> braket);
    ~Active() { }

    /// Prints active tensor prefactor, indices and delta (equivalent indices).
    void print(const std::string& indent = "") const;
    /// Return const index list.
    const std::list<std::shared_ptr<const Index>> index() const;

    /// Compares active tensors. Comparison is rdm order specific now. TODO could be made more general.
    bool operator==(const Active& o) const;

    /// This generate does get_block, sort_indices, and the merged (fock) multiplication for Gamma summation.
    std::string generate(const std::string indent, const std::string tag, const std::list<std::shared_ptr<const Index>> index, const std::list<std::shared_ptr<const Index>> merged = std::list<std::shared_ptr<const Index>>(), const std::string mlab = "", const bool use_blas = false) const;
    /// Returns vector of int cooresponding to RDM numbers in Gamma. RDM0 is not included.
    std::vector<int> required_rdm() const;

};

}

#endif
