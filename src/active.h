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


//
// implements the active part
//


#ifndef __ACTIVE_H
#define __ACTIVE_H

#include "op.h"
#include "rdm.h"
#include "rdm00.h"
#include "rdmI0.h"

namespace smith {

/// A class for active tensors.
class Active {
  protected:
    /// List of RDMs.
    std::list<std::shared_ptr<RDM>> rdm_;
    /// This function calls RDM::reduce_one and RDM::reduce_done functions and does sort to apply Wick's theorem to this RDM. Uses anticommutator property to rearrange indices.
    void reduce(std::shared_ptr<RDM> in);

    /// if have bra
    bool bra_;
    /// if have ket
    bool ket_;

    /// TODO double check if needed.
    mutable int count__;

    /// Map from ket reindexing.
    std::map<int, int> num_map_;


  public:
    /// Make active object from const list index and braket.
    Active(const std::list<std::shared_ptr<const Index>>& in, std::pair<bool, bool> braket);
    ~Active() { }

    /// Prints active tensor prefactor, indices and delta (equivalent indices).
    void print(const std::string& indent = "") const;
    /// Return const index list.
    const std::list<std::shared_ptr<const Index>> index() const;

    /// Used in listtensor absorb_ket.
    std::list<std::shared_ptr<RDM>> rdm() { return rdm_; }

    /// Compares active tensors. Comparison is rdm order specific now. TODO could be made more general.
    bool operator==(const Active& o) const;

    /// This generate does get_block, sort_indices, and the merged (fock) multiplication for Gamma summation.
    std::string generate(const std::string indent, const std::string tag, const std::list<std::shared_ptr<const Index>> index, const std::list<std::shared_ptr<const Index>> merged = std::list<std::shared_ptr<const Index>>(), const std::string mlab = "", const bool use_blas = false) const;
    std::string generate_sources(const std::string indent, const std::string tag, const std::list<std::shared_ptr<const Index>> index, const std::list<std::shared_ptr<const Index>> merged = std::list<std::shared_ptr<const Index>>(), const std::string mlab = "", const bool use_blas = false) const;
    /// Returns vector of int cooresponding to RDM numbers in Gamma. RDM0 is not included for non derivative trees.
    std::vector<std::string> required_rdm(const bool merged = false) const;


    /// Map from ket reindexing.
    std::map<int, int> num_map() const { return num_map_; }

    /// Merge two Active's
    void merge(std::shared_ptr<const Active> o, const double fac);

};

}

#endif
