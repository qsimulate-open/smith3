//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: rdm.h
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

#ifndef __SRC_RDM_H
#define __SRC_RDM_H

#include <map>
#include "index.h"

namespace smith {

/// Abstract base class for reduced density matrices (RDMs).
class RDM {
  protected:
    /// Prefactor for RDM.
    double fac_;
    /// Operators that constitute RDM.
    std::list<std::shared_ptr<const Index>> index_;
    /// Kronecker's delta, map with two index pointers.
    std::map<std::shared_ptr<const Index>, std::shared_ptr<const Index>> delta_;

    /// Inherits bra from diagram, done in active ctor.
    bool bra_;
    /// Inherits ket from diagram, done in active ctor.
    bool ket_;

    /// Generate entire task code for Gamma RDM summation.
    virtual std::string generate_not_merged(std::string indent, const std::string tlab, const std::list<std::shared_ptr<const Index>>& loop, std::vector<std::string> in_tensors) = 0;
    /// Generates entire task code for Gamma RDM summation with merged object (additional tensor, here fock tensor) multiplication.
    virtual std::string generate_merged(std::string indent, const std::string itag, const std::list<std::shared_ptr<const Index>>& index, const std::list<std::shared_ptr<const Index>>& merged, const std::string mlab, std::vector<std::string> in_tensors, const bool use_blas) = 0;

    /// Makes if statement in delta cases ie index equivalency check line.
    virtual std::string make_delta_if(std::string& indent, std::vector<std::string>& close) = 0;
    /// Replaces tensor labels to more general labels in(x), where x is a counter for in tensors. RDM tensors numbered before merged (fock) tensor. Eg, rdm1 is mapped to in(0), rdm2 -> in(1), and in merged case with max rdm2, f1 -> in(2).
    virtual void map_in_tensors(std::vector<std::string> in_tensors, std::map<std::string,std::string>& inlab) = 0;

    /// Loops over delta indices in Gamma summation.
    virtual std::string make_sort_loops(const std::string itag, std::string& indent, const std::list<std::shared_ptr<const Index>>& index, std::vector<std::string>& close) = 0;

    // for task summation line
    /// Generates odata (Gamma) part of for summation ie LHS in equations gamma += rdm or gamma += rdm * f1
    virtual std::string make_odata(const std::string itag, std::string& indent, const std::list<std::shared_ptr<const Index>>& index) = 0;

    /// Do blas multiplication of Gamma and fock tensors...not implemented yet for subtask code!
    virtual std::string make_blas_multiply(std::string indent, const std::list<std::shared_ptr<const Index>>& loop, const std::list<std::shared_ptr<const Index>>& index) = 0;
    /// Used for blas multiplication of RDM and merged (fock) tensors. NB not implemented yet for subtask code!
    virtual std::pair<std::string, std::string> get_dim(const std::list<std::shared_ptr<const Index>>& di, const std::list<std::shared_ptr<const Index>>& index) const = 0;


  public:
    /// Make RDM object from list of indices, delta indices and factor.
    RDM(const std::list<std::shared_ptr<const Index>>& in,
        const std::map<std::shared_ptr<const Index>, std::shared_ptr<const Index>>& in2, std::pair<bool, bool> braket,
        const double& f = 1.0)
      : fac_(f), index_(in), delta_(in2), bra_(braket.first), ket_(braket.second) { }
    virtual ~RDM() { }


    /// Prints RDM with prefactors and braket.
    void print(const std::string& indent = "") const;

    /// Sort indices so that it will be 0+0 1+1 ... (spin ordering is arbitrary).
    void sort();

    /// Returns the factor
    double factor() const { return fac_; }
    /// Returns a reference to the factor
    double& fac() { return fac_; }

    /// Returns bra, false for zero order reference.
    bool bra() const { return bra_; }
    /// Returns ket, false for zero order reference.
    bool ket() const { return ket_; }
    /// Returns bra and ket pair.
    std::pair<bool, bool> braket() { return std::make_pair(bra_,ket_); }

    /// Set bra, needed if ket is absorbed for diagram to reuse modified rdms.
    void set_bra(bool b) { bra_ = b; }
    /// Set ket, needed if ket is absorbed for diagram to reuse modified rdms.
    void set_ket(bool b) { ket_ = b; }

    /// Returns a reference of index_.
    std::list<std::shared_ptr<const Index>>& index() { return index_; }
    /// Returns a const reference of index_.
    const std::list<std::shared_ptr<const Index>>& index() const { return index_; }

    /// Returns a const reference of delta_.
    const std::map<std::shared_ptr<const Index>, std::shared_ptr<const Index>>& delta() const { return delta_; }
    /// Returns a reference of delta_.
    std::map<std::shared_ptr<const Index>, std::shared_ptr<const Index>>& delta() { return delta_; }

    /// Returns if this is in the final form..ie aligned as a0+ a0 a1+ a1..Member function located in active.cc
    bool done() const;
    /// Checks if there is an annihilation operator with creation operators to the right hand side.
    bool reduce_done(const std::list<int>& done) const;

    /// Returns an integer representing rdm rank value, ie (index size)/2
    int rank() const { assert(index_.size()%2 == 0); return index_.size()/2; }

    /// Compares for equivalency based on prefactor, indices, delta, and braket.
    bool operator==(const RDM& o) const;
    bool identical(std::shared_ptr<const RDM> o) const;

    // virtual public functions
    /// Application of Wick's theorem and is controlled by const Index::num_. See active.cc. One index is going to be annihilated. done is updated inside the function.
    virtual std::list<std::shared_ptr<RDM>> reduce_one(std::list<int>& done) const = 0;

    /// Generate Gamma summation task, for both non-merged and merged case (RDM * f1 tensor multiplication).
    virtual std::string generate(std::string indent, const std::string itag, const std::list<std::shared_ptr<const Index>>& index, const std::list<std::shared_ptr<const Index>>& merged, const std::string mlab, std::vector<std::string> in_tensors, const bool use_blas) = 0;

    /// Copies this rdm, needs to be virtual as creates derived rdms.
    virtual std::shared_ptr<RDM> copy() const = 0;

    /// Reverse order and dagger information for list. Only used in derivative cases.
    virtual std::list<std::shared_ptr<const Index>> conjugate() = 0;


};

}

#endif
