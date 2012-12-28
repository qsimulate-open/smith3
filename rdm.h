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

#include <list>
#include <memory>
#include <string>
#include <cassert>
#include <map>
#include "index.h"

namespace smith {

class RDM {
  protected:
    /// Prefactor for RDM.
    double fac_;
    /// Operators that constitute RDM.
    std::list<std::shared_ptr<Index> > index_;
    /// Kronecker's delta.
    std::map<std::shared_ptr<Index>, std::shared_ptr<Index> > delta_;
   
    /// Generate get block - source data to be added to target (move block).
    std::string make_get_block(std::string indent, std::string tag, std::string lbl);
    /// Generate unique_ptr with hash TODO I think this is obsolete now.
    std::string make_scratch_area(std::string indent, std::string tag, std::string lbl);
    /// Generate sort_indices which makes array. This version has no addition (or factor multiplication-0111).
    std::string make_sort_indices(std::string indent, std::string tag, const std::list<std::shared_ptr<Index> >& loop);
    /// Do blas multiplication of Gamma and fock tensors...not implemented yet for subtask code!
    std::string make_blas_multiply(std::string indent, const std::list<std::shared_ptr<Index> >& loop, const std::list<std::shared_ptr<Index> >& index);
    /// Used for blas multiplication of RDM and merged (fock) tensors. NB not implemented yet for subtask code!
    std::pair<std::string, std::string> get_dim(const std::list<std::shared_ptr<Index> >& di, const std::list<std::shared_ptr<Index> >& index) const;
    /// If delta case, also makes index loops then checks to see if merged-or-delta indices are in loops..
    std::string make_merged_loops(std::string& indent, const std::string tag, std::vector<std::string>& close);
    /// Generates RDM and merged (fock) tensors multipication.
    std::string multiply_merge(const std::string itag, std::string& indent,  const std::list<std::shared_ptr<Index> >& merged);
    /// Generates odata (Gamma) part of for summation ie LHS in equations gamma += rdm or gamma += rdm * f1
    std::string make_odata(const std::string itag, std::string& indent, const std::list<std::shared_ptr<Index> >& index);
    /// Loops over delta indices in Gamma summation.
    std::string make_sort_loops(const std::string itag, std::string& indent, const std::list<std::shared_ptr<Index> >& index, std::vector<std::string>& close);
    /// Makes if statement in delta cases ie index equivalency check line.
    std::string make_delta_if(std::string& indent, std::vector<std::string>& close);
    /// Adds merged (fock) tensor with indices, used by muliply_merge member.
    std::string fdata_mult(const std::string itag, const std::list<std::shared_ptr<Index> >& merged);
  
    /// Used to replace tensor labels to more general labels in(x). RDM tensors are always numbered before merged(fock) tensor. For example, rdm1 will be mapped to in(0), rdm2 -> in(1), and in merged case with max rdm2, f1 -> in(2).
    void map_in_tensors(std::vector<std::string> in_tensors, std::map<std::string,std::string>& inlab);

    /// Generate code for Gamma RDM summation.
    std::string generate_not_merged(std::string indent, const std::string tlab, const std::list<std::shared_ptr<Index> >& loop, std::vector<std::string> in_tensors);
    /// Generates code for Gamma RDM summation with merged object (additional tensor, here fock tensor) multiplication.
    std::string generate_merged(std::string indent, const std::string itag, const std::list<std::shared_ptr<Index> >& index, const std::list<std::shared_ptr<Index> >& merged, const std::string mlab, std::vector<std::string> in_tensors, const bool use_blas);

  public:
    RDM(const std::list<std::shared_ptr<Index> >& in,
        const std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& in2,
        const double& f = 1.0)
      : fac_(f), index_(in), delta_(in2) { };
    ~RDM() {};


    /// Prints RDM with indentation and prefactors, located in active.cc
    void print(const std::string& indent = "") const;
    /// Sort indices so that it will be 0+0 1+1 ... (spin ordering is arbitrary).
    void sort();

    /// Copies RDM located in active.cc
    std::shared_ptr<RDM> copy() const;

    /// Returns private members.
    double factor() const { return fac_; };
    double& fac() { return fac_; };
    /// Returns a reference of index_.
    std::list<std::shared_ptr<Index> >& index() { return index_; };
    /// Returns a const reference of index_.
    const std::list<std::shared_ptr<Index> >& index() const { return index_; };
    /// Returns a const reference of delta_.
    const std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& delta() const { return delta_; };
    /// Returns a reference of delta_.
    std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& delta() { return delta_; };

    /// Returns if this is in the final form..ie aligned as a0+ a0 a1+ a1..Member function located in active.cc
    bool done() const;
    /// Checks if there is an annihilation operator with creation operators to the right hand side.
    bool reduce_done(const std::list<int>& done) const;

    /// Application of Wick's theorem and is controlled by Index::num_. See active.cc. One index is going to be annihilated. done is updated inside the function.
    std::list<std::shared_ptr<RDM> > reduce_one(std::list<int>& done) const;

    /// Compares RDMs, based on indices.
    bool operator==(const RDM& o) const;

    /// Generate Gamma summation code, for both non-merged and merged case (RDM * f1 tensor multiplication).
    std::string generate(std::string indent, const std::string itag, const std::list<std::shared_ptr<Index> >& index, const std::list<std::shared_ptr<Index> >& merged, const std::string mlab, std::vector<std::string> in_tensors, const bool use_blas);

    /// Returns an integer representing rdm rank value.
    int rank() const { assert(index_.size()%2 == 0); return index_.size()/2; }; 
};

}

#endif
