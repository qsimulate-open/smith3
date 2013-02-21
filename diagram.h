//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: diagram.h
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


#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "active.h"
#include "operator.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <map>

namespace smith {

/// This class is used on a collection of operators.
class Diagram {
  protected:
    /// A list of operators.
    std::list<std::shared_ptr<Operator>> op_;

    /// A constant factor.
    double fac_;
    /// A scalar to be defined later on BAGEL side.
    std::string scalar_;
    /// The active part.
    std::shared_ptr<Active> rdm_;

    /// If this Diagram has a daggered counterpart (often the case for residual equations).
    bool dagger_;

  public:
    /// Construct diagram from operator list, prefactor and dagger information.
    Diagram(std::list<std::shared_ptr<Operator>> op) : op_(op), fac_(1.0), dagger_(false) { };
    /// Construct diagram from operator list, prefactor, scalar, and dagger information.
    Diagram(std::list<std::shared_ptr<Operator>> op, std::string s) : op_(op), fac_(1.0), scalar_(s), dagger_(false) { };
    /// Construct diagram with prefactor and dagger information.
    Diagram() : fac_(1.0), dagger_(false) { };
    // copy constructor is complicated but preserves the same topology as this.
    ~Diagram() {};

    /// Returns a shared_ptr of a diagram that has the same topology as this.
    std::shared_ptr<Diagram> copy() const;

    /// Generate all combination of diagrams (related to general indices).
    std::list<std::shared_ptr<Diagram>> get_all() const;

    // Get functions.
    /// Return the diagram (term) prefactor.
    double& fac() { return fac_; };
    /// Return the prefactor for const diagram.
    const double fac() const { return fac_; };
    /// Return scalar name reference.
    std::string& scalar() { return scalar_; }; 
    /// Returns rdm pointer.
    std::shared_ptr<Active> rdm() { return rdm_; };
    /// If diagram is transposed.
    bool dagger() const { return dagger_; };

    /// Careful, returns a const reference of op_ operator.
    const std::list<std::shared_ptr<Operator>>& op() const { return op_; };
    /// Set operator for private members.
    void set_op(const std::list<std::shared_ptr<Operator>>& o) { op_ = o; };
    /// Set factor for private members.
    void set_fac(const double a) { fac_ = a; };

    /// Refresh the indices for each operator in diagram (ie calls operators refresh_indices function).
    void refresh_indices();

    /// Processes the active part. Performs Wick's in constructor of an Active object.
    void active();

    /// Daggered Diagram added to the sum.
    void add_dagger() { dagger_ = true; };

    /// Permute indices in operators. return false when finished.
    bool permute(const bool proj); 
    /// If diagrams are same, based on size, indices, and spin. 
    bool identical(std::shared_ptr<Diagram> o) const;

    /// Print function for diagram, CAUTION: it also refreshes the indices.
    void print();
    /// This print version does not refresh indices. Prints factor, scalar and operators.
    void print() const;

    /// The number of daggered indices.
    int num_dagger() const;
    /// The number of general indices.
    int num_general() const;
    /// Returns if this diagram has a consistent set of dagger and undaggered indices.
    bool consistent_indices() const;

    /// This function performs one contraction ** IN PLACE ** called from equation.
    bool reduce_one_noactive(const int skip);

    /// Returns if this diagram is still valid.
    bool valid() const;
    /// Returns if this diagram is fully contracted and sorted.
    bool done() const;
    /// Returns if this diagram is fully contracted (looking up only nonactive parts).
    bool done_noactive() const;

    /// Gathers active indices.
    std::list<std::shared_ptr<Index>> active_indices() const; 
};

}

#endif
