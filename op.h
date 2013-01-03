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


#ifndef __TWOOP_H
#define __TWOOP_H

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <list>
#include <tuple>
#include "index.h"


namespace smith {

/// Base class for spin-summed operators. Operator class works on Index.
class Op {
  protected:
    /// Related to tensor info. 
    std::string label_;
    // this op_ is very important (and does not seem clear...).
    // Here is the convention
    // get<0>  :  Index object
    // get<1>  :  Operator info.
    //              -1: no operator (i.e., already contracted)
    //               0: operator
    //               2: active operator 
    /// Tuple with index object pointer, operator info, and spin info. Operator info is defined as -1: no operator (i.e., already contracted), 0: operator, 2: active operator.
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> > op_;

    /// Spin operator info.
    std::vector<std::shared_ptr<Spin> > rho_;

    /// First excitation index.
    std::shared_ptr<Index> a_;
    /// First excitation index partner.
    std::shared_ptr<Index> b_;
    /// Second excitation index.
    std::shared_ptr<Index> c_;
    /// Second excitation index partner.
    std::shared_ptr<Index> d_;

    /// This is permutation count. 
    std::vector<int> perm_;

  public:
    /// Create two-body tensor operator.
    Op(const std::string lab, const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td);
    /// Create one-body tensor operator.
    Op(const std::string lab, const std::string& ta, const std::string& tb);
    /// Create one-body tensor operator with spin information.
    Op(const std::string lab, std::shared_ptr<Index> ta, std::shared_ptr<Index> tb, std::shared_ptr<Spin> ts = std::make_shared<Spin>());
    /// Create operator with label.
    Op(const std::string lab = "") : label_(lab) { };
    virtual ~Op() {};

    /// Returns if this operator is completely contracted.
    bool contracted() const;
    /// Returns if this operator is a general operator (i.e., Hamiltonian), checks for non-zero count of general operators.
    bool general() const; 
    /// Counts number of general operators.
    int num_general() const;

    /// Counts number of nondaggered active operators.
    int num_active_nodagger() const;
    /// Counts number of daggered active operators.
    int num_active_dagger() const;
    /// Change general operator to active operator...
    void mutate_general(int& i);
    /// Counts number of nondaggered operators.
    int num_nodagger() const;
    /// Counts number of daggered operators.
    int num_dagger() const;

    /// Creates a new Op pointer.
    std::shared_ptr<Op> copy() const;
    /// Makes a possible permutation of indices.
    std::pair<bool, double> permute(const bool proj);

    /// Checks label, and first two operator tuple fields (index and operator contraction info). **NOTE** that spin info (third op field) is not checked.
    bool identical(std::shared_ptr<Op> o) const;

    /// Returns operator name.
    std::string label() const { return label_; };

    /// Set spin.
    void set_rho(const int i, std::shared_ptr<Spin> a) { rho_[i] = a; };
    /// Returns spin.
    std::vector<std::shared_ptr<Spin> >& rho() { return rho_; };
    /// Returns const spin.
    std::shared_ptr<Spin> rho(const int i) const { return rho_.at(i); };
    /// Returns a const pointer to spin.
    const std::shared_ptr<Spin>* rho_ptr(const int i) const { return &rho_.at(i); };
    /// Returns a pointer to spin.
    std::shared_ptr<Spin>* rho_ptr(const int i) { return &rho_.at(i); };

    /// Return const operator reference.
    const std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() const { return op_; };
    /// Return operator reference.
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() { return op_; };


    /// CAUTION:: this function returns the first daggered operator (not an active operator) **AND** deletes the corresponding entry from this->op_, by marking as contracted.
    std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* > first_dagger_noactive();

    /// Perform a contraction, skipping first "skip" from equation. Returns the factor, new, and old spin. Called from Diagram::reduce_one_noactive.
    std::tuple<double, std::shared_ptr<Spin>, std::shared_ptr<Spin> >
      contract(std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* >& dat, const int skip);

    /// Returns if you can contract two labels. Labels (type) must be same or one must be of type general in order to do contraction.
    bool contractable(std::string a, std::string b) { return a == b || a == "g" || b == "g"; };

    /// Returns which index to be kept when contraction is performed.
    std::shared_ptr<Index>* survive(std::shared_ptr<Index>* a, std::shared_ptr<Index>* b);

    /// Function to update Index and Spin and check if contracted.  Should be called from Diagram objects.
    void refresh_indices(std::map<std::shared_ptr<Index>, int>& dict,
                         std::map<std::shared_ptr<Index>, int>& done,
                         std::map<std::shared_ptr<Spin>, int>& spin);


    /// Printing this object out.
    void print() const;
};

}

#endif
