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

// base class for spin-summed operators

class Op {
  protected:
    // tensor info
    std::string label_;
    // this op_ is very important (and does not seem clear...).
    // Here is the convention
    // get<0>  :  Index object
    // get<1>  :  Operator info.
    //              -1: no operator (i.e., already contracted)
    //               0: operator
    //               2: active operator 
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> > op_;

    // operator info
    std::vector<std::shared_ptr<Spin> > rho_;

    std::shared_ptr<Index> a_;
    std::shared_ptr<Index> b_;
    std::shared_ptr<Index> c_;
    std::shared_ptr<Index> d_;

    // perm_count
    std::vector<int> perm_;

  public:
    Op(const std::string lab, const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td);
    Op(const std::string lab, const std::string& ta, const std::string& tb);
    // creating tensor..
    Op(const std::string lab, std::shared_ptr<Index> ta, std::shared_ptr<Index> tb, std::shared_ptr<Spin> ts = std::make_shared<Spin>());
    Op(const std::string lab = "") : label_(lab) { };
    virtual ~Op() {};

    // returns if this operator is completely contracted
    bool contracted() const;
    // returns if this operator is a general operator (i.e., Hamiltonian)
    bool general() const; 
    int num_general() const;
    int num_active_nodagger() const;
    int num_active_dagger() const;
    void mutate_general(int& i);
    int num_nodagger() const;
    int num_dagger() const;

    std::shared_ptr<Op> copy() const;
    std::pair<bool, double> permute(const bool proj);

    bool identical(std::shared_ptr<Op> o) const;

    std::string label() const { return label_; };

    void set_rho(const int i, std::shared_ptr<Spin> a) { rho_[i] = a; };
    std::vector<std::shared_ptr<Spin> >& rho() { return rho_; };
    std::shared_ptr<Spin> rho(const int i) const { return rho_.at(i); };
    const std::shared_ptr<Spin>* rho_ptr(const int i) const { return &rho_.at(i); };
    std::shared_ptr<Spin>* rho_ptr(const int i) { return &rho_.at(i); };

    const std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() const { return op_; };
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() { return op_; };

    // CAUTION:: this function returns the first daggered operator
    //           **AND** deletes the corresponding entry from this->op_.
    std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* > first_dagger_noactive();

    // returns if you can contract two labels
    bool contractable(std::string a, std::string b) { return a == b || a == "g" || b == "g"; };

    // returns which index to be kept when contraction is performed
    std::shared_ptr<Index>* survive(std::shared_ptr<Index>* a, std::shared_ptr<Index>* b);

    // perform a contraction (skipping first "skip" possibilities) and returns the factor to be multiplied.
    // Should be called from Diagram objects
    // fac, new, old
    std::tuple<double, std::shared_ptr<Spin>, std::shared_ptr<Spin> >
      contract(std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* >& dat, const int skip);

    // function to update num_ fields in Index and Spin. Should be called from Diagram objects
    void refresh_indices(std::map<std::shared_ptr<Index>, int>& dict,
                         std::map<std::shared_ptr<Index>, int>& done,
                         std::map<std::shared_ptr<Spin>, int>& spin);

    // printing this object out
    void print() const;
};



#endif
