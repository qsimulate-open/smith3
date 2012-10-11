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
#include "op.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <map>

namespace smith {

class Diagram {
  protected:
    // a list of operators
    std::list<std::shared_ptr<Op> > op_;
    // a constant factor
    double fac_;
    // a scalar to be defined later on bagel side
    std::string sclr_;
    // active part
    std::shared_ptr<Active> rdm_;

    // if this Diagram has a daggered counterpart (often the case for residual equations).
    bool dagger_;

  public:
    //Diagram(std::list<std::shared_ptr<Op> > op) : op_(op), fac_(1.0), dagger_(false) { };
    //Diagram() : fac_(1.0), dagger_(false) { };
    Diagram(std::list<std::shared_ptr<Op> > op) : op_(op), fac_(1.0), sclr_(""), dagger_(false) { };
    Diagram() : fac_(1.0), sclr_(""), dagger_(false) { };
    // copy constructor is complicated but preserves the same topology as this.
    ~Diagram() {};

    // returns a shared_ptr of a diagram that has the same topology as this.
    std::shared_ptr<Diagram> copy() const;

    // generate all combination of diagrams (related to general indices)
    std::list<std::shared_ptr<Diagram> > get_all() const;

    // get functions
    double& fac() { return fac_; };
    const double fac() const { return fac_; };
    std::string& sclr() { return sclr_; }; 
    std::shared_ptr<Active> rdm() { return rdm_; };
    bool dagger() const { return dagger_; };

    // careful returns a const reference
    const std::list<std::shared_ptr<Op> >& op() const { return op_; };
    // set functions for private members
    void set_op(const std::list<std::shared_ptr<Op> >& o) { op_ = o; };
    void set_fac(const double a) { fac_ = a; };
    void set_sclr(std::string s) { sclr_ = s; };

    // refresh the indices
    void refresh_indices();

    // processes the active part
    void active();

    // daggered Diagram added to the sum
    void add_dagger() { dagger_ = true; };

    // permute indices in operators. return false when finished
    bool permute(const bool proj); 
    bool identical(std::shared_ptr<Diagram> o) const;

    // printing function
    // CAUTION: it also refreshes the indices
    void print();
    // this version does not refresh indices
    void print() const;

    // the number of daggered indices
    int num_dagger() const;
    // the number of general indices
    int num_general() const;
    // returns if this diagram has a consistent set of dagger and undaggered indices
    bool consistent_indices() const;

    // this function performs one contraction ** IN PLACE **
    bool reduce_one_noactive(const int skip);

    // returns if this diagram is still valid
    bool valid() const;
    // returns if this diagram is fully contracted and sorted
    bool done() const;
    // returns if this diagram is fully contracted (looking up only nonactive parts) 
    bool done_noactive() const;

    // gathers active indices
    std::list<std::shared_ptr<Index> > active_indices() const; 
};

}

#endif
