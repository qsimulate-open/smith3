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
#include "op.h"
#include "tensor.h"

class Tensor;

class RDM {
  friend class Tensor;
  protected:
    // prefactor
    double fac_;
    // operators that constitute RDM
    std::list<std::shared_ptr<Index> > index_;
    // kronecker's delta
    std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > > delta_;

  public:
    RDM(const std::list<std::shared_ptr<Index> >& in,
        const std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > >& in2,
        const double& f = 1.0)
      : index_(in), delta_(in2), fac_(f) { };
    ~RDM() {};

    void print(const std::string& indent = "") const;
    // sort indices so that it will be 0+0 1+1 ... (spin ordering is arbitrary)
    void sort();

    std::shared_ptr<RDM> copy() const;

    // returns private members
    double factor() const { return fac_; };
    double& fac() { return fac_; };
    // returns a reference of index_
    std::list<std::shared_ptr<Index> >& index() { return index_; };
    // returns a const reference of index_
    const std::list<std::shared_ptr<Index> >& index() const { return index_; };
    // returns a const reference of delta_
    const std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > >& delta() const { return delta_; };
    // returns a reference of delta_
    std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > >& delta() { return delta_; };

    // returns if this is in the final form
    bool done() const;
    bool reduce_done(const std::list<int>& done) const;

    // One index is going to be annihilated. done is updated inside the function
    std::list<std::shared_ptr<RDM> > reduce_one(std::list<int>& done) const;

    // generate a code
    std::string str(std::string target, std::shared_ptr<Tensor> merged) const;
};


class Active {
  friend class Tensor;
  protected:
    std::list<std::shared_ptr<RDM> > rdm_;
    void reduce(std::shared_ptr<RDM> in);
    // label of this tensor.
    std::string label_;
    // a list of indices
    std::list<std::shared_ptr<Index> > index_;

    mutable int count__;

  public:
    Active(const std::list<std::shared_ptr<Index> >& in);
    ~Active() {};

    void print(const std::string& indent = "") const;
    const std::list<std::shared_ptr<Index> > index() const;

    std::string generate(std::shared_ptr<Tensor> merged) const;
    std::string generate_get_block(const std::string, const std::string, const bool move, const std::shared_ptr<Tensor> gammat) const;

};

#endif
