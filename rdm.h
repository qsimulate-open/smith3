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
    // prefactor
    double fac_;
    // operators that constitute RDM
    std::list<std::shared_ptr<Index> > index_;
    // kronecker's delta
    std::map<std::shared_ptr<Index>, std::shared_ptr<Index> > delta_;

  public:
    RDM(const std::list<std::shared_ptr<Index> >& in,
        const std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& in2,
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
    const std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& delta() const { return delta_; };
    // returns a reference of delta_
    std::map<std::shared_ptr<Index>, std::shared_ptr<Index> >& delta() { return delta_; };

    // returns if this is in the final form
    bool done() const;
    bool reduce_done(const std::list<int>& done) const;

    // One index is going to be annihilated. done is updated inside the function
    std::list<std::shared_ptr<RDM> > reduce_one(std::list<int>& done) const;

    // generate a code
    std::string generate(std::string indent, const std::string tlab,const std::list<std::shared_ptr<Index> >& loop) const;

    int rank() const { assert(index_.size()%2 == 0); return index_.size()/2; }; 
};

}

#endif
