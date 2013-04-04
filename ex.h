//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: ex.h
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


#ifndef __EX_H
#define __EX_H

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
#include "operator.h"


namespace smith {

/// Derived class for excitation operators.  The Ex operators are not directly associated with a tensor. 
class Ex : public Operator {

  public:
    /// Create one-body operator. daggered index, partner index
    Ex(const std::string& oa, const std::string& ob);
    /// Create two-body operator. Projection operator should be replaced via this.
    Ex(const std::string& oa, const std::string& ob, const std::string& oc, const std::string& od);
    virtual ~Ex() {};

    /// Print out operator.
    void print() const override;

    bool is_ex() const { return true; }

    /// Makes a possible permutation of indices. Cannot permute if there are active daggered and no-daggered operators.
    std::pair<bool, double> permute(const bool proj) override;

    /// checks first two operator tuple fields (index and operator contraction info). **NOTE** that spin info (third op field) is not checked.
    bool identical(std::shared_ptr<Operator> o) const override;

    /// Creates a new Ex pointer.
    std::shared_ptr<Operator> copy() const override;

    /// Excitation operators have no label as not associated with tensor.
    std::string label() const override { return ""; };

};

}
#endif
