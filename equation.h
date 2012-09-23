//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: equation.h
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


#ifndef __EQUATION_H
#define __EQUATION_H

#include <list>
#include "diagram.h"

namespace smith {

class Equation {
  protected:
    // list of diagrams
    std::list<std::shared_ptr<Diagram> > diagram_;
    // internal function used by factorize()
    void duplicates_(const bool);

    std::string name_;

  public:
    Equation(std::shared_ptr<Diagram>, std::string nam);
    ~Equation() {};

    // merging two set of Equation.
    void merge(const std::shared_ptr<Equation> o) {
      diagram_.insert(diagram_.end(), o->diagram_.begin(), o->diagram_.end());
    };

    // print function
    void print();
    // active parts are processed
    void active();
    // identifys the same terms
    void duplicates();
    // refresh indices in each diagram
    void refresh_indices();

    // returns the name of this equation
    std::string name() const { return name_; };

    std::list<std::shared_ptr<Diagram> > diagram() { return diagram_; };

};

}

#endif

