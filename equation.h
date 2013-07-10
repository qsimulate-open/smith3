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

/// This class enables work on a collection of diagrams.
class Equation {
  protected:
    /// List of diagrams.
    std::list<std::shared_ptr<Diagram>> diagram_;

    /// Internal function used by to find double permutations. Also if bool is true, equations can become daggered, and where permutations cause an additional get_block to be added to the generated code.
    void duplicates_(const bool);

    /// Name of theory. Generated code for BAGEL will have this name, set in main.cc.
    std::string name_;


  public:
    /// Construct equation from diagram and name. Contract operators in diagram.
    Equation(std::shared_ptr<Diagram>, std::string nam);
    ~Equation() { }


    /// Merging two sets of Equation.  Done by adding the diagrams of new equation (in merge arguement) to the original equation.
    void merge(const std::shared_ptr<Equation> o) {
      diagram_.insert(diagram_.end(), o->diagram_.begin(), o->diagram_.end());
    };

    /// Prunes equation to those terms containing target indices.
    void term_select(std::list<std::string> t);

    /// Check for targets.
    bool ex_targets() const;

    /// Print function. This triggers Diagram::refresh_indices().
    void print();
    /// The active parts are processed.
    void active();
    /// Identifies terms which are the same (via permuations) and updating term prefactor. Also can add a dagger to diagram.
    void duplicates();
    /// Refresh indices in each diagram.
    void refresh_indices();

    /// Returns the name of this equation.
    std::string name() const { return name_; }

    /// Returns list of diagram pointers.
    std::list<std::shared_ptr<Diagram>> diagram() { return diagram_; }

};

}

#endif

