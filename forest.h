//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: forest.h
// Copyright (C) 2013 Matthew MacLeod
//
// Author: Matthew MacLeod <matthew.macleod@northwestern.edu>
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


#ifndef __FOREST_H
#define __FOREST_H

#include <string>
#include <memory>
#include <list>
#include "tree.h"
#include "equation.h"
#include "listtensor.h"

namespace smith {


/// This class contains a collection of trees. 
class Forest {
  protected:
    std::list<std::shared_ptr<Tree>> trees_;   

    /// This list contains the unique gamma tensors in all trees.
    std::list<std::shared_ptr<Tensor>> gamma_;
    
    /// Name for generated code.
    std::string forest_name_;

    /// Add dependency tasks.
    std::string add_depend(const std::shared_ptr<const Tensor> o) const;

    /// When we generate, a counter is used to generate a list of tasks.
    mutable int num_;

  
  public:
    Forest(std::list<std::shared_ptr<Tree>> o) : trees_(o), forest_name_(trees_.front()->tree_name()) { }
    ~Forest() { }

    /// Function runs from top level (main.cc) adds unique gamma to gamma_ list.
    void filter_gamma();
    /// Returns the unique Gamma tensors.
    std::list<std::shared_ptr<Tensor>> gamma() const { return gamma_; }

    /// Returns name of generated code.
    std::string name() const { return forest_name_; }


    // code generation //
    /// Driver for code generation goes through trees and generates task and task list files.
    std::pair<std::string,std::string> generate_code() const;
    /// Generates headers and residual target task.
    std::pair<std::string,std::string> generate_headers() const;
    /// Generates code for all unique gamma.
    std::pair<std::string,std::string> generate_gammas() const;
    /// Generates the algorithm to be used in BAGEL.
    std::pair<std::string,std::string> generate_algorithm() const;

    /// Returns num_. Should be greater than zero, otherwise throws an error.
    int num() const {
      if (num_ < 0) throw std::logic_error("it seems that the logic is broken - Forest::num_ is not initialized.");
      return num_;
    }

};

}


#endif
