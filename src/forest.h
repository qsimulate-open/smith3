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
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __FOREST_H
#define __FOREST_H

#include "tree.h"

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
    /// This is a zero level task for a tree.
    mutable int i0;
    /// Task counter.
    mutable int icnt;

    /// Intermediate tensors
    mutable std::vector<std::shared_ptr<Tensor>> itensors_;

    /// Returns the main body of the CASPT2 driver
    static std::string caspt2_main_driver_();
    /// Returns the main body of the MS-MRCI driver
    static std::string msmrci_main_driver_();

  public:
    Forest(std::list<std::shared_ptr<Tree>> o) : trees_(o), forest_name_(trees_.front()->tree_name()) { }

    /// Function runs from top level (main.cc) adds unique gamma to gamma_ list.
    void filter_gamma();
    /// Returns the unique Gamma tensors.
    std::list<std::shared_ptr<Tensor>> gamma() const { return gamma_; }

    /// Returns name of generated code.
    std::string name() const { return forest_name_; }
    /// Returns intermediate tensors.
    std::vector<std::shared_ptr<Tensor>> itensors() const { return itensors_; }

    // code generation //
    /// Driver for code generation goes through trees and generates task and task list files.
    OutStream generate_code() const;
    /// Generates headers and residual target task.
    OutStream generate_headers() const;
    /// Generates code for all unique gamma.
    OutStream generate_gammas() const;
    /// Generates the algorithm to be used in BAGEL.
    OutStream generate_algorithm() const;

    /// Returns num_. Should be greater than zero, otherwise throws an error.
    int num() const {
      if (num_ < 0) throw std::logic_error("it seems that the logic is broken - Forest::num_ is not initialized.");
      return num_;
    }

};

}


#endif
