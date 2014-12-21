//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: residual.h
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


#ifndef __RESIDUAL_H
#define __RESIDUAL_H

#include <string>
#include <memory>
#include "equation.h"
#include "listtensor.h"
#include "tree.h"

namespace smith {

/// Derived class for residual tree. Solve for T amplitudes iteratively.
class Residual : public Tree {
  protected:


  public:
    /// Construct tree of equation pointers and set tree label.
    Residual(const std::shared_ptr<Equation> eq, std::string lab = "") : Tree(eq, lab) { }
    Residual(const std::shared_ptr<ListTensor> l, std::string lab, const bool t) : Tree(l, lab, t) { }
    virtual ~Residual() { }


    /// Return label of tree.
    std::string label() const override { return label_; }

    OutStream create_target(const std::string, const int) const override;
    std::shared_ptr<Tensor> create_tensor(std::list<std::shared_ptr<const Index>>) const override;

    OutStream generate_task(const std::string, const int ip, const int ic, const std::vector<std::string>, const std::string scalar = "", const int i0 = 0, bool der = false) const override;
    OutStream generate_compute_header(const int, const std::list<std::shared_ptr<const Index>> ti, const std::vector<std::shared_ptr<Tensor>>, const bool = false) const override;
    OutStream generate_compute_footer(const int, const std::list<std::shared_ptr<const Index>> ti, const std::vector<std::shared_ptr<Tensor>> ) const override;
    OutStream generate_bc(const std::string, const std::shared_ptr<BinaryContraction>) const override;


};

}


#endif
