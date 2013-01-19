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
#include "op.h"


namespace smith {

/// Class for excitation operators, derived from Op base class.  The ex operators are not directly associated with a tensor. 
class Ex : public Op {



  public:
    /// Create one-body operator. label, daggered index, partner index
    Ex(const std::string lab, const std::string& oa, const std::string& ob);
    ~Ex() {};





};

}
#endif
