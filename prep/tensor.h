//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: main.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the SMITH3 package.
//
// The SMITH3 package is free software; you can redistribute it and/or modify
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

//////////////////////////////////////////////////////////////////////////////////
//
//  This program generates a simple CASPT2 theory generator
//  
//  Here the CASPT2 test includes two excitation operators, xxaa and xxxa 
//
//  compile: 
//  g++ -std=c++11 generate-main-caspt2-double.cc -o generate-main-caspt2-double
//  run:
//  ./generate-main-caspt2-double > main.cc 
//
//////////////////////////////////////////////////////////////////////////////////

#ifndef __SRC_INDICES_H
#define __SRC_INDICES_H

#include <iostream>
#include <tuple>
#include <string>
#include <memory>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <sstream>
#include <initializer_list>

class Indices {
  protected:
    std::list<std::string> indices_;

  public:
    Indices(const std::initializer_list<std::string> o) {
      for (auto& i : o) indices_.push_back(i);
    };

    std::string str() const {
      std::stringstream ss;
      for (auto i = indices_.begin(); i != indices_.end(); ++i) {
        if (i != indices_.begin()) ss << ", ";
        ss << "\"" << *i << "\"";
      }
      return ss.str();
    };

    std::string generate_tensor(std::string tag) const {
      std::stringstream ss;
      ss << "  shared_ptr<Op> " << tag << "(new Op(" << str() << "));" << std::endl; 
      return ss.str();
    };

};

#endif
