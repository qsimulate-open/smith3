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

#ifndef __SRC_TENSOR_H
#define __SRC_TENSOR_H

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

namespace SMITH3 {
namespace Prep {

class Tensor {
  protected:
    std::string tag_;
    std::string base_;
    std::list<std::string> indices_;

  public:
    Tensor(const std::string b, const std::string postfix, const std::initializer_list<std::string> o)
     : tag_(b+postfix), base_(b) {
      for (auto& i : o) indices_.push_back(i);
    };

    std::string str_index() const {
      std::stringstream ss;
      for (auto i = indices_.begin(); i != indices_.end(); ++i) {
        if (i != indices_.begin()) ss << ", ";
        ss << "\"" << *i << "\"";
      }
      return ss.str();
    };

    std::string generate() const {
      std::stringstream ss;
      ss << "  shared_ptr<Op> " << tag_ << "(new Op(\"" << base_ << "\"" << (indices_.empty() ? "" : ", ") << str_index() << "));" << std::endl; 
      return ss.str();
    };

    std::string base() const { return base_; };
    std::string tag() const { return tag_; };
};

}}

#endif
