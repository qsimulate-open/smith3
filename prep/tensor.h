//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tensor.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
    bool alpha_;

  public:
    Tensor(const std::string b, const std::string postfix, const std::initializer_list<std::string> o)
     : tag_(b+postfix), base_(b), alpha_(false) {
      for (auto& i : o) indices_.push_back(i);
    }

    Tensor(const std::string postfix, const std::initializer_list<std::string> o, const bool a = false)
     : tag_("ex_" + postfix), base_(""), alpha_(a) {
      for (auto& i : o) indices_.push_back(i);
    }

    std::string str_index() const {
      std::stringstream ss;
      for (auto i = indices_.begin(); i != indices_.end(); ++i) {
        if (i != indices_.begin()) ss << ", ";
        std::string j = *i;
        std::transform(j.begin(), j.end(), j.begin(), ::toupper);
        ss << "_" << j;
      }
      return ss.str();
    }

    std::string generate() const {
      std::stringstream ss;
      std::string al = alpha_ ? ", true" : "";
      if (!base_.empty()) {
        ss << "  shared_ptr<Operator> " << tag_ << " = make_shared<Op>(\"" << base_ << "\"" << (indices_.empty() ? "" : ", ") << str_index() << al << ");" << std::endl;
      } else {
        ss << "  shared_ptr<Operator> " << tag_ << " = make_shared<Op>(" << str_index() << al << ");" << std::endl;
      }
      return ss.str();
    }

    std::string base() const { return base_; }
    std::string tag() const { return tag_; }

    bool external() const {
      std::list<std::string> l1 = {"a", "a", "c", "c"};
      std::list<std::string> l2 = {"c", "c", "a", "a"};
      return l1 == indices_ || l2 == indices_;
    }
};

}}

#endif
