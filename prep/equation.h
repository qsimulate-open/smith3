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

#ifndef __EQUATION_H
#define __EQUATION_H

#include "diagram.h"
#include <sstream>
#include <memory>
#include <list>
#include <vector>
#include <string>

namespace SMITH3 {
namespace Prep {

class Equation {
  protected:
    std::string method_;

    std::list<std::shared_ptr<Diagram>> diagram_;
    std::string label_;
    double fac_;
    std::string tree_type_;
    std::string tree_name_;
    std::pair<bool, bool> braket_;
    bool ci_derivative_;

  public:
    Equation(const std::string m, const std::string l, const std::initializer_list<std::vector<std::shared_ptr<Tensor>>> in, const double d = 1.0,
             const std::pair<bool,bool> brkt = std::make_pair(false,false))
     : Equation(m, l, in, d, "", brkt) { }

    Equation(const std::string m, const std::string l, const std::initializer_list<std::vector<std::shared_ptr<Tensor>>> in, const std::string scalar,
             const std::pair<bool,bool> brkt = std::make_pair(false,false))
     : Equation(m, l, in, 1.0, scalar, brkt) { }
    // end ctor

    Equation(const std::string m, const std::string l, const std::initializer_list<std::vector<std::shared_ptr<Tensor>>> in, const double d, const std::string scalar,
             const std::pair<bool,bool> bk = std::make_pair(false,false))
     : method_(m), label_(l), fac_(d), tree_type_(""), braket_(bk) {

      if (bk.first == false && bk.second == false) {
        ci_derivative_ = false;
      } else {
        ci_derivative_ = true;
      }


      std::list<int> max;
      for (auto& i : in) max.push_back(i.size());

      std::list<std::list<std::shared_ptr<Tensor>>> out;
      std::list<int> current(in.size(), 0);
      std::list<int> start = current;
      do {
        // set the current vector
        std::list<std::shared_ptr<Tensor>> cc;
        auto inp = in.begin();
        for (auto i = current.begin(); i != current.end(); ++i, ++inp)
          cc.push_back((*inp)[*i]);
        // diagonal cc/aa will be removed from the CASPT2 residual equation for efficiency.
//      if ((method_ != "CASPT2") || label_[0] != 'r' || !(cc.back()->external() && cc.front()->tag() == "proje" && (*(++cc.begin()))->external())) {
        if ((method_ != "CASPT2" && method_ != "RelCASPT2") || label_[0] != 'r' || !(cc.back()->external() && cc.front()->tag() == "proje" && (*(++cc.begin()))->external())) {
          out.push_back(cc);
        }

        // get the next vector
        auto m = max.rbegin();
        for (auto i = current.rbegin(); i != current.rend(); ++i, ++m) {
          if (++*i == *m) {
            *i = 0;
          } else {
            break;
          }
        }
      } while (current != start);

      // construct Diagrams
      int cnt = 0;
      for (auto& i : out) {
        std::stringstream ss; ss << label_ << cnt;
        if (!ci_derivative_) {
          if (d == 1.0) {
            diagram_.push_back(std::make_shared<Diagram>(i, ss.str(), scalar));
          } else {
            diagram_.push_back(std::make_shared<Diagram>(i, ss.str(), d, scalar));
          }
        } else {
          if (d == 1.0) {
            diagram_.push_back(std::make_shared<Diagram>(i, ss.str(), scalar, bk));
          } else {
            diagram_.push_back(std::make_shared<Diagram>(i, ss.str(), d, scalar, bk));
          }
        }
        ++cnt;
      }

    }


    std::string tree_label() const { return "t" + label_; }

    void set_tree_type(std::string ttype, std::string tname = "") { tree_type_=ttype; tree_name_ = tname.empty() ? ttype : tname; }

    void merge(std::shared_ptr<Equation> o) {
      diagram_.insert(diagram_.end(), o->diagram_.begin(), o->diagram_.end());
    }

    std::string generate() const {
      std::stringstream ss;
      for (auto& i : diagram_) ss << i->construct_str();
      for (auto& i : diagram_) ss << i->diagram_str();
      for (auto& i : diagram_) ss << i->equation_str();

      for (auto i = diagram_.begin(); i != diagram_.end(); ++i)
        if (i != diagram_.begin())
          ss << "  " << diagram_.front()->eqn_label() << "->merge(" << (*i)->eqn_label() << ");" << std::endl;
      if (ci_derivative_) ss << "  " << diagram_.front()->eqn_label() << "->absorb_ket();" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->duplicates();" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->active();" << std::endl;
      if (method_ != "CASPT2" && method_ != "RelCASPT2" && method_ != "MSCASPT2") {
        ss << "  " << diagram_.front()->eqn_label() << "->reorder_tensors();" << std::endl;
        ss << "  " << diagram_.front()->eqn_label() << "->simplify();" << std::endl;
      }

      if (!tree_type_.empty() && tree_type_ == "residual") {
          ss << "  auto " << tree_label() << " = make_shared<Residual>(e" << diagram_.front()->label() << ", \"" << tree_name_ << "\");" << std::endl;
      } else if (!tree_type_.empty() && tree_type_ == "energy") {
          ss << "  auto " << tree_label() << " = make_shared<Energy>(e" << diagram_.front()->label() << ", \"" << tree_name_ << "\");" << std::endl;
      } else {
          throw std::logic_error("prep/equation.cc error, tree must be of derived type");
      }

      ss << std::endl;
      return ss.str();
    };


};

}}


#endif
