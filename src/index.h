//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: index.h
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


#ifndef __INDEX_H
#define __INDEX_H

#include <string>
#include <sstream>
#include <memory>
#include <list>
#include <iostream>
#include <cassert>

namespace smith {

/// Used to describe spin orbital attribute. Index uses this.
class Spin {
  protected:
    /// Spin number.
    int num_;
    /// Partial spin
    bool alpha_;
  public:
    /// Construct spin, initialize number information.
    Spin(const bool al) : num_(0), alpha_(al) { }
    ~Spin() { }

    /// Returns spin.
    int num() const { return num_; }
    /// Sets spin.
    void set_num(const int i) { num_ = i; }
    /// Returns alpha.
    bool alpha() const { return alpha_; }
    /// Sets alpha
    void set_alpha(const bool a) { alpha_ = a; }

    /// Returns spin in parenthesis (number).
    std::string str() const {
      std::stringstream ss;
      ss << "(" << num_ << (alpha_ ? "*" : "") << ")";
      return ss.str();
    }
};

/// A class for tensor indices. Can refer to orbital attributes: Index defined by label (space), spin, electron number and if is transposed (daggered). Also can refer to cI index.
class Index_Core {
  protected:
    /// Index label, related to closed, active, or virtual (c, x, and a, respectively).
    std::string label_;
    /// Index number (if orbital index, electron).
    int num_;
    /// If transposed, ie daggered. Important in Wick's theorem and equations.
    bool dagger_;

  public:
    /// Make index object from label and dagger info. Initialize label, number(0), dagger.
    Index_Core(std::string lab, bool dag) : label_(lab), num_(0), dagger_(dag) {}
    /// Make a copy of the index
    Index_Core(const Index_Core& o) : label_(o.label_), num_(o.num_), dagger_(o.dagger_) { }
    /// Make copy of the index but with reversed dagger info
    Index_Core(const Index_Core& o, bool b) : label_(o.label_), num_(o.num_), dagger_(!b) { }
    /// Make copy of index but with altered number
    Index_Core(const Index_Core& o, int i) : label_(o.label_), num_(i), dagger_(o.dagger_) { }

    /// Return index number.
    int num() const { return num_; }
    /// Return if should be transposed.
    bool dagger() const { return dagger_; }
    /// Set index number.
    void set_num(const int i) { num_ = i; }
    /// Return index label (orbital type).
    const std::string label() const { return label_; }
    /// Set index type, default is a (virtual).
    void set_label(const std::string& a) { label_ = a; }
};

class Index {

  protected:
    /// Index
    std::shared_ptr<Index_Core> core_;
    /// Spin of index.
    mutable std::shared_ptr<Spin> spin_; // TODO mutable should be removed

  public:
    /// Make index object from label and dagger info. Initialize label, number(0), dagger.
    Index(std::string lab, bool dag) { core_ = std::make_shared<Index_Core>(lab, dag); }
    Index(const Index& o) : spin_(o.spin_) { core_ = std::make_shared<Index_Core>(*o.core_); }
    /// Make copy of the index but with reversed dagger info
    Index(const Index& o, bool b) : spin_(o.spin_) { core_ = std::make_shared<Index_Core>(*o.core_, b); }
    /// Make copy of index but with altered number
    Index(const Index& o, int i) : spin_(o.spin_) { core_ = std::make_shared<Index_Core>(*o.core_, i); }
    Index(std::shared_ptr<Index_Core> c) : core_(c) { }
    ~Index() { }

    /// Return index number.
    int num() const { return core_->num(); }
    /// Return if should be transposed.
    bool dagger() const { return core_->dagger(); }
    /// Set index number.
    void set_num(const int i) { core_->set_num(i); }
    /// Return index label (orbital type).
    std::string label() const { return core_->label(); }
    /// Set index type, default is a (virtual).
    void set_label(const std::string& a) { core_->set_label(a); }

    /// If active.  Checks label if active (x).
    bool active() const { return label() == "x"; }

    /// Returns true if index number is same for both indices.
    bool same_num(const std::shared_ptr<const Index>& o) const { return o->num() == num(); }
    /// Returns true if label is same for both indices.
    bool same_label(const std::shared_ptr<const Index>& o) const { return o->label() == label(); }

    /// Returns string with index label_, and if argument is true: dagger info (nothing or if daggered, +) and spin info.
    std::string str(const bool opr = true) const {
      std::stringstream ss;
      ss << label() << num();
      if (dagger() && opr) ss << "+";
      if (opr && spin_)    ss << spin_->str();
      if (spin_ && spin_->alpha()) ss << "*";
      return ss.str();
    }
    /// Returns string with index label_ and num_.
    std::string str_gen() const {
      std::stringstream ss;
      ss << label() << abs(num());
      return ss.str();
    }
    /// Prints index using str().
    void print() const { std::cout << str() << std::endl; };
    /// Sets spin.

    void set_spin(const std::shared_ptr<Spin> s) const { spin_ = s; }
    /// Returns spin.
    std::shared_ptr<Spin> spin() { assert(spin_); return spin_; }
    /// Returns const spin.
    const std::shared_ptr<Spin> spin() const { assert(spin_); return spin_; }

    /// Returns true if spin is same for both indices.
    bool same_spin(const std::shared_ptr<const Index>& o) const { return o->spin() == spin(); }

    /// Clone Index with label_, num_ and dagger_ info. Note that this does not set spin.
    std::shared_ptr<Index> clone() const {
      return std::make_shared<Index>(core_);
    }

    /// Check if indices are equal by comparing num() and label(). Be careful that this does not check dagger! Should not check, actually.
    bool identical(std::shared_ptr<const Index> o) const {
      return num() == o->num() && label() == o->label();
    }

    /// Gives orbital space name (closed_, virt_, active_) based on index label_.
    std::string generate() const {
      std::string out;
      if (label() == "c") {
        out = "closed_";
      } else if (label() == "a") {
        out = "virt_";
      } else if (label() == "x") {
        out = "active_";
      } else if (label() == "ci") {
        out = "ci_";
      } else {
        throw std::runtime_error("unkonwn index type in Index::generate()");
      }
      return out;
    }

    /// Gives index range name ([0], [1], [2] for closed, active, virtual orbital spaces, respectively and [3] for ci range) based on index label.
    std::string generate_range(const std::string postfix = "") const {
      std::string out = "range" + postfix;
      if (label() == "c") {
        out += "[0]";
      } else if (label() == "x") {
        out += "[1]";
      } else if (label() == "a") {
        out += "[2]";
      } else if (label() == "ci") {
        out += "[3]";
      } else {
        throw std::runtime_error("unkonwn index type in Index::generate_range()");
      }
      return out;
    }

};

}

#endif
