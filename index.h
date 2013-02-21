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


#ifndef __INDEX_H
#define __INDEX_H

#include <string>
#include <sstream>
#include <cassert>

namespace smith {

/// Used to describe spin orbital attribute. Index uses this. 
class Spin {
  protected:
    /// Spin number.
    int num_;
  public:
    /// Construct spin, initialize number information.
    Spin() : num_(0) {};
    ~Spin() {};

    /// Returns spin.
    int num() const { return num_; };
    /// Sets spin.
    void set_num(const int i) { num_ = i; };
   
    /// Returns spin in parenthesis (number).
    std::string str() const {
      std::stringstream ss;
      ss << "(" << num_ << ")";
      return ss.str();
    };
};

/// A class for orbital attributes. Index defined by label (space), spin, electron number and if is transposed (daggered).
class Index {
  protected:
    /// Index label, related to closed, active, or virtual (c, x, and a, respectively).
    std::string label_;
    /// Index number (electron).
    int num_;
    /// If transposed, ie daggered. Important in Wick's theorem and equations.
    bool dagger_;
    /// Spin of index.
    std::shared_ptr<Spin> spin_;
    /// If target, index will not be summed over in code generation.
    bool target_;
  

  public:
    /// Make index object from label and dagger info. Initialize label, number(0), dagger, and target(false).
    Index(std::string lab, bool dag) : label_(lab), num_(0), dagger_(dag), target_(false){};
    ~Index() {};

    /// Return index number. 
    int num() const { return num_; };
    /// Return if should be transposed.
    bool dagger() const { return dagger_; };
    /// Set index number.
    void set_num(const int i) { num_ = i; };
    /// Return index label (orbital type).
    const std::string label() const { return label_; };
    /// Set index type, default is a (virtual).
    void set_label(const std::string& a) { label_ = a; };

    /// If active.  Checks label if active (x).
    bool active() const { return label_ == "x"; };
  
    /// Return target index status.
    bool target() const { return target_; };
    /// set target as true
    void mark_target() { target_ = true; };
  

    /// Sets spin.
    void set_spin(const std::shared_ptr<Spin> s) { spin_ = s; };
    /// Returns spin.
    std::shared_ptr<Spin> spin() { assert(spin_); return spin_; };
    /// Returns const spin.
    const std::shared_ptr<Spin> spin() const { assert(spin_); return spin_; };

    /// Returns true if spin is same for both indices.
    bool same_spin(const std::shared_ptr<Index>& o) const { return o->spin() == spin(); };
    /// Returns true if index number is same for both indices.
    bool same_num(const std::shared_ptr<Index>& o) const { return o->num() == num(); };

    /// Returns string with index label_, and if argument is true: dagger info (nothing or if daggered, +) and spin info.
    std::string str(const bool& opr = true) const {
      std::stringstream ss;
      ss << label_ << num_;
      if (dagger_ && opr) ss << "+";
      if (opr) {
        assert(spin_);
        ss << spin_->str();
      }
      return ss.str();
    };
    /// Returns string with index label_ and num_. 
    std::string str_gen() const {
      std::stringstream ss;
      ss << label_ << abs(num_);
      return ss.str();
    };
    /// Prints index using str().
    void print() const { std::cout << str() << std::endl; };


    /// Clone Index with label_, num_ and dagger_ info. Note that this does not set spin.
    std::shared_ptr<Index> clone() { 
      std::shared_ptr<Index> out(new Index(label_, dagger_));
      out->set_num(num_);
      return out;
    };

    /// Check if indices are equal by comparing num() and label(). Be careful that this does not check dagger! Should not check, actually. 
    bool identical(std::shared_ptr<Index> o) const {
      return num() == o->num() && label() == o->label();
    };

    /// Gives orbital space name (closed_, virt_, active_) based on index label_.
    std::string generate() const {
      std::string out;
      if (label_ == "c") {
        out = "closed_";
      } else if (label_ == "a") {
        out = "virt_";
      } else if (label_ == "x") {
        out = "active_";
      } else {
        throw std::runtime_error("unkonwn index type in Index::generate()");
      }
      return out;
    };

    /// Gives orbital space range name ([0], [1], [2] for closed, active, and virtual spaces, respectively) based on index label.
    std::string generate_range(const std::string postfix = "") const {
      std::string out = "range" + postfix;
      if (label_ == "c") {
        out += "[0]";
      } else if (label_ == "x") {
        out += "[1]";
      } else if (label_ == "a") {
        out += "[2]";
      } else {
        throw std::runtime_error("unkonwn index type in Index::generate_range()");
      }
      return out;
    };

};

/// Global function to list indices in reverse order.
static std::string list_keys(const std::list<std::shared_ptr<Index>>& index) {
  std::stringstream tt; 
  for (auto iter = index.rbegin(); iter != index.rend(); ++iter) {
    if (iter != index.rbegin()) tt << ", ";
    tt << (*iter)->str_gen() << ".key()";
  }
  return tt.str();
}

}

#endif
