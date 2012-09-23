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

class Spin {
  protected:
    int num_;
  public:
    Spin() : num_(0) {};
    ~Spin() {};

    int num() const { return num_; };
    void set_num(const int i) { num_ = i; };

    std::string str() const {
      std::stringstream ss;
      ss << "(" << num_ << ")";
      return ss.str();
    };
};


class Index {
  protected:
    std::string label_;
    int num_;
    bool dagger_;
    std::shared_ptr<Spin> spin_;

  public:
    Index(std::string lab, bool dag) : label_(lab), num_(0), dagger_(dag) {};
    ~Index() {};

    int num() const { return num_; };
    bool dagger() const { return dagger_; };
    void set_num(const int i) { num_ = i; };
    const std::string label() const { return label_; };
    void set_label(const std::string& a) { label_ = a; };

    bool active() const { return label_ == "x"; };

    void set_spin(const std::shared_ptr<Spin> s) { spin_ = s; };
    std::shared_ptr<Spin> spin() { assert(spin_); return spin_; };
    const std::shared_ptr<Spin> spin() const { assert(spin_); return spin_; };

    bool same_spin(const std::shared_ptr<Index>& o) const { return o->spin() == spin(); };
    bool same_num(const std::shared_ptr<Index>& o) const { return o->num() == num(); };

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
    std::string str_gen() const {
      std::stringstream ss;
      ss << label_ << abs(num_);
      return ss.str();
    };

    std::shared_ptr<Index> clone() { // note that this does not set spin.
      std::shared_ptr<Index> out(new Index(label_, dagger_));
      out->set_num(num_);
      return out;
    };

    void print() const { std::cout << str() << std::endl; };

    // be careful that this does not check dagger! Should not check, actually. 
    bool identical(std::shared_ptr<Index> o) const {
      return num() == o->num() && label() == o->label();
    };

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
};

// global function
static std::string list_keys(const std::list<std::shared_ptr<Index> >& index) {
  std::stringstream tt; 
  for (auto iter = index.rbegin(); iter != index.rend(); ++iter) {
    if (iter != index.rbegin()) tt << ", ";
    tt << (*iter)->str_gen() << ".key()";
  }
  return tt.str();
}

}

#endif
