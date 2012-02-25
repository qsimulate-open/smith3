//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: cost.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef _smith_cost_h
#define _smith_cost_h

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <tuple>
#include <cassert>
#include <iostream>
#include "indexmap.h"

class PCost {

  protected:
    // a vector of integers with length cat
    std::vector<int> pcost_;
    // mapping information. Maybe this is not the most beautiful way..
    IndexMap indmap_;

  public:
    PCost(const std::vector<int>& pcst): pcost_(pcst) { };
    PCost() { pcost_.resize(indmap_.size()); };
    ~PCost() { };

    bool operator<(const PCost& other)  const { return pcost_total() < other.pcost_total(); };
    bool operator>(const PCost& other)  const { return pcost_total() > other.pcost_total(); };
    bool operator==(const PCost& other) const { return pcost() == other.pcost(); };
    bool operator!=(const PCost& other) const { return !(*this == other);};

    const double pcost_total() const {
      double out = 0.0;
      auto j = indmap_.begin();
      for (auto i = pcost_.begin(); i != pcost_.end(); ++i, ++j)
        out += ::log(static_cast<double>(j->second.second))* *i;
      return out;
    }
    const std::vector<int> pcost() const { return pcost_;};

    void add(std::vector<int>& o) {
      for (auto i = pcost_.begin(), j = o.begin(); i != pcost_.end(); ++i, ++j) *i += *j; 
    };

    int pcost(const int i) const { return pcost_[i]; };

    const std::string show() const;

};


class Cost {

  protected:
    std::vector<PCost> cost_;

  public:
    Cost(const std::vector<PCost>& cst) : cost_(cst) { };
    Cost() { };
    ~Cost() {};

    bool operator<(const Cost& other) const {
      std::vector<PCost> otherc = other.cost();
      std::vector<PCost> myc = cost();
      for (auto i = myc.begin(), j = otherc.begin(); i != myc.end(); ++i, ++j) {
        if (j == otherc.end()) return false;
        if      (*i < *j)      return true;
        else if (*i > *j)      return false;
      }
      return true;
    };

    bool operator==(const Cost& other) const { return cost()==other.cost(); };
    bool operator!=(const Cost& other) const { return !(*this == other);};
    bool operator>(const Cost& other)  const { return !(*this < other);};

    const std::vector<PCost> cost() const {return cost_;};

    void add_pcost(const PCost& p) { cost_.push_back(p); };
//  void add_pcost(int i, int j, int k) { PCost a(i, j, k); cost_.push_back(a); };

    const std::string show() const;

    void sort_pcost();

};


#endif

