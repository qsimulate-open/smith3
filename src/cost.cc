//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: cost.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <algorithm>
#include <sstream>
#include "cost.h"

using namespace std;
using namespace smith;

string PCost::show() const {
  stringstream out;
  auto j = pcost_.begin();
  for (auto i = indmap_.begin(); i != indmap_.end(); ++i, ++j)
    out << i->first <<  *j;
  return out.str();
}


void Cost::sort_pcost() {
  sort(cost_.rbegin(), cost_.rend());
}


string Cost::show() const {
  string out;
  for (auto& i : cost_) out += i.show() + " ";
  return out;
}


