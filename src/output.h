//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: output.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SMITH_OUTPUT_H
#define __SMITH_OUTPUT_H

#include <sstream>

namespace smith {

struct OutStream {
  std::stringstream ss; //name.h
  std::stringstream tt; //name_tasks.h
  std::stringstream cc; //name_gen.cc
  std::stringstream dd; //name_tasks.cc
  std::stringstream ee; //name.cc
  std::stringstream gg; //name_gamma.cc

  OutStream() { }
  OutStream(const OutStream& o) {
    ss << o.ss.str();
    tt << o.tt.str();
    cc << o.cc.str();
    dd << o.dd.str();
    ee << o.ee.str();
    gg << o.gg.str();
  }
  OutStream& operator=(const OutStream& a) {
    ss.str(std::string()); tt.str(std::string()); cc.str(std::string()); dd.str(std::string()); ee.str(std::string()); gg.str(std::string());
    ss.clear(); tt.clear(); cc.clear(); dd.clear(); ee.clear(); gg.clear();
    ss << a.ss.str(); tt << a.tt.str(); cc << a.cc.str(); dd << a.dd.str(); ee << a.ee.str(); gg << a.gg.str();
    return *this;
  }
};

namespace {
OutStream& operator<<(OutStream& o, const OutStream& a) {
  o.ss << a.ss.str();
  o.tt << a.tt.str();
  o.cc << a.cc.str();
  o.dd << a.dd.str();
  o.ee << a.ee.str();
  o.gg << a.gg.str();
  return o;
}
}

}

#endif
