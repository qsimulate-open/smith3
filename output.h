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
// The SMITH3 package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
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

  OutStream() { }
  OutStream(const OutStream& o) {
    ss << o.ss.str();
    tt << o.tt.str();
    cc << o.cc.str();
    dd << o.dd.str();
    ee << o.ee.str();
  }
  OutStream& operator=(const OutStream& a) {
    ss.str(std::string()); tt.str(std::string()); cc.str(std::string()); dd.str(std::string()); ee.str(std::string());
    ss.clear(); tt.clear(); cc.clear(); dd.clear(); ee.clear();
    ss << a.ss.str(); tt << a.tt.str(); cc << a.cc.str(); dd << a.dd.str(); ee << a.ee.str();
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
  return o;
}
}

} 

#endif
