//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: constants.h
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


#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cassert>

namespace smith {

static std::string header(const std::string& filename) {
  std::stringstream ss;
  ss << "//" << std::endl;
  ss << "// Newint - Parallel electron correlation program." << std::endl;
  ss << "// Filename: " << filename << ".h" << std::endl;
  ss << "// Copyright (C) 2012 Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// Author: Shiozaki group <shiozaki@northwestern.edu>" << std::endl;
  ss << "// Maintainer: Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// This file is part of the BAGEL package." << std::endl;
  ss << "//" << std::endl;
  ss << "// The Newint package is free software; you can redistribute it and/or modify" << std::endl;
  ss << "// it under the terms of the GNU Library General Public License as published by" << std::endl;
  ss << "// the Free Software Foundation; either version 3, or (at your option)" << std::endl;
  ss << "// any later version." << std::endl;
  ss << "//" << std::endl;
  ss << "// The Newint package is distributed in the hope that it will be useful," << std::endl;
  ss << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
  ss << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
  ss << "// GNU Library General Public License for more details." << std::endl;
  ss << "//" << std::endl;
  ss << "// You should have received a copy of the GNU Library General Public License" << std::endl;
  ss << "// along with the Newint package; see COPYING.  If not, write to" << std::endl;
  ss << "// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA." << std::endl;
  ss << "//" << std::endl;
  ss << "" << std::endl;
  ss << "" << std::endl;
  return ss.str();
};

static std::string prefac__(const double& factor_) {
  const double thresh = 1.0e-10;
  // bruteforce way...
  const int large = 1024;
  int i, j;
  for (i = 1; i != large; ++i) {
    j = static_cast<int>(factor_*i + (factor_<0 ? -thresh : thresh));
    if (fabs(factor_*i-j) < thresh)
      break;
  }
  std::stringstream ss;
  ss << j << "," << i;
  return ss.str();
};


}

#endif

