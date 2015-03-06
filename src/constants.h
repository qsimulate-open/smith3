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
#include <algorithm>

namespace smith {
namespace {

std::string header(const std::string& filename) {
  std::stringstream ss;
  ss << "//" << std::endl;
  ss << "// BAGEL - Parallel electron correlation program." << std::endl;
  ss << "// Filename: " << filename << std::endl;
  ss << "// Copyright (C) 2014 Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// Author: Shiozaki group <shiozaki@northwestern.edu>" << std::endl;
  ss << "// Maintainer: Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// This file is part of the BAGEL package." << std::endl;
  ss << "//" << std::endl;
  ss << "// The BAGEL package is free software; you can redistribute it and/or modify" << std::endl;
  ss << "// it under the terms of the GNU Library General Public License as published by" << std::endl;
  ss << "// the Free Software Foundation; either version 3, or (at your option)" << std::endl;
  ss << "// any later version." << std::endl;
  ss << "//" << std::endl;
  ss << "// The BAGEL package is distributed in the hope that it will be useful," << std::endl;
  ss << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
  ss << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
  ss << "// GNU Library General Public License for more details." << std::endl;
  ss << "//" << std::endl;
  ss << "// You should have received a copy of the GNU Library General Public License" << std::endl;
  ss << "// along with the BAGEL package; see COPYING.  If not, write to" << std::endl;
  ss << "// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA." << std::endl;
  ss << "//" << std::endl;
  ss << "" << std::endl;
  ss << "" << std::endl;
  return ss.str();
}

std::string prefac__(const double& factor_) {
  const double thresh = 1.0e-10;
  // bruteforce way...
  const int large = 1024;
  int i;
  for (i = 1; i != large; ++i) {
    if (fabs(factor_*i-std::round(factor_*i)) < thresh)
      break;
  }
  std::stringstream ss;
  ss << std::round(factor_*i) << "," << i;
  return ss.str();
}

std::string target_name__(std::string label) {
  std::string out;
  if (label == "residual")      out = "r";
  else if (label == "source")   out = "s";
  else if (label == "density")  out = "den2";
  else if (label == "density1") out = "den1";
  else if (label == "density2") out = "Den1";
  else if (label == "deci")     out = "deci";
  else if (label == "norm")     out = "n";
  else throw std::logic_error("unrecognized label in constant.h static string target_name");
  return out;
}

std::string merge__(std::vector<std::string> array, std::string name = "") {
  std::stringstream ss;
  std::vector<std::string> done;
  for (auto& label : array) {
    size_t found = label.find("dagger");
    if (found != std::string::npos) {
      std::string tmp(label.begin(), label.begin() + found);
      label = tmp;
    }
    // we only register once
    if (std::find(done.begin(), done.end(), label) != done.end()) continue;
    done.push_back(label);

    // some tweaks
    if (label == "f1" || label == "v2" || label == "h1")
      label = label + "_";
    else if (label != array.front() && label.find("Gamma") != std::string::npos)
      label = label + "_()";

    ss << (label != array.front() ? ", " : "") << ((label == "proj") ? target_name__(name) : label);
  }
  return ss.str();
}

}
}

#endif

