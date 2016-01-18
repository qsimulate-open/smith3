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
  ss << "// BAGEL - Brilliantly Advanced General Electronic Structure Library" << std::endl;
  ss << "// Filename: " << filename << std::endl;
  ss << "// Copyright (C) 2014 Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// Author: Shiozaki group <shiozaki@northwestern.edu>" << std::endl;
  ss << "// Maintainer: Shiozaki group" << std::endl;
  ss << "//" << std::endl;
  ss << "// This file is part of the BAGEL package." << std::endl;
  ss << "//" << std::endl;
  ss << "// This program is free software; you can redistribute it and/or modify" << std::endl;
  ss << "// it under the terms of the GNU General Public License as published by" << std::endl;
  ss << "// the Free Software Foundation, either version 3 of the License, or" << std::endl;
  ss << "// (at your option) any later version." << std::endl;
  ss << "//" << std::endl;
  ss << "// This program is distributed in the hope that it will be useful," << std::endl;
  ss << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
  ss << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
  ss << "// GNU General Public License for more details." << std::endl;
  ss << "//" << std::endl;
  ss << "// You should have received a copy of the GNU General Public License" << std::endl;
  ss << "// along with this program.  If not, see <http://www.gnu.org/licenses/>." << std::endl;
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

bool same_tensor__(std::string a, std::string b) {
  auto strip = [](const std::string& label) {
    size_t found = label.find("dagger");
    if (found != std::string::npos) {
      std::string tmp(label.begin(), label.begin() + found);
      return tmp;
    } else {
      return label;
    }
  };
  return strip(a) == strip(b);
}

int count_distinct_tensors__(const std::vector<std::string>& labels) {
  int out = 0;
  std::vector<std::string> done;
  for (auto& i : labels) {
    bool found = false;
    for (auto& j : done)
      if (same_tensor__(i,j))
        found = true;
    if (!found) {
      ++out;
      done.push_back(i);
    }
  }
  return out;
}

static const std::string DataType = "double";
//static const std::string DataType = "std::complex<double>";
static const double fac2 = (DataType == "double" ? 2.0 : 1.0);
static const std::string GEMM = (DataType == "double" ? "dgemm_" : "zgemm3m_");
static const std::string SCAL = (DataType == "double" ? "dscal_" : "zscal_");
static const std::string DOT = (DataType == "double" ? "ddot_" : "zdotu_");
static const std::string MatType = (DataType == "double" ? "Matrix" : "ZMatrix");

}
}

#endif

