#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <sstream>
#include <string>

namespace SMITH3 {
namespace Prep {

static std::string header() {
  std::stringstream mm;
  mm << "//" << std::endl;
  mm << "// SMITH3 - generates spin-free multireference electron correlation programs." << std::endl;
  mm << "// Filename: main.cc" << std::endl;
  mm << "// Copyright (C) 2012 Toru Shiozaki" << std::endl;
  mm << "//" << std::endl;
  mm << "// Author: Toru Shiozaki <shiozaki@northwestern.edu>" << std::endl;
  mm << "// Maintainer: Shiozaki group" << std::endl;
  mm << "//" << std::endl;
  mm << "// This file is part of the SMITH3 package." << std::endl;
  mm << "//" << std::endl;
  mm << "// The SMITH3 package is free software; you can redistribute it and/or modify" << std::endl;
  mm << "// it under the terms of the GNU Library General Public License as published by" << std::endl;
  mm << "// the Free Software Foundation; either version 2, or (at your option)" << std::endl;
  mm << "// any later version." << std::endl;
  mm << "//" << std::endl;
  mm << "// The SMITH3 package is distributed in the hope that it will be useful," << std::endl;
  mm << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
  mm << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
  mm << "// GNU Library General Public License for more details." << std::endl;
  mm << "//" << std::endl;
  mm << "// You should have received a copy of the GNU Library General Public License" << std::endl;
  mm << "// along with the SMITH3 package; see COPYING.  If not, write to" << std::endl;
  mm << "// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA." << std::endl;
  mm << "//" << std::endl;
  mm << "" << std::endl;
  mm << "" << std::endl;
  mm << "// This program is supposed to perform Wick's theorem for multireference problems." << std::endl;
  mm << "// Spin averaged quantities assumed." << std::endl;
  mm << "" << std::endl;
  mm << "#include <iostream>" << std::endl;
  mm << "#include <fstream>" << std::endl;
  mm << "#include <list>" << std::endl;
  mm << "#include <string>" << std::endl;
  mm << "#include \"equation.h\"" << std::endl;
  mm << "#include \"tree.h\"" << std::endl;
  mm << "" << std::endl;
  mm << "using namespace std;" << std::endl;
  mm << "using namespace smith;" << std::endl;
  mm << "" << std::endl;
  mm << "int main() {" << std::endl;
  return mm.str();
};


static std::string footer(const std::string res, const std::string energy) {
  std::stringstream mm;
  mm << "" <<  std::endl; 
  mm << "  ofstream fs(" << res << "->tree_name() + \".h\");" << std::endl;
  mm << "  ofstream es(" << res << "->tree_name() + \"_tasks.h\");" << std::endl;
  mm << "  pair<string, std::string> tmp = " << res << "->generate_task_list(false, " << energy << ");" << std::endl;
  mm << "  fs << tmp.first;" << std::endl;
  mm << "  es << tmp.second;" << std::endl;
  mm << "  fs.close();" << std::endl;
  mm << "  es.close();" << std::endl;
  mm << "  cout << std::endl;" << std::endl;
  mm << "" <<  std::endl; 
  mm << "  // output" << std::endl;
  mm << "  cout << std::endl << \"   *** Residual ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << res << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Energy  ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << energy << "->print();" << std::endl;
  mm << "  cout << std::endl << std::endl;" << std::endl;
  mm << "" <<  std::endl; 
  mm << "  return 0;" << std::endl;
  mm << "}" << std::endl;
  return mm.str();
};

}}

#endif
