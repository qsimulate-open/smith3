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
  mm << "#include <fstream>" << std::endl;
  mm << "#include \"forest.h\"" << std::endl;
  mm << "#include \"residual.h\"" << std::endl;
  mm << "#include \"energy.h\"" << std::endl;
  mm << "" << std::endl;
  mm << "using namespace std;" << std::endl;
  mm << "using namespace smith;" << std::endl;
  mm << "" << std::endl;
  mm << "int main() {" << std::endl;
  return mm.str();
};


static std::string footer(const std::string res, const std::string energy, const std::string correction, const std::string density, const std::string density1, const std::string density2, const std::string dedci) {
  std::stringstream mm;

  mm << "  list<shared_ptr<Tree>> trees = {" << res << ", " << energy << ", " << correction << ", " << density << ", " <<  density1 << ", " << density2 << ", " << dedci <<  "};" << std::endl;
  mm << "  auto fr = make_shared<Forest>(trees);" << std::endl;

  mm << "" <<  std::endl;
  mm << "  fr->filter_gamma();" << std::endl;
  mm << "  list<shared_ptr<Tensor>> gm = fr->gamma();" << std::endl;
  mm << "  const list<shared_ptr<Tensor>> gamma = gm;" << std::endl;

  mm << "" <<  std::endl;
  mm << "  auto tmp = fr->generate_code();" << std::endl;

  mm << "" <<  std::endl;
  mm << "  ofstream fs(fr->name() + \".h\");" << std::endl;
  mm << "  ofstream es(fr->name() + \"_tasks.h\");" << std::endl;
  mm << "  ofstream cs(fr->name() + \"_gen.cc\");" << std::endl;
  mm << "  ofstream ds(fr->name() + \"_tasks.cc\");" << std::endl;
  mm << "  ofstream gs(fr->name() + \".cc\");" << std::endl;
  mm << "  fs << tmp.ss.str();" << std::endl;
  mm << "  es << tmp.tt.str();" << std::endl;
  mm << "  cs << tmp.cc.str();" << std::endl;
  mm << "  ds << tmp.dd.str();" << std::endl;
  mm << "  gs << tmp.ee.str();" << std::endl;
  mm << "  fs.close();" << std::endl;
  mm << "  es.close();" << std::endl;
  mm << "  cs.close();" << std::endl;
  mm << "  ds.close();" << std::endl;
  mm << "  gs.close();" << std::endl;
  mm << "  cout << std::endl;" << std::endl;
  mm << "" <<  std::endl;
  mm << "  // output" << std::endl;
  mm << "  cout << std::endl << \"   ***  Residual  ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << res << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Energy E2 ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << energy << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Correlated norm <1|1> ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << correction << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Correlated one-body density matrix d2 ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << density << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Correlated one-body density matrix d1 ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << density1 << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  Correlated two-body density matrix D1 ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << density2 << "->print();" << std::endl;
  mm << "  cout << std::endl << \"   ***  CI derivative  ***\" << std::endl << std::endl;" << std::endl;
  mm << "  " << dedci << "->print();" << std::endl;
  mm << "  cout << std::endl << std::endl;" << std::endl;
  mm << "" <<  std::endl;
  mm << "  return 0;" << std::endl;
  mm << "}" << std::endl;
  return mm.str();
};

}}

#endif
