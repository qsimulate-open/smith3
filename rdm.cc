//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: active_gen.cc
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


#include "constants.h"
#include "active.h"
#include <sstream>
#include <algorithm>
#include <iomanip>

using namespace std;


// need to do here both what generate_get_block and a generate_sort_indices usually does
// ok to do this here b/c these blocks need special attention because of deltas
string RDM::generate(string indent, const string tlab, const list<shared_ptr<Index> >& loop) const {
  assert(!index_.empty() && !loop.empty());
  stringstream tt;
  indent += "  ";
  const string itag = "i";

  // first let's do the generate_get_block for an rdm..
  tt << indent << "std::vector<size_t> i0hash = {";
  for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
    if (iter != index_.rbegin()) tt << ", ";
    tt << (*iter)->str_gen() << ".key()";
  }
  tt << "};" << endl;
  tt << indent << "std::unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;

  // now do the sort
  vector<string> close;

  // in case delta_ is not empty
  if (!delta_.empty()) {
    // first delta loops for blocks
    tt << indent << "if (";
    for (auto d = delta_.begin(); d != delta_.end(); ++d) {
      tt << d->first->str_gen() << " == " << d->second->str_gen() << (d != --delta_.end() ? "&&" : "");
    }
    tt << ") {" << endl;
    close.push_back(indent + "}");
    indent += "  ";

    // start sort loops
    for (auto& i : loop) {
      const int inum = i->num();
      bool found = false;
      for (auto& d : delta_)
        if (d.first->num() == inum) found = true;
      if (!found) { 
        tt << indent << "for (int " << itag << inum << " = 0; " << itag << inum << " != " << i->str_gen() << ".size(); ++" << itag << inum << ") {" << endl;
        close.push_back(indent + "}");
        indent += "  ";
      }
    }

    // make odata part of summation for target
    tt  << indent << "odata[";
    for (auto ri = loop.rbegin(); ri != loop.rend(); ++ri) {
      int inum = (*ri)->num();
      for (auto& d : delta_)
        if (d.first->num() == inum) inum = d.second->num();
      tt << itag << inum << "+" << (*ri)->str_gen() << ".size()" << (ri != --loop.rend() ? "*(" : "");
    }
    for (auto ri = ++loop.begin(); ri != loop.end(); ++ri)
      tt << ")";
    tt << "]" << endl;

    // make data part of summation
    tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * data[";
    for (auto riter = index_.rbegin(); riter != index_.rend(); ++riter)
      tt << itag << (*riter)->num() << "+" << (*riter)->str_gen() << ".size()" << (riter != --index_.rend() ? "*(" : "");
    for (auto riter = ++index_.begin(); riter != index_.end(); ++riter)
      tt << ")";
    tt << "];" << endl;

    // close loops
    for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
      tt << *iter << endl;

  // if delta_ is empty call sort_indices
  } else {
    // TODO please code appropriate code here (loop up the operator generators)
    tt << "****** please write a code in rdm.cc" << endl;
  }
  return tt.str();
}




