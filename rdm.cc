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
using namespace smith;

// TODO replace by a standard function as soon as I get this working
namespace smith {
  static string prefac__(const double& factor_) {
    stringstream tt;
    if (fabs(factor_-1.0) < 1.0e-10) { tt << "1,1";
    } else if (fabs(factor_+1.0) < 1.0e-10) { tt << "-1,1";
    } else if (fabs(factor_-2.0) < 1.0e-10) { tt << "2,1";
    } else if (fabs(factor_+2.0) < 1.0e-10) { tt << "-2,1";
    } else if (fabs(factor_-4.0) < 1.0e-10) { tt << "4,1";
    } else if (fabs(factor_+4.0) < 1.0e-10) { tt << "-4,1";
    } else if (fabs(factor_-8.0) < 1.0e-10) { tt << "8,1";
    } else if (fabs(factor_+8.0) < 1.0e-10) { tt << "-8,1";
    } else if (fabs(factor_-16.0) < 1.0e-10) { tt << "16,1";
    } else if (fabs(factor_+16.0) < 1.0e-10) { tt << "-16,1";
    } else if (fabs(factor_-32.0) < 1.0e-10) { tt << "32,1";
    } else if (fabs(factor_+32.0) < 1.0e-10) { tt << "-32,1";
    } else if (fabs(factor_-0.5) < 1.0e-10) { tt << "1,2";
    } else if (fabs(factor_+0.5) < 1.0e-10) { tt << "-1,2";
    } else if (fabs(factor_-0.25) < 1.0e-10) { tt << "1,4";
    } else if (fabs(factor_+0.25) < 1.0e-10) { tt << "-1,4";
    } else {
      tt << "this case is not yet considered " << factor_ << " in RDM::generate()";
      throw runtime_error(tt.str());
    }
    return tt.str();
  }
}


// do a special sort_indices for rdm summation with possible delta cases
string RDM::generate(string indent, const string tlab, const list<shared_ptr<Index> >& loop) const {
  assert(!index_.empty() && !loop.empty());
  stringstream tt;
  indent += "  ";
  const string itag = "i";

  // now do the sort
  vector<string> close;

  // in case delta_ is not empty
  if (!delta_.empty()) {
    // first delta loops for blocks
    tt << indent << "if (";
    for (auto d = delta_.begin(); d != delta_.end(); ++d) {
      tt << d->first->str_gen() << " == " << d->second->str_gen() << (d != --delta_.end() ? " && " : "");
    }
    tt << ") {" << endl;
    close.push_back(indent + "}");
    indent += "  ";

    tt << indent << "std::vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
    tt << indent << "std::unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;

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
      const string tmp = "+" + (*ri)->str_gen() + ".size()*(";
      tt << itag << inum << (ri != --loop.rend() ? tmp : "");
    }
    for (auto ri = ++loop.begin(); ri != loop.end(); ++ri)
      tt << ")";
    tt << "]" << endl;

    // make data part of summation
    tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * data[";
    for (auto riter = index_.rbegin(); riter != index_.rend(); ++riter) {
      const string tmp = "+" + (*riter)->str_gen() + ".size()*(";
      tt << itag << (*riter)->num() << (riter != --index_.rend() ? tmp : "");
    }
    for (auto riter = ++index_.begin(); riter != index_.end(); ++riter)
      tt << ")";
    tt << "];" << endl;

    // close loops
    for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
      tt << *iter << endl;

  // if delta_ is empty call sort_indices
  } else {
    // loop up the operator generators
    tt << indent << "std::vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
    tt << indent << "std::unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;
   
    // call sort_indices here
    vector<int> done;
    tt << indent << "sort_indices<";
    for (auto i = loop.rbegin(); i != loop.rend(); ++i) {
      int cnt = 0;
      for (auto j = index_.rbegin(); j != index_.rend(); ++j, ++cnt) {
        if ((*i)->identical(*j)) break;
      }
      if (cnt == index_.size()) throw logic_error("should not happen.. RDM::generate");
      done.push_back(cnt);
    }
    // then fill out others
    for (int i = 0; i != index_.size(); ++i) {
      if (find(done.begin(), done.end(), i) == done.end())
        done.push_back(i);
    }
    // write out
    for (auto& i : done) 
      tt << i << ",";

    // add factor information
    tt << 0 << ",1," << prefac__(fac_);
   
    // add source data dimensions
    tt << ">(data, " << tlab << "data, " ;
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
        tt << (*iter)->str_gen() << "->size()";
    }
    tt << ");" << endl;

  }

  return tt.str();
}




