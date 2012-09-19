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
string RDM::generate(string indent) const {
  stringstream tt;
  //first let's do the generate_get_block...
  //problem: below lbl needs to be the rdm name, but isn't this bad b/c label is above in hierarchy
  tt << indent << "std::unique_ptr<double[]> " << "i0" << "data = "
                  << "lbl" << "->" << "get" << "_block(" << "i0" << "hash);" << endl;
 
 
   

  return tt.str();
}




