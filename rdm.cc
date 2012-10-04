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


// do a special sort_indices for rdm summation with possible delta cases
string RDM::generate(string indent, const string tlab, const list<shared_ptr<Index> >& loop) const {
// temp turn off check for testing task file
  //if (!index_.empty() && !loop.empty()) {
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

     tt << indent << "vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
     tt << indent << "unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;

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
     tt << indent << "vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
     tt << indent << "unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;
   
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
     tt << "1,1," << prefac__(fac_);
   
     // add source data dimensions
     tt << ">(data, " << tlab << "data, " ;
     for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
       if (iter != index_.rbegin()) tt << ", ";
         tt << (*iter)->str_gen() << "->size()";
     }
     tt << ");" << endl;

   } 
   return tt.str();
// temp turn off check for testing task file
#if 0 
  } else if (index_.empty()) {
    throw logic_error ("TODO fix if index empty");
  } else if (loop.empty()) {
    throw logic_error ("TODO fix if loop empty");
  }  
#endif
}


string RDM::generate_mult(string indent, const string tag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, const string mlab) {
  stringstream tt;
  //indent += "  ";
  const string itag = "i";

  // now do the sort
  vector<string> close;

  // in case delta_ is not empty // //
  if (!delta_.empty()) {
    // first delta loops for blocks
    tt << indent << "if (";
    for (auto d = delta_.begin(); d != delta_.end(); ++d) {
      tt << d->first->str_gen() << " == " << d->second->str_gen() << (d != --delta_.end() ? " && " : "");
    }
    tt << ") {" << endl;
    close.push_back(indent + "}");
    indent += "  ";
   
    tt <<  make_get_block(indent);

    // start sort loops
    for (auto& i : index) {
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

    // extra for loops for merged 
    tt << make_merged_loops(indent, itag, index, merged, close);
    // make odata part of summation for target
    tt << make_odata(itag, indent, index);
    // mulitiply data and merge on the fly
    tt << multiply_merge(itag, indent, merged);
    // close loops
    for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
      tt << *iter << endl;

  // if delta_ is empty do sort_indices with merged multiplication // // 
  } else {
    
    tt << make_get_block(indent);

    // add loops over odata
    for (auto& i : index) {
      const int inum = i->num();
      tt << indent << "for (int " << itag << inum << " = 0; " << itag << inum << " != " << i->str_gen() << ".size(); ++" << itag << inum << ") {" << endl;
      close.push_back(indent + "}");
      indent += "  ";
    }
     
    // extra for loops for merged 
    tt << make_merged_loops(indent, itag, index, merged, close);
    // make odata part of summation for target
    tt << make_odata(itag, indent, index);
    // mulitiply data and merge on the fly
    tt << multiply_merge(itag, indent, merged);
    // close loops
    for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
      tt << *iter << endl;

  } 
  return tt.str();
}

// protected functions  //////
std::string RDM::make_get_block(std::string ident) {
  stringstream tt;
  tt << ident << "vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
  tt << ident << "unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;
  return tt.str();
}


std::string RDM::make_merged_loops(string& indent, const string itag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, vector<string>& close) {
  stringstream tt;
  bool mfound = false;
  for (auto& j : merged) {
    for (auto& i : index) {
      if (j->num() == i->num()) mfound = true; 
    }
    if (!mfound) {
      // add extra for loop(s) for merged
      const int jnum = j->num();
      tt << indent << "for (int " << itag << jnum << " = 0; " << itag << jnum << " != " << j->str_gen() << ".size(); ++" << itag << jnum << ") {" << endl;
      close.push_back(indent + "}");
      indent += "  ";
    } else {
      assert(false);
    }
  }
 return tt.str();
}


std::string RDM::multiply_merge(const string itag, string& indent,  const list<shared_ptr<Index> >& merged) {
  stringstream tt;
    // make data part of summation
    tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * data[";
    for (auto riter = index_.rbegin(); riter != index_.rend(); ++riter) {
      const string tmp = "+" + (*riter)->str_gen() + ".size()*(";
      tt << itag << (*riter)->num() << (riter != --index_.rend() ? tmp : "");
    }
    for (auto riter = ++index_.begin(); riter != index_.end(); ++riter)
      tt << ")";
    tt << "]";
    // multiply merge
    tt << " * " << "fdata" << "[";
    for (auto mi = merged.rbegin(); mi != merged.rend()  ; ++mi) { 
      const string tmp = "+" + (*mi)->str_gen() + ".size()*(";
      tt << itag << (*mi)->num() << (mi != --merged.rend() ? tmp : "");
    }
    for (auto mi = ++merged.begin(); mi != merged.end()  ; ++mi)  
      tt << ")";
    tt << "];" << endl;
  return tt.str();
}


std::string RDM::make_odata(const string itag, string& indent, const list<shared_ptr<Index> >& index) {
  stringstream tt;
  tt  << indent << "odata[";
  for (auto ri = index.rbegin(); ri != index.rend(); ++ri) {
    int inum = (*ri)->num();
    for (auto& d : delta_)
      if (d.first->num() == inum) inum = d.second->num();
    const string tmp = "+" + (*ri)->str_gen() + ".size()*(";
    tt << itag << inum << (ri != --index.rend() ? tmp : "");
  }
  for (auto ri = ++index.begin(); ri != index.end(); ++ri)
    tt << ")";
  tt << "]" << endl;
  return tt.str();
}


