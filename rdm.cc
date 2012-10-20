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


bool RDM::operator==(const RDM& o) const {
  bool out = true;
  // compare all rdms of active objects
  out &= fac_ == o.factor();
  out &= index_.size() == o.index().size();
  out &= delta_.size() == o.delta().size();
  if (index_.size() == o.index().size()) { 
    for (auto i = index_.begin(), j = o.index().begin(); i != index_.end(); ++i, ++j)
      out &= (*i)->identical(*j);
  } else {
    out &= false;
  }
  return out;
}


string RDM::generate(string indent, const string tag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, const string mlab) {
  return merged.empty() ? generate_not_merged(indent, tag, index) : generate_merged(indent, tag, index, merged, mlab);
}


string RDM::generate_not_merged(string indent, const string tag, const list<shared_ptr<Index> >& index) {
  stringstream tt;
  tt << indent << "{" << endl;
  const string lindent = indent;

  indent += "  ";
  const string itag = "i";

  // now do the sort
  vector<string> close;

  // in case delta_ is not empty
  if (!delta_.empty()) {

    // first delta if statement 
    tt << make_delta_if(indent, close);

    if (!index_.empty())
      tt << make_get_block(indent);
   
    // loops over delta indices
    tt << make_sort_loops(itag, indent, index, close); 

    // make odata part of summation for target
    tt << make_odata(itag, indent, index);

    // make data part of summation
    if (index_.empty()) {
      tt << indent << "  += " << setprecision(1) << fixed << factor() << ";" << endl;
    } else {
      tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * data[";
      for (auto riter = index_.rbegin(); riter != index_.rend(); ++riter) {
        const string tmp = "+" + (*riter)->str_gen() + ".size()*(";
        tt << itag << (*riter)->num() << (riter != --index_.rend() ? tmp : "");
      }
      for (auto riter = ++index_.begin(); riter != index_.end(); ++riter)
        tt << ")";
      tt << "];" << endl;
    }
 
    // close loops
    for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
      tt << *iter << endl;

  // if delta_ is empty call sort_indices
  } else {
    // loop up the operator generators

    tt <<  make_get_block(indent);
 
    // call sort_indices here
    vector<int> done;
    tt << indent << "sort_indices<";
    for (auto i = index.rbegin(); i != index.rend(); ++i) {
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
    tt << ">(data, " << tag << "data, " ;
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
        tt << (*iter)->str_gen() << ".size()";
    }
    tt << ");" << endl;
  } 

  tt << lindent << "}" << endl;

  return tt.str();
}


string RDM::generate_merged(string indent, const string tag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, const string mlab) {
  stringstream tt;
  //indent += "  ";
  const string itag = "i";
  const string lindent = indent;

  // now do the sort
  vector<string> close;

  // first delta loops for blocks
  if (!delta_.empty()) {
    tt << indent << "if (";
    for (auto d = delta_.begin(); d != delta_.end(); ++d) {
      tt << d->first->str_gen() << " == " << d->second->str_gen() << (d != --delta_.end() ? " && " : "");
    }
    tt << ") {" << endl;
    close.push_back(indent + "}");
  } else {
    tt << indent << "{" << endl;
  }
 
  indent += "  ";
  tt <<  make_get_block(indent);

  // loops for index and merged 
  tt << make_merged_loops(indent, itag, index, merged, close);
  // make odata part of summation for target
  tt << make_odata(itag, indent, index);
  // mulitiply data and merge on the fly
  tt << multiply_merge(itag, indent, merged);
  // close loops
  for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
    tt << *iter << endl;

  if (delta_.empty()) tt << lindent << "}" << endl;

  return tt.str();
}


// protected functions start //////
string RDM::make_get_block(string indent) {
  stringstream tt;
  tt << indent << "std::vector<size_t> i0hash = {" << list_keys(index_) << "};" << endl;
  tt << indent << "std::unique_ptr<double[]> data = rdm" << rank() << "->get_block(i0hash);" << endl;
  return tt.str();
}


string RDM::make_merged_loops(string& indent, const string itag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, vector<string>& close) {
  stringstream tt;

  // gather all the loop indices
  list<shared_ptr<Index> > loop;
  for (auto& i : index_) {
    bool found = false;
    for (auto& j : delta_) {
      // second index in deltas will not be looped
      if (j.first->num() == i->num() || j.second->num() == i->num()) {
        found = true;
        break;
      }
    }
    if (!found) loop.push_back(i);
  }
  for (auto& j : delta_)
    loop.push_back(j.second);

  // generate loops
  for (auto& i : loop) {
    const int inum = i->num();
    tt << indent << "for (int " << itag << inum << " = 0; " << itag << inum << " != " << i->str_gen() << ".size(); ++" << itag << inum << ") {" << endl;
    close.push_back(indent + "}");
    indent += "  ";
  }

  return tt.str();
}


string RDM::multiply_merge(const string itag, string& indent, const list<shared_ptr<Index> >& merged) {
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
    int inum = (*mi)->num();
    for (auto& i : delta_)
      if (i.first->num() == inum) inum = i.second->num(); 
    const string tmp = "+" + (*mi)->str_gen() + ".size()*(";
    tt << itag << inum << (mi != --merged.rend() ? tmp : "");
  }
  for (auto mi = ++merged.begin(); mi != merged.end()  ; ++mi)  
    tt << ")";
  tt << "];" << endl;
  return tt.str();
}


string RDM::make_odata(const string itag, string& indent, const list<shared_ptr<Index> >& index) {
  stringstream tt;

  tt  << indent << "odata[";
  if (index.empty()) {
    tt << "0" ;
  } else { 
    for (auto ri = index.rbegin(); ri != index.rend(); ++ri) {
      int inum = (*ri)->num();
      for (auto& d : delta_)
        if (d.first->num() == inum) inum = d.second->num();
      const string tmp = "+" + (*ri)->str_gen() + ".size()*(";
      tt << itag << inum << (ri != --index.rend() ? tmp : "");
    }
  }
  for (auto ri = ++index.begin(); ri != index.end(); ++ri)
    tt << ")";
    tt << "]" << endl;
  return tt.str();
}

string RDM::make_sort_loops(const string itag, string& indent, const list<shared_ptr<Index> >& loop, vector<string>&  close) {
  stringstream tt;
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
  return tt.str();
}

string RDM::make_delta_if(string& indent, vector<string>& close) {
  stringstream tt;

  tt << indent << "if (";
  for (auto d = delta_.begin(); d != delta_.end(); ++d) {
    tt << d->first->str_gen() << " == " << d->second->str_gen() << (d != --delta_.end() ? " && " : "");
  }
  tt << ") {" << endl;
  close.push_back(indent + "}");
  indent += "  ";

  return tt.str();
}
