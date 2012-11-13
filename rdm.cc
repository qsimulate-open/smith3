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


string RDM::generate(string indent, const string tag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, const string mlab, const bool use_blas) {
  return merged.empty() ? generate_not_merged(indent, tag, index) : generate_merged(indent, tag, index, merged, mlab, use_blas);
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

    if (!index_.empty() && rank() != 0) {
      stringstream zz;
      zz << "rdm" << rank();
      string rlab = zz.str();
      tt << make_get_block(indent, "i0", rlab);
    }    

    // loops over delta indices
    tt << make_sort_loops(itag, indent, index, close); 

    // make odata part of summation for target
    tt << make_odata(itag, indent, index);

    // make data part of summation
    if (index_.empty()) {
      tt << "  += " << setprecision(1) << fixed << factor() << ";" << endl;
    } else {
      tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * i0data[";
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

    if (rank() != 0) {
      stringstream zz;
      zz << "rdm" << rank();
      string rlab = zz.str();
      tt << make_get_block(indent, "i0", rlab);
    }
 
    // do sort_indices here
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
    tt << ">(i0data, " << tag << "data, " ;
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
        tt << (*iter)->str_gen() << ".size()";
    }
    tt << ");" << endl;
  } 

  tt << lindent << "}" << endl;

  return tt.str();
}


string RDM::generate_merged(string indent, const string tag, const list<shared_ptr<Index> >& index, const list<shared_ptr<Index> >& merged, const string mlab, const bool use_blas) {
  stringstream tt;
  //indent += "  ";
  const string itag = "i";
  const string lindent = indent;
  // now do the sort
  vector<string> close;


  if (rank() == 0)
    tt << indent << "// rdm0 merged case" << endl;

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
  stringstream zz;
  zz << "rdm" << rank();
  string rlab = zz.str();

  if (!use_blas) {
    if (rank() !=0)
      tt <<  make_get_block(indent, "i0", rlab);
    // loops for index and merged 
    tt << make_merged_loops(indent, itag, close);
    // make odata part of summation for target
    tt << make_odata(itag, indent, index);
    // mulitiply data and merge on the fly
    tt << multiply_merge(itag, indent, merged);
  } else {
    if (rank() != 0) {
      tt << make_get_block(indent, "i0", rlab);
      tt << make_scratch_area(indent,"i0", rlab);
      tt << make_sort_indices(indent, "i0", merged);
      tt << endl;

      tt << make_blas_multiply(indent, merged, index);
      tt << endl;

      if (!index.empty()) {
        // determine mapping
        list<shared_ptr<Index> > source;
        vector<int> done;

        // compare rdm and merged indices
        for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
          bool found = false;
          for (auto& j : merged)
            if ((*i)->identical(j)) found = true;
          if (!found) source.push_back(*i);
        } 

        // complete source
        for (auto i = index.rbegin(); i != index.rend(); ++i) {
          bool found = false;
          for (auto& j : source)
            if ((*i)->identical(j)) found = true;
          if (!found) source.push_back(*i);
        }

        // go through odata target indices
        for (auto j = source.begin(); j != source.end(); ++j) {
          // check delta mapping
          if (!delta_.empty()) {
            bool matched_first = false;
            bool matched_second = false;
            shared_ptr<Index> first_partner;
            shared_ptr<Index> second_partner;
            for (auto d = delta_.begin(); d != delta_.end(); ++d) {
              if ((d->first)->identical(*j)) {
                matched_first = true;
                first_partner = d->second;
              }
              if ((d->second)->identical(*j)) {
                matched_second = true;
                second_partner = d->first;
              }
            }
            int cnt = 0;
            if (!matched_first && !matched_second) {
              for (auto i = index.rbegin(); i != index.rend(); ++i, ++cnt) 
                if ((*i)->identical(*j)) break;
            } else if (matched_first) {
              for (auto i = index.rbegin(); i != index.rend(); ++i, ++cnt) 
                if ((*i)->identical(*j) || first_partner->identical(*j)) break;
            } else if (matched_second) {
              for (auto i = index.rbegin(); i != index.rend(); ++i, ++cnt) 
                if ((*i)->identical(*j) || second_partner->identical(*j)) break;
            }  
            if (cnt == index.size()) throw logic_error("should not happen.. RDM odata target, delta case");
            done.push_back(cnt);
          } else {
            int cnt = 0;
            for (auto i = index.rbegin(); i != index.rend(); ++i, ++cnt) {
              if ((*i)->identical(*j)) break;
            }
            if (cnt == index.size()) throw logic_error("should not happen.. RDM odata target, non-delta case");
            done.push_back(cnt);
          }
        }
        tt << indent << "sort_indices<";
        for (auto& i : done) { 
          tt << i << ",";
        }
        tt << "1,1,1,1>(odata_sorted, odata";
        for (auto i = source.begin(); i != source.end(); ++i) tt << ", " << (*i)->str_gen() << ".size()";
        tt << ");" << endl;
        tt << endl;
      }

    } else {
      // for rdm0 case 
      // add dscal
      tt << indent << "dscal_("; 
      for (auto i = merged.rbegin(); i != merged.rend(); ++i)
        tt << (i != merged.rbegin() ? "*" : "") << (*i)->str_gen() << ".size()";
      tt << ", " << setprecision(1) << fixed << factor() << ", " << "fdata_sorted.get(), 1);" << endl; 
  
      // sort back to target layout (odata layout)
      list<shared_ptr<Index> > source;
      vector<int> done;
      for (auto i = index.rbegin(); i != index.rend(); ++i) {
        bool found = false;
        for (auto& j : merged)
          if ((*i)->identical(j)) found = true;
        if (!found) source.push_back(*i);
      } 
      // complete source
      for (auto i = index.rbegin(); i != index.rend(); ++i) {
        bool found = false;
        for (auto& j : source)
          if ((*i)->identical(j)) found = true;
        if (!found) source.push_back(*i);
      }
  
  
      for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
        // also check if in deltas
        bool matched_first = false;
        bool matched_second = false;
        shared_ptr<Index> first_partner;
        shared_ptr<Index> second_partner;
        for (auto d = delta_.begin(); d != delta_.end(); ++d) {
          if ((d->first)->identical(*i)) {
            matched_first = true;
            first_partner = d->second;
          }
          if ((d->second)->identical(*i)) {
            matched_second = true;
            second_partner = d->first;
          }
        }
        int cnt = 0;
        if (!matched_first && !matched_second) {
          for (auto j = index.rbegin(); j != index.rend(); ++j, ++cnt) 
            if ((*i)->identical(*j)) break;
        } else if (matched_first) {
          for (auto j = index.rbegin(); j != index.rend(); ++j, ++cnt) 
            if ((*i)->identical(*j) || first_partner->identical(*j)) break;
        } else if (matched_second) {
          for (auto j = index.rbegin(); j != index.rend(); ++j, ++cnt) 
            if ((*i)->identical(*j) || second_partner->identical(*j)) break;
        }  
        if (cnt == index.size()) 
          throw logic_error("should not happen.. RDM::generate");
        done.push_back(cnt);
      }
      // then fill out others
      for (int i = 0; i != index.size(); ++i) {
        if (find(done.begin(), done.end(), i) == done.end())
          done.push_back(i);
      }
      // write out
      tt << indent << "sort_indices<";
      for (auto& i : done) 
        tt << i << ",";
      // add factor information
      tt << "1,1,1,1";
      // add source data dimensions
      tt << ">(fdata_sorted, odata, " ;
      for (auto iter = source.begin(); iter != source.end(); ++iter) {
        if (iter != source.begin()) tt << ", ";
          tt << (*iter)->str_gen() << ".size()";
      }
      tt << ");" << endl;
    }
  } 
  // close loops
  for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
    tt << *iter << endl;

  if (delta_.empty()) tt << lindent << "}" << endl;

  return tt.str();
}


// protected functions start //////
string RDM::make_get_block(string indent, string tag, string lbl) {
  stringstream tt;
  tt << indent << "std::vector<size_t> "<< tag << "hash = {" << list_keys(index_) << "};" << endl;
  tt << indent << "std::unique_ptr<double[]> " << tag << "data = " << lbl << "->get_block(i0hash);" << endl;
  return tt.str();
}

string RDM::make_scratch_area(string indent, string tag, string lbl) {
  stringstream tt;
  tt << indent << "std::unique_ptr<double[]> " << tag << "data_sorted(new double["
                << lbl << "->get_size(" << tag << "hash)]);" << endl;
  return tt.str();
}

string RDM::make_blas_multiply(string dindent, const list<shared_ptr<Index> >& loop, const list<shared_ptr<Index> >& index) {
  stringstream tt;
    
  pair<string,string> t1 = get_dim(loop, index);
  // call dgemv  
  if (t1.second != "") {
    tt << dindent << "dgemv_(\"T\", ";
    string tt1 = t1.first == "" ? "1" : t1.first;
    string tt2 = t1.second== "" ? "1" : t1.second;
    tt << tt1 << ", " << tt2 << ", " << endl;
    tt << dindent << "       " << setprecision(1) << fixed << factor() << ", i0data_sorted, " << tt1 << ", fdata_sorted, 1.0, "  << endl
       << dindent << "       1.0, odata_sorted, 1.0);" << endl;
  } else {
    // add check.. 
#if 0
    throw logic_error("ddot needed in f1 merged");
#else  
    // todo need to check this when general case is available
    tt << dindent << "odata[0] = " << setprecision(1) << fixed << factor() <<  " * " << "ddot_(" << t1.first << ", fdata_sorted, 1, i0data_sorted, 1);" << endl;
#endif
  }
  return tt.str();
}

pair<string, string> RDM::get_dim(const list<shared_ptr<Index> >& di, const list<shared_ptr<Index> >& index) const {
  vector<string> s, t;

  // get shared (f1) indices, equal to number of rows in matrix
  for (auto& i : di) {
     s.push_back((i)->str_gen() + ".size()");
  }
  stringstream ss;
  for (auto i = s.begin(); i != s.end(); ++i) {
    if (i != s.begin()) ss << "*";
    ss << *i;
  }

  // get number of columns in matrix
  for (auto i = index.rbegin(); i != index.rend(); ++i) {
    bool shared = false;
    for (auto& j : di) {
      if ((*i)->identical(j)) {
        shared = true;
        break;
      }
    }
    if (!shared) 
      t.push_back((*i)->str_gen() + ".size()");
  }
  stringstream tt;
  for (auto i = t.begin(); i != t.end(); ++i) {
    if (i != t.begin()) tt << "*";
    tt << *i;
  }
  return make_pair(ss.str(), tt.str());
}


string RDM::make_sort_indices(string indent, string tag, const list<shared_ptr<Index> >& loop) {
  stringstream tt;
    vector<int> done;
    // then fill out others
    for (int i = 0; i != index_.size(); ++i) {
      if (find(done.begin(), done.end(), i) == done.end())
        done.push_back(i);
    }
    // write out
    tt << indent << "sort_indices<";
    for (auto& i : done) 
      tt << i << ",";

    // add factor information
    tt << "0,1,1,1";
 
    // add source data dimensions
    tt << ">(i0data, " << tag << "data_sorted, " ;
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
        tt << (*iter)->str_gen() << ".size()";
    }
    tt << ");" << endl;
  return tt.str();
}


string RDM::make_merged_loops(string& indent, const string itag, vector<string>& close) {
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
  if (rank() == 0) {
    tt << "  += " << setprecision(1) << fixed << factor();
    tt << fdata_mult(itag, merged);
  } else { 
    // make data part of summation
    tt << indent << "  += (" << setprecision(1) << fixed << factor() << ") * i0data[";
    for (auto riter = index_.rbegin(); riter != index_.rend(); ++riter) {
      const string tmp = "+" + (*riter)->str_gen() + ".size()*(";
      tt << itag << (*riter)->num() << (riter != --index_.rend() ? tmp : "");
    }
    for (auto riter = ++index_.begin(); riter != index_.end(); ++riter)
      tt << ")";
    tt << "]";
    // multiply merge
    tt << fdata_mult(itag, merged);
  }
  return tt.str();
}

string RDM::fdata_mult(const string itag, const list<shared_ptr<Index> >& merged) {
  stringstream tt;
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
  // factor will now be added on same line in rdm0 case
  tt << "]";
  if (rank() > 0) tt << endl; 
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
