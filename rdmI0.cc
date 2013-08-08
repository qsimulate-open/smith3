//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: rdmI0.cc
// Copyright (C) 2013 Matthew MacLeod
//
// Author: Matthew MacLeod <matthew.macleod@northwestern.edu>
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
#include <map>

using namespace std;
using namespace smith;

shared_ptr<RDM> RDMI0::copy() const {
  // first clone all the indices
  list<shared_ptr<Index>> in;
  for (auto& i : index_) in.push_back(i->clone());

  map<shared_ptr<Spin>, shared_ptr<Spin>> dict;
  auto j = in.begin();
  for (auto i = index_.begin(); i != index_.end(); ++i, ++j) {
    // get the original spin
    shared_ptr<Spin> o = (*i)->spin();
    if (dict.find(o) == dict.end()) {
      (*j)->set_spin(o);
    } else {
      shared_ptr<Spin> s(new Spin());
      s->set_num(o->num());
      dict.insert(make_pair(o,s));
      (*j)->set_spin(s);
    }
  }

  // lastly clone all the delta functions
  map<shared_ptr<const Index>, shared_ptr<const Index>> d;
  for (auto& i : delta_) d.insert(make_pair(i.first->clone(), i.second->clone()));

  list<shared_ptr<const Index>> inc;
  for (auto& i : in) inc.push_back(i);

  shared_ptr<RDM> out(new RDMI0(inc, d, make_pair(bra_, ket_)));
  out->fac() = fac_;
  return out;
}

//
// An application of "Wick's theorem"
// This function is controlled by const Index::num_. Not a great code, it could have been driven by pointers... lazy me.
//
list<shared_ptr<RDM>> RDMI0::reduce_one(list<int>& done) const {
  // first find non-daggered operator which is not aligned
  list<shared_ptr<RDM>> out;

  for (auto i = index_.begin(); i != index_.end(); ++i) {

    // looking for a non-daggered index that is not registered in done
    if ((*i)->dagger() || find(done.begin(), done.end(), (*i)->num()) != done.end())
      continue;

    // again this function is controlled by numbers... sorry...
    const int inum = (*i)->num();

    for (auto j = i; j != index_.end(); ++j) {
      if (!(*j)->dagger() || j==i) continue;

      // if you find daggered object in the right hand side...
      shared_ptr<RDM> tmp = this->copy();

      // find the indices to be deleted.
      vector<list<shared_ptr<const Index>>::iterator> rml;
      int cnt0 = -1;
      int cnt = 0;
      for (auto k = tmp->index().begin(); k != tmp->index().end(); ++k, ++cnt) {
        if ((*k)->same_num(*i) || (*k)->same_num(*j)) {
          cnt0 = cnt0 >= 0 ? cnt-cnt0 : cnt;
          rml.push_back(k);
        }
      }
      assert(rml.size() == 2);

      tmp->delta().insert(make_pair(*rml[0],*rml[1]));
      // Please note that this procedure does not change the sign (you can prove it in 30sec)
      tmp->fac() *= ((cnt0-1)&1 ? -1.0 : 1.0);
      if ((*i)->same_spin(*j)) {
        tmp->fac() *= 2.0;
      } else {
        // this case we need to replace a spin
        const shared_ptr<Spin> s0 = (*rml[0])->spin();
        const shared_ptr<Spin> s1 = (*rml[1])->spin();
        for (auto& k : tmp->index()) {
          if (k->spin() == s0)
            k->set_spin(s1);
        }
      }

      // erasing indices which are push-backed in delta
      tmp->index().erase(rml[0]);
      tmp->index().erase(rml[1]);
      out.push_back(tmp);
    }
    done.push_back((*i)->num());
    break;
  }
  return out;
}


string RDMI0::generate(string indent, const string tag, const list<shared_ptr<const Index>>& index, const list<shared_ptr<const Index>>& merged, const string mlab, vector<string> in_tensors, const bool use_blas) {
  return merged.empty() ? generate_not_merged(indent, tag, index, in_tensors) : generate_merged(indent, tag, index, merged, mlab, in_tensors, use_blas);
}


string RDMI0::generate_not_merged(string indent, const string tag, const list<shared_ptr<const Index>>& index, vector<string> in_tensors) {
  stringstream tt;
  tt << indent << "{" << endl;
  const string lindent = indent;

  indent += "  ";
  const string itag = "i";

  list<shared_ptr<const Index>> dindex = index;
  // todo replace all index_

  // now do the sort
  vector<string> close;

  // in case delta_ is not empty
  if (!delta_.empty()) {
#define debug_rdm_tasks
#ifdef debug_rdm_tasks
    if (rank() == 0) tt << indent <<  "// rdm0 non-merged case" << endl;
#endif

    // first delta if statement
    tt << make_delta_if(indent, close);

    if (!index_.empty() && rank() != 0) {
      stringstream zz;
      zz << "rdm" << rank();
      string rlab = zz.str();

      // make map of in_tensors
      map<string, string> inlab;
      map_in_tensors(in_tensors, inlab);

      tt << make_get_block(indent, "i0", inlab[rlab], index);
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

    // make map of in_tensors
    map<string, string> inlab;
    map_in_tensors(in_tensors, inlab);

    if (rank() != 0) {
      stringstream zz;
      zz << "rdm" << rank();
      string rlab = zz.str();
      tt << make_get_block(indent, "i0", inlab[rlab], index);
    }

    // do sort_indices here
    vector<int> done;
    tt << indent << "sort_indices<";
    for (auto i = index.rbegin(); i != index.rend(); ++i) {
      int cnt = 0;
      for (auto j = dindex.rbegin(); j != dindex.rend(); ++j, ++cnt) {
        if ((*i)->identical(*j)) break;
      }
      if (cnt == index.size()) throw logic_error("should not happen in sort RDMI0::generate");
      done.push_back(cnt);
    }
    // then fill out others
    for (int i = 0; i != index.size(); ++i) {
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
    for (auto iter = index.rbegin(); iter != index.rend(); ++iter) {
      if (iter != index.rbegin()) tt << ", ";
        tt << (*iter)->str_gen() << ".size()";
    }
    tt << ");" << endl;
  }

  tt << lindent << "}" << endl;

  return tt.str();
}


string RDMI0::generate_merged(string indent, const string tag, const list<shared_ptr<const Index>>& index, const list<shared_ptr<const Index>>& merged, const string mlab, vector<string> in_tensors, const bool use_blas) {
  stringstream tt;
  //indent += "  ";
  const string itag = "i";
  const string lindent = indent;
  // now do the sort
  vector<string> close;

#if 1
  if (rank() == 0) tt << indent << "// rdm0 merged case" << endl;
#endif

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

  // make map of in_tensors
  map<string, string> inlab;
  map_in_tensors(in_tensors, inlab);

  if (!use_blas) {
    if (rank() != 0) {
      tt <<  make_get_block(indent, "i0", inlab[rlab], index);
    }
    // loops for index and merged
    tt << make_merged_loops(indent, itag, close);
    // make odata part of summation for target
    tt << make_odata(itag, indent, index);
    // mulitiply data and merge on the fly
    tt << multiply_merge(itag, indent, merged);
  } else {
    if (rank() != 0) {
      tt << make_get_block(indent, "i0", inlab[rlab], index);
      tt << make_sort_indices(indent, "i0", merged);
      tt << endl;

      tt << make_blas_multiply(indent, merged, index);
      tt << endl;

      if (!index.empty()) {
        // determine mapping
        list<shared_ptr<const Index>> source;
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
            shared_ptr<const Index> first_partner;
            shared_ptr<const Index> second_partner;
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
            if (cnt == index.size()) throw logic_error("should not happen.. RDMI0 odata target, delta case");
            done.push_back(cnt);
          } else {
            int cnt = 0;
            for (auto i = index.rbegin(); i != index.rend(); ++i, ++cnt) {
              if ((*i)->identical(*j)) break;
            }
            if (cnt == index.size()) throw logic_error("should not happen.. RDMI0 odata target, non-delta case");
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
      list<shared_ptr<const Index>> source;
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
        shared_ptr<const Index> first_partner;
        shared_ptr<const Index> second_partner;
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
          throw logic_error("should not happen for merged.. RDMI0::generate");
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


////// previously protected functions start //////
void RDMI0::map_in_tensors(std::vector<string> in_tensors, map<string,string>& inlab) {
  // make map of in_tensors
  int icnt = 0;
  for (auto& i : in_tensors) {
    stringstream zz; zz << "in(" << icnt << ")";
    inlab[i] = zz.str();
    ++icnt;
  }
}


string RDMI0::make_get_block(string indent, string tag, string lbl, const list<shared_ptr<const Index>>& index) {
  stringstream tt;
  tt << indent << "std::unique_ptr<double[]> " << tag << "data = " << lbl << "->get_block(";
  for (auto i = index.rbegin(); i != index.rend(); ++i) {
    if (i != index.rbegin()) tt << ", ";
    tt << (*i)->str_gen();
  }
  tt << ");" << endl;
  return tt.str();
}


string RDMI0::make_blas_multiply(string dindent, const list<shared_ptr<const Index>>& loop, const list<shared_ptr<const Index>>& index) {
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
    // TODO need to check this when general case is available
    tt << dindent << "odata[0] = " << setprecision(1) << fixed << factor() <<  " * " << "ddot_(" << t1.first << ", fdata_sorted, 1, i0data_sorted, 1);" << endl;
#endif
  }
  return tt.str();
}

pair<string, string> RDMI0::get_dim(const list<shared_ptr<const Index>>& di, const list<shared_ptr<const Index>>& index) const {
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


string RDMI0::make_sort_indices(string indent, string tag, const list<shared_ptr<const Index>>& loop) {
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


string RDMI0::make_merged_loops(string& indent, const string itag, vector<string>& close) {
  stringstream tt;

  // gather all the loop indices
  list<shared_ptr<const Index>> loop;
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


string RDMI0::multiply_merge(const string itag, string& indent, const list<shared_ptr<const Index>>& merged) {
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


string RDMI0::fdata_mult(const string itag, const list<shared_ptr<const Index>>& merged) {
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


string RDMI0::make_odata(const string itag, string& indent, const list<shared_ptr<const Index>>& index) {
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


string RDMI0::make_sort_loops(const string itag, string& indent, const list<shared_ptr<const Index>>& loop, vector<string>&  close) {
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


string RDMI0::make_delta_if(string& indent, vector<string>& close) {
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
