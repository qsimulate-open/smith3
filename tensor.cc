//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tensor.cc
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


#include "tensor.h"
#include "constants.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace smith;

Tensor::Tensor(const shared_ptr<Op> op) : factor_(1.0), scalar_("")  {
  // scalar quantity..defined on bagel side
  // label
  label_ = op->label();
  // op
  for (auto& i : op->op()) {
    shared_ptr<Index> in = *get<0>(i);
    index_.push_back(in);
  }

}

static int ig = 0;

Tensor::Tensor(const shared_ptr<Active> activ) : factor_(1.0), scalar_("") {
  // scalar quantity..defined on bagel side
  // label
  stringstream ss; ss << "Gamma" << ig; ++ig;
  label_ = ss.str();
  // op
  index_ = activ->index();
  active_ = activ;
}


string Tensor::str() const {
  stringstream ss;
  if (factor_ != 1.0 || !scalar_.empty())
    ss << " " << fixed << setw(4) << setprecision(2) << factor_ << (scalar_.empty() ? "" : " "+scalar_) << " ";
  ss << label() << "(";
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    // we don't need the spin part here
    if (i != index_.begin()) ss << ", ";
    ss << (*i)->str(false);
  }
  ss << ")";

  if (merged_) ss << " << " << merged_->str();
  if (alias_) ss << " (" << label_ << ")";
  return ss.str();
}


// returns if all the indices are of active orbitals
bool Tensor::all_active() const {
  bool out = true;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    out &= (*i)->active();
  }
  return out;
}


// adds all-active tensor to Active_;
void Tensor::merge(shared_ptr<Tensor> a) {
  assert(active_);
  merged_ = a;
  list<list<shared_ptr<Index> >::iterator> remove;
  // remove Index that belongs to a
  for (auto& i : a->index()) {
    const int n = i->num();
    for (auto j = index_.begin(); j != index_.end(); ++j) {
      // careful, this code is driven by numbers
      if ((*j)->num() == n) {
        remove.push_back(j);
        break;
      }
    }
  }
  for (auto& i : remove) index_.erase(i);
}


bool Tensor::operator==(const Tensor& o) const {
  bool out = true;
  // if comparing Gammas we don't need to have similar labels or factors
  if (label_.find("Gamma") != string::npos && o.label().find("Gamma") != string::npos) {
  } else {
    out &= label() == o.label();
    out &= factor_ == o.factor();
  }
  // TODO for the time being, active is not assumed to be factorized
  // both objects should be active for tensors to be equal
  if (active() && o.active()) {
     // should return true or false
     out &= (*active() == *o.active());        
  }
  // index
  out &= index_.size() == o.index().size();
  if (out) {
    for (auto i = index_.begin(), j = o.index().begin(); i != index_.end(); ++i, ++j) {
      out &= (*i)->identical(*j);
    }
  }
  return out;
}


string Tensor::constructor_str(string indent) const {
  stringstream ss;
  ss << indent << "std::vector<IndexRange> " << label() << "_index";
  if (index_.empty()) {
    ss << ";" << endl;
  } else {
    ss << " = {";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << (i != index_.rbegin() ? ", this->" : "this->") << (*i)->generate();
    ss << "};" << endl;
  }
  ss << indent << "std::shared_ptr<Tensor<T> > " << label() << "(new Tensor<T>(" << label() << "_index, false));";
  return ss.str();
}

string Tensor::generate_get_block(const string cindent, const string lab, const bool move, const bool noscale) const {
  string lbl = label();
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  bool trans = false;
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
    trans = true;
  }

  stringstream tt;
  tt << cindent << "std::vector<size_t> " << lab << "hash = {";
  if (!trans) {
    tt << list_keys(index_);
  } else {
    assert(!(index_.size() & 1));
    for (auto iter = index_.rbegin(); iter != index_.rend(); ++iter) {
      if (iter != index_.rbegin()) tt << ", ";
      string i0 = (*iter)->str_gen() + ".key()";
      ++iter;
      tt << (*iter)->str_gen() << ".key()" << ", " << i0;
    }
  } 
  // for scalar.
  if (index_.empty() && merged_)
    tt << "0lu" ; 
  tt << "};" << endl;
  {
    tt << cindent << "std::unique_ptr<double[]> " << lab << "data = "
                  << lbl << "->" << (move ? "move" : "get") << "_block(" << lab << "hash);" << endl;
  }
  if (!scalar_.empty() && !noscale) {
    tt << cindent << "dscal_("; 
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      tt << (i != index_.rbegin() ? "*" : "") << (*i)->str_gen() << ".size()";
    // update scalar_ name directly
    tt << ", -" << scalar_ << "_, " << lab << "data.get(), 1);" << endl; 
  }
  return tt.str();
}


string Tensor::generate_scratch_area(const string cindent, const string lab, const bool zero) const {
  string lbl = label();
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
  }

  stringstream ss;
  ss << cindent << "std::unique_ptr<double[]> " << lab << "data_sorted(new double["
                << lbl << "->get_size(" << lab << "hash)]);" << endl;
  if (zero) {
    ss << cindent << "std::fill_n(" << lab << "data_sorted.get(), " << lbl << "->get_size(" << lab << "hash), 0.0);" << endl;
  }
  return ss.str();
}

string Tensor::generate_sort_indices(const string cindent, const string lab, const list<shared_ptr<Index> >& loop, const bool op) const {
  stringstream ss;
  if (!op) ss << generate_scratch_area(cindent, lab);

  vector<int> map(index_.size());
  ss << cindent << "sort_indices<";
  // determine mapping
  // first loop indices. order as in loop
  vector<int> done;
  for (auto i = loop.rbegin(); i != loop.rend(); ++i) {
    // count
    int cnt = 0;
    for (auto j = index_.rbegin(); j != index_.rend(); ++j, ++cnt) {
      if ((*i)->identical(*j)) break;
    }
    if (cnt == index_.size()) {
      //cout << "Potential problem Tensor::generate_sort_indices..did not find loop index in index_" << endl;
      throw logic_error("should not happen.. Tensor::generate_sort_indices");
    }
    done.push_back(cnt);
  }
  // then fill out others
  for (int i = 0; i != index_.size(); ++i) {
    if (find(done.begin(), done.end(), i) == done.end())
      done.push_back(i);
  }
  // if trans, transpose here!
  size_t found = label_.find("dagger");
  const bool trans = found != string::npos;
  if (trans && index_.size() & 1) throw logic_error("transposition not possible with 3-index objects");
  if (trans) {
    vector<int> tmp;
    for (int i = 0; i != done.size(); i += 2) {
      tmp.push_back(done[i+1]);
      tmp.push_back(done[i]);
    }
    done = tmp;
  }

  // then write them out.
  for (auto i = done.begin(); i != done.end(); ++i)
    ss << *i << ",";

  string target_label = op ? "odata" : lab + "data_sorted";

  ss << (op ? 1 : 0) << ",1," << prefac__(factor_);
  ss << ">(" << lab << "data, " << target_label;
  if (!trans) {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << ", " << (*i)->str_gen() << ".size()";
  } else {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
      string tmp = ", " + (*i)->str_gen() + ".size()";
      ++i;
      ss << ", " << (*i)->str_gen() << ".size()" << tmp;
    }
  }
  ss << ");" << endl;
  return ss.str();
}


string Tensor::generate_sort_indices_target(const string cindent, const string lab, const list<shared_ptr<Index> >& loop,
                                            const shared_ptr<Tensor> a, const shared_ptr<Tensor> b) const {
  stringstream ss;
  vector<int> map(index_.size());
  ss << cindent << "sort_indices<";

  // determine mapping
  // first obtain the ordering of indices from dgemm
  list<shared_ptr<Index> > source;
  {
    list<shared_ptr<Index> > aind = a->index();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto& j : loop)
        if ((*i)->identical(j)) found = true;
      if (!found) source.push_back(*i);
    }
    aind = b->index();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto& j : loop)
        if ((*i)->identical(j)) found = true;
      if (!found) source.push_back(*i);
    }
  }

  for (auto j = index_.rbegin(); j != index_.rend(); ++j) {
    // count
    int cnt = 0;
    // note that source is already reversed!
    for (auto i = source.begin(); i != source.end(); ++i, ++cnt) {
      if ((*i)->identical(*j)) break;
    }
    if (cnt == index_.size()) throw logic_error("should not happen.. Tensor::generate_sort_indices_target");
    ss << cnt << ",";
  }

  ss << "1,1," << prefac__(factor_);
  ss << ">(" << lab << "data_sorted, " << lab << "data";
  for (auto i = source.begin(); i != source.end(); ++i) ss << ", " << (*i)->str_gen() << ".size()";
  ss << ");" << endl;
  return ss.str();
}


pair<string, string> Tensor::generate_dim(const list<shared_ptr<Index> >& di) const {
  vector<string> s, t;
  // first indices which are not shared
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    bool shared = false;
    for (auto& j : di) {
      if ((*i)->identical(j)) {
        shared = true;
        break;
      }
    }
    if (shared) {
      t.push_back((*i)->str_gen() + ".size()");
    } else {
      s.push_back((*i)->str_gen() + ".size()");
    }
  }

  stringstream ss, tt;
  for (auto i = s.begin(); i != s.end(); ++i) {
    if (i != s.begin()) ss << "*";
    ss << *i;
  }
  for (auto i = t.begin(); i != t.end(); ++i) {
    if (i != t.begin()) tt << "*";
    tt << *i;
  }
  return make_pair(ss.str(), tt.str());
}


string Tensor::generate_active(string indent, const string tag, const bool use_blas) const {
  assert(label_.find("Gamma") != string::npos);
  stringstream tt;
  if (!merged_) {
    tt << active()->generate(indent, tag, index());
  } else {
    
    tt << indent <<"// associated with merged" << endl;
    // add merge loops
    vector<string> mclose;
    list<shared_ptr<Index> >& merged = merged_->index();
    // the m->generate() makes active_ this label must be identical on bagel side (src/smith/spinfreebase.h)
    for (auto& m : merged) {
      tt << indent << "for (auto& " << m->str_gen() << " : "  <<  m->generate()  << ") {" << endl;
      mclose.push_back(indent + "}");
      indent += "  ";
    }
    
    // add fdata 
    tt << indent << "std::vector<size_t> fhash = {" << list_keys(merged) << "};" << endl;
    tt << indent << "std::unique_ptr<double[]> fdata = " << merged_->label() << "->get_block(fhash);" << endl;
   
    if (use_blas) {
      tt << indent << "std::unique_ptr<double[]> fdata_sorted(new double["<< merged_->label() << "->get_size(fhash)]);" << endl;

      // make sort_indices for merged op
      vector<int> done;
      tt << indent << "sort_indices<";
      for (auto i = merged.begin(); i != merged.end(); ++i) {
        int cnt = 0;
        for (auto j = merged.rbegin(); j != merged.rend(); ++j, ++cnt) {
          if ((*i)->identical(*j)) break;
        }
        if (cnt == merged.size()) throw logic_error("should not happen... fdata sort indices");
        done.push_back(cnt);
      }
      // then fill out others
      for (int i = 0; i != merged.size(); ++i) {
        if (find(done.begin(), done.end(), i) == done.end())
          done.push_back(i);
      }
      // write out
      for (auto& i : done) 
        tt << i << ",";

      tt << "0,1,1,1";
   
      // add source data dimensions
      tt << ">(fdata, fdata_sorted, " ;
      for (auto iter = merged.rbegin(); iter != merged.rend(); ++iter) {
        if (iter != merged.rbegin()) tt << ", ";
          tt << (*iter)->str_gen() << ".size()";
      }
      tt << ");" << endl;
    } 
    
    // generate  merged and/or rdm 
    tt << active()->generate(indent, tag, index(), merged_->index(), merged_->label(), use_blas);

    // close merge for loops
    for (auto iter = mclose.rbegin(); iter != mclose.rend(); ++iter)
      tt << *iter << endl;

  }
  return tt.str();
}


string Tensor::generate_loop(string& indent, vector<string>& close) const {
  stringstream tt;
  for (auto iter = index_.begin(); iter != index_.end(); ++iter, indent += "  ") {
    string cindex = (*iter)->str_gen();
    tt << indent << "for (auto& " << cindex << " : " << (*iter)->generate() << ") {" << endl;
    close.push_back(indent + "}");
  }
  return tt.str();
}


string Tensor::generate_gamma(const int ic, const bool enlist, const bool use_blas) const {
  assert(label_.find("Gamma") != string::npos);

  stringstream tt;
  vector<int> rdmn = active()->required_rdm();

  tt << "template <typename T>" << endl;
  if (!enlist) {
    tt << "class Task" << ic << " : public Task<T> {" <<  "  // associated with gamma" << endl;
  } else {
    tt << "class Task" << ic << " : public EnergyTask<T> {" <<  "  // associated with gamma" << endl;
  }
  tt << "  protected:" << endl;
  tt << "    IndexRange closed_;" << endl;
  tt << "    IndexRange active_;" << endl;
  tt << "    IndexRange virt_;" << endl;
  tt << "    std::shared_ptr<Tensor<T> > " << label() << ";" << endl;
  for (auto& i: rdmn)
      tt << "    std::shared_ptr<Tensor<T> > rdm" << i << ";" << endl;
  if (merged_) 
    tt << "    std::shared_ptr<Tensor<T> > " << merged_->label() << ";" << endl;
  tt << endl;
  // loops
  tt << "    void compute_() {" << endl;

  string indent ="      ";
  vector<string> close;
  tt << generate_loop(indent, close);
  // generate gamma get block, true does a move_block
  tt << generate_get_block(indent, "o", true, true); // first true means move, second true means we don;t scale
  if (merged_) {
    if (use_blas && !index_.empty()) tt << generate_scratch_area(indent,"o",true);
  }
  // now generate codes for rdm
  tt << generate_active(indent, "o", use_blas);

  // generate gamma put block
  tt << indent << label() << "->put_block(ohash, odata);" << endl;
  // close the loops
  for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
    tt << *iter << endl;
  tt << "    };  " << endl;
  tt << "" << endl;
  // done with protected part
  tt << "" << endl;
  tt << "  public:" << endl;
  if (!enlist) {
    tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {" << endl;
  } else {
    tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {" << endl;
  }
  tt << "      closed_ = i[0];" << endl;
  tt << "      active_ = i[1];" << endl;
  tt << "      virt_   = i[2];" << endl;
  tt << "      " << label() << "  = t[0];" << endl;
  //this is related to what rdms we have
  {
    int itmp = 1;
    for (auto i = rdmn.begin(); i != rdmn.end(); ++i, ++itmp) {
      tt << "      rdm" << (*i) << "    = t[" << itmp << "];" << endl;
    }
    // related to merged tensor
    if (merged_) tt << "      "<< merged_->label() << "      = t[" <<  itmp << "];" << endl;
  }
  tt << "    };" << endl;
  tt << "    ~Task" << ic << "() {};" << endl;
  tt << "};" << endl << endl;
  return tt.str();
}


