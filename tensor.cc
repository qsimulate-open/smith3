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

#define debug_tasks

using namespace std;
using namespace smith;


Tensor::Tensor(const shared_ptr<Operator> op) : factor_(1.0), scalar_("")  {
  // scalar quantity..defined on bagel side
  // label
  label_ = op->label();
  // op
  for (auto& i : op->op()) {
    shared_ptr<const Index> in = *get<0>(i);
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
  list<list<shared_ptr<const Index>>::iterator> remove;
  // remove const Index that belongs to a
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
  ss << indent << "std::shared_ptr<Tensor<T>> " << label() << "(new Tensor<T>(" << label() << "_index, false));";
  return ss.str();
}


string Tensor::generate_get_block(const string cindent, const string lab, const string tlab, const bool move, const bool noscale) const {
  string lbl = label();
  if (lbl == "proj") lbl = "r";
  size_t found = lbl.find("dagger");
  if (found != string::npos) {
    string tmp(lbl.begin(), lbl.begin()+found);
    lbl = tmp;
  }

  stringstream tt;
  // for scalar. 
  if (index_.empty() && merged_) {
#ifdef debug_tasks
   tt << cindent << "// scalar" << endl;
#endif
  }

  { 

#ifdef debug_tasks // if needed, eg debug
    tt  << cindent << "// tensor label: " << lbl << endl;
#endif
    tt << cindent << "std::unique_ptr<double[]> " << lab << "data = "
                  << tlab << (move ? "->move" : "->get") << "_block(";
    if (found != string::npos) {
      for (auto i = index_.begin(); i != index_.end(); ++i) {
        if (i != index_.begin()) tt << ", ";
        tt << (*i)->str_gen();
      }
      tt << ");" << endl;
    } else { 
      for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
        if (i != index_.rbegin()) tt << ", ";
        tt << (*i)->str_gen();
      }
      tt << ");" << endl;
    }
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


string Tensor::generate_scratch_area(const string cindent, const string lab, const string tensor_lab, const bool zero) const {
  const string lbl = tensor_lab;
  size_t found = label_.find("dagger");

  stringstream ss;
  // using new move/get/put block interface
  ss << cindent << "std::unique_ptr<double[]> " << lab << "data_sorted(new double[" << lbl << "->get_size(";
  if (found != string::npos) {
    for (auto i = index_.begin(); i != index_.end(); ++i) {
      if (i != index_.begin()) ss << ", ";
      ss << (*i)->str_gen();
    }
    ss << ")]);" << endl;
  } else {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
      if (i != index_.rbegin()) ss << ", ";
      ss << (*i)->str_gen();
    }
    ss << ")]);" << endl;
  }
  if (zero) {
    ss << cindent << "std::fill_n(" << lab << "data_sorted.get(), " << lbl << "->get_size(";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
      if (i != index_.rbegin()) ss << ", ";
       ss << (*i)->str_gen();
    }
    ss << "), 0.0);" << endl;
  }
  return ss.str();
}

string Tensor::generate_sort_indices(const string cindent, const string lab, const string tensor_lab, const list<shared_ptr<const Index>>& loop, const bool op) const {
  stringstream ss;
  if (!op) ss << generate_scratch_area(cindent, lab, tensor_lab, false);

  vector<int> map(index_.size());
  // determine mapping
  // first loop indices. order as in loop
  vector<int> done;

  // if trans, transpose here!
  size_t found = label_.find("dagger");
  const bool trans = found != string::npos;
  if (trans && index_.size() & 1) throw logic_error("transposition not possible with 3-index objects");
  if (trans) {
    for (auto i = loop.rbegin(); i != loop.rend(); ++i) {
      int cnt = 0;
      for (auto j = index_.begin(); j != index_.end(); ++j, ++cnt) {
        if ((*i)->identical(*j)) break;
      }
      if (cnt == index_.size()) {
        throw logic_error("should not happen..trans Tensor::generate_sort_indices");
      }
      done.push_back(cnt);
    }
    // then fill out others
    for (int i = 0; i != index_.size(); ++i) {
      if (find(done.begin(), done.end(), i) == done.end())
        done.push_back(i);
    }
  } else {
    for (auto i = loop.rbegin(); i != loop.rend(); ++i) {
      int cnt = 0;
      for (auto j = index_.rbegin(); j != index_.rend(); ++j, ++cnt) {
        if ((*i)->identical(*j)) break;
      }
      if (cnt == index_.size()) {
        throw logic_error("should not happen.. Tensor::generate_sort_indices");
      }
      done.push_back(cnt);
    }
    // then fill out others
    for (int i = 0; i != index_.size(); ++i) {
      if (find(done.begin(), done.end(), i) == done.end())
        done.push_back(i);
    }
  }

  // then write them out.
  ss << cindent << "sort_indices<";
  for (auto i = done.begin(); i != done.end(); ++i)
    ss << *i << ",";

  string target_label = op ? "odata" : lab + "data_sorted";

  ss << (op ? 1 : 0) << ",1," << prefac__(factor_);
  ss << ">(" << lab << "data, " << target_label;
  if (!trans) {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << ", " << (*i)->str_gen() << ".size()";
  } else {
#if 0
    for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
      string tmp = ", " + (*i)->str_gen() + ".size()";
      ++i;
      ss << ", " << (*i)->str_gen() << ".size()" << tmp;
    }
#else
    for (auto i = index_.begin(); i != index_.end(); ++i) {
      ss << ", " << (*i)->str_gen() << ".size()";
    }
#endif
  }
  ss << ");" << endl;
  return ss.str();
}


string Tensor::generate_sort_indices_target(const string cindent, const string lab, const list<shared_ptr<const Index>>& loop,
                                            const shared_ptr<Tensor> a, const shared_ptr<Tensor> b) const {
  stringstream ss;
  vector<int> map(index_.size());
  ss << cindent << "sort_indices<";

  // determine mapping
  // first obtain the ordering of indices from dgemm
  list<shared_ptr<const Index>> source;
  {
    list<shared_ptr<const Index>> aind = a->index();
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


pair<string, string> Tensor::generate_dim(const list<shared_ptr<const Index>>& di) const {
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


string Tensor::generate_active(string indent, const string tag, const int ninptensors, const bool use_blas) const {
  assert(label_.find("Gamma") != string::npos);
  stringstream tt;
  if (!merged_) {
    tt << active()->generate(indent, tag, index());
  } else {
    
#ifdef debug_tasks
    tt << indent <<"// associated with merged" << endl;
#endif

    // add fdata 
    list<shared_ptr<const Index>>& merged = merged_->index();
    // fdata tensor should be last to mirror gamma footer
    tt << indent << "std::unique_ptr<double[]> fdata = in("<< ninptensors-1 << ")->get_block(";
    for (auto j = merged.rbegin(); j != merged.rend(); ++j) {
      if (j != merged.rbegin()) tt << ", ";
      tt << (*j)->str_gen();
    }
    tt << ");" << endl;
   
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
    
    // generate merged and/or rdm 
    tt << active()->generate(indent, tag, index(), merged_->index(), merged_->label(), use_blas);


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


string Tensor::generate_gamma(const int ic, const bool use_blas) const {
  // TODO split generate_gamma function up to header/body/footer once working (too large now)
  assert(label_.find("Gamma") != string::npos);
  stringstream tt;
  
  // determine number of task loops to be separeted, if merged combine
  int nindex;
  list<shared_ptr<const Index>> merged; 
  
  if (merged_) {
    merged = merged_->index();
    nindex = index_.size() + merged.size();
  } else {
    nindex = index_.size();
  }   

  // determine number of tensors
  vector<int> rdmn = active()->required_rdm();
  int ninptensors;
  if (merged_) 
    ninptensors = rdmn.size()+1;
  else
    ninptensors = rdmn.size();


  //////////// gamma header ////////////
  tt << "template <typename T>" << endl;
#ifdef debug_tasks
  tt << "class Task" << ic << " : public Task<T> {" <<  "  // associated with gamma" << endl;
#else
  tt << "class Task" << ic << " : public Task<T> {" << endl;
#endif
  tt << "  protected:" << endl;
  tt << "    class Task_local : public SubTask<" << nindex << "," << ninptensors << ",T> {" << endl;
  tt << "      protected:" << endl;
  tt << "        const std::array<std::shared_ptr<const IndexRange>,3> range_;" << endl;
  tt << endl;

  tt << "        const Index& b(const size_t& i) const { return this->block(i); }" << endl;
  tt << "        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }" << endl;
  tt << "        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }" << endl;
  tt << endl;
  tt << "      public:" << endl;
  tt << "        Task_local(const std::array<const Index," << nindex << ">& block, const std::array<std::shared_ptr<const Tensor<T>>," << ninptensors <<  ">& in, std::shared_ptr<Tensor<T>>& out," << endl;
  tt << "                   std::array<std::shared_ptr<const IndexRange>,3>& ran)" << endl;
  tt << "          : SubTask<" << nindex << "," << ninptensors << ",T>(block, in, out), range_(ran) { }" << endl;
  tt << endl;
  tt << endl;
  tt << "        void compute() override {" << endl;

  //////////// gamma body  ////////////
  string indent ="          ";
  // map indices
  int bcnt = 0;
  for (auto i = index_.begin(); i != index_.end(); ++i, bcnt++) 
    tt << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  if (merged_) {
    for (auto i = merged.begin(); i != merged.end(); ++i, bcnt++) 
      tt << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  }    
  
#ifdef debug_tasks // debug purposes
   tt << "//   std::shared_ptr<Tensor<T> > " << label() << ";" << endl;
   for (auto& i: rdmn)
     tt << "//   std::shared_ptr<Tensor<T> > rdm" << i << ";" << endl;
   if (merged_)
     tt << "//   std::shared_ptr<Tensor<T> > " << merged_->label() << ";" << endl;
   tt << endl;
#endif
  
  // generate gamma get block, true does a move_block
  tt << generate_get_block(indent, "o", "out()", true, true); // first true means move, second true means we don't scale
  if (merged_) {
    if (use_blas && !index_.empty()) tt << generate_scratch_area(indent, "o", "out", true);
  }
  // now generate codes for rdm
  tt << generate_active(indent, "o", ninptensors, use_blas);

  // generate gamma put block
  tt << indent << "out()->put_block(odata";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) 
    tt << ", " << (*i)->str_gen();
  tt << ");" << endl;

  //////////// gamma footer  ////////////
  tt << "        }  " << endl;
  tt << "    };  " << endl;
  tt << "" << endl;
  tt << "    std::vector<std::shared_ptr<Task_local>> subtasks_;" << endl;
  tt << "" << endl;

  tt << "    void compute_() override {" << endl;
  tt << "      for (auto& i : subtasks_) i->compute();" << endl;
  tt << "    }" << endl << endl; 

  tt << "  public:" << endl;
  tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {" << endl;
  tt << "      std::array<std::shared_ptr<const Tensor<T>>," << ninptensors << "> in = {{";

  // write out tensors in increasing order
  for (auto i = 1;  i < ninptensors + 1; ++i) 
    tt << "t[" << i << "]" << (i == ninptensors ? "" : ", ");
  tt << "}};" << endl << endl;

  
  // over original outermost indices
  if (!index_.empty()) {
    tt << "      subtasks_.reserve("; 
    for (auto i =index_.begin(); i != index_.end(); ++i) {
      if (i != index_.begin()) tt << "*";
      tt << (*i)->generate_range() << "->nblock()";
    } 
    if (merged_){
      for (auto i = merged.begin(); i != merged.end(); ++i) {
        tt << "*" << (*i)->generate_range() << "->nblock()";
      } 
    }
    tt << ");" << endl;
  }
  // loops 
  for (auto i = index_.begin(); i != index_.end(); ++i, indent += "  ") 
    tt << indent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  if (merged_) {
    for (auto i = merged.begin(); i != merged.end(); ++i, indent += "  ") 
      tt << indent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  }
  // add subtasks
  tt << indent  << "subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index," << nindex << ">{{";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    if (i != index_.rbegin()) tt << ", ";
    tt << (*i)->str_gen();
  }
  if (merged_) {
    tt << (index_.empty() ? "" : ", ");
    for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
      tt << (*i)->str_gen();
      if (i != --merged.rend()) tt << ", ";
    }
  }
  tt << "}}, in, t[0], range)));" << endl;

  tt << "    };" << endl;
  tt << "    ~Task" << ic << "() {};" << endl;
  tt << "};" << endl << endl;

  return tt.str();
}


