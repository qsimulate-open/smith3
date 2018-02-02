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
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <iomanip>
#include "tensor.h"
#include "constants.h"

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


Tensor::Tensor(const shared_ptr<Active> activ, const list<shared_ptr<const Index>>& in, map<int, int> m) : factor_(1.0), scalar_(""), der_(in), num_map_(m) {
  // scalar quantity..defined on bagel side
  // label
  stringstream ss; ss << "Gamma" << ig; ++ig;
  label_ = ss.str();
  // op
  index_ = activ->index();
  // add extra index, eg ci0
  for (auto& i : in) index_.push_back(i);
  active_ = activ;
}


// print
string Tensor::str() const {
  stringstream ss;
  if (factor_ != 1.0 || !scalar_.empty())
    ss << " " << fixed << setw(4) << setprecision(2) << factor_ << (scalar_.empty() ? "" : " " + scalar_) << " ";
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


// adds all-active tensor to Active_;
void Tensor::merge(shared_ptr<Tensor> a) {
  assert(active_);
  merged_ = a;
  list<list<shared_ptr<const Index>>::iterator> remove;
  // remove const Index that belongs to a
  for (auto& i : a->index()) {
    const int n = i->num();
    for (auto j = index_.begin(); j != index_.end(); ++j) {
      // Careful, this code is driven by numbers
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


string Tensor::constructor_str(const bool diagonal) const {
  stringstream ss;
  string indent = "";
  if (diagonal) {
    indent += "  ";
    ss << "  shared_ptr<Tensor> " << label() << ";" << endl;
    ss << "  if (diagonal) {" << endl;
  }
  ss << indent << "  vector<IndexRange> " << label() << "_index";
  if (index_.empty()) {
    ss << ";" << endl;
  } else {
    ss << " = {";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << (i != index_.rbegin() ? ", " : "") << (*i)->generate();
    ss << "};" << endl;
  }
  ss << indent << "  " << (diagonal ? "" : "auto ") << label() << " = make_shared<Tensor>(" << label() << "_index);";
  if (diagonal)
    ss << endl << "  }";
  return ss.str();
}

string Tensor::constructor_str_ci(const bool diagonal) const {
  stringstream ss;
  string indent = "";
  // ci0 index is now propagated to Itensor_1...5
  for (int j = 0; j != 5; ++j) {
    if (j != 0)
    ss << endl;
    ss << indent << "  vector<IndexRange> " << label() << "_" << j << "_index";
    int no;
    if (j > 3) no = 3;
    else no = j%4;

    if (no) {
      ss << " = {";
      for (int i = 0; i != no; ++i) {
        ss << "active_, active_, ";
      }
      ss << "active_, active_";
      ss << "};" << endl;
    } else {
      ss << ";" << endl;
    }
    ss << indent << "  " << (diagonal ? "" : "auto ") << label() << "_" << j << " = make_shared<Tensor>(" << label() << "_" << j << "_index);";
  }
  return ss.str();
}


string Tensor::generate_get_block(const string cindent, const string lab, const string tlab, const bool move, const bool noscale, const int number, const bool merged, const list<shared_ptr<const Index>>& mergedlist) const {
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
    if (!move)
      tt << cindent << "std::unique_ptr<" << DataType << "[]> " << lab << "data = " << tlab << "->get_block(";
    else
      tt << cindent << "std::unique_ptr<" << DataType << "[]> " << lab << "data(new " << DataType << "[" << tlab << "->get_size(";
    string listind = "";
    if (found != string::npos) {
      int no = 0;
      for (auto i = index_.begin(); i != index_.end(); ++i) {
        if (number==-1 || no==(number+1)) break;
        if (i != index_.begin()) listind += ", ";
        listind += (*i)->str_gen();
        ++no;
      }
    } else {
      int no = 0;
      for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
        if (number==-1 || no==(number+1)) break;
        if (i != index_.rbegin()) listind += ", ";
        listind += (*i)->str_gen();
        ++no;
      }
    }
    if (merged) {
      for (auto i = mergedlist.rbegin(); i != mergedlist.rend(); ++i) {
        if (index_.size()!=0 || i != mergedlist.rbegin())
          listind += ", ";
        listind += (*i)->str_gen();
      }
    }
    if (!move) {
      tt << listind << ");" << endl;
    } else {
      tt << listind << ")]);" << endl;
      tt << cindent << "std::fill_n(" << lab << "data.get(), " << tlab << "->get_size(" << listind << "), 0.0);" << endl;
    }
  }
  if (!scalar_.empty() && !noscale) {
    assert(!move);
    tt << cindent << SCAL << "(";
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      tt << (i != index_.rbegin() ? "*" : "") << (*i)->str_gen() << ".size()";
    // update scalar_ name directly
    tt << ", " << scalar_ << "_, " << lab << "data.get(), 1);" << endl;
  }
  return tt.str();
}


string Tensor::generate_scratch_area(const string cindent, const string lab, const string tensor_lab, const bool zero) const {
  const string lbl = tensor_lab;
  size_t found = label_.find("dagger");

  stringstream ss;
  // using new move/get/put block interface
  ss << cindent << "std::unique_ptr<" << DataType << "[]> " << lab << "data_sorted(new " << DataType << "[" << lbl << "->get_size(";
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
  for (auto& i : done)
    ss << i << ",";

  string target_label = op ? "odata" : lab + "data_sorted";

  ss << (op ? 1 : 0) << ",1," << prefac__(factor_);
  ss << ">(" << lab << "data, " << target_label;
  if (!trans) {
    for (auto i = index_.rbegin(); i != index_.rend(); ++i)
      ss << ", " << (*i)->str_gen() << ".size()";
  } else {
    for (auto& i : index_) {
      ss << ", " << i->str_gen() << ".size()";
    }
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
    // if a is a daggered tensor, we reverse
    if (a->label().find("dagger") != string::npos) aind.reverse();
    for (auto i = aind.rbegin(); i != aind.rend(); ++i) {
      bool found = false;
      for (auto& j : loop)
        if ((*i)->identical(j)) found = true;
      if (!found) source.push_back(*i);
    }
    aind = b->index();
    // if b is a daggered tensor, we reverse
    if (b->label().find("dagger") != string::npos) aind.reverse();
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


string Tensor::generate_active_sources(string indent, const string tag, const int ninptensors, const bool use_blas, const shared_ptr<Tensor> source) const {
  assert(label_.find("Gamma") != string::npos);
  stringstream dd;
  if (!merged_) {
    dd << active()->generate_sources(indent, tag, index());
  } else {

#ifdef debug_tasks
    dd << indent <<"// associated with merged" << endl;
#endif

    // add fdata
    list<shared_ptr<const Index>>& merged = merged_->index();
    // fdata tensor should be last to mirror gamma footer
    dd << indent << "std::unique_ptr<" << DataType << "[]> fdata = in(1)->get_block(";
    for (auto j = merged.rbegin(); j != merged.rend(); ++j) {
      if (j != merged.rbegin()) dd << ", ";
      dd << (*j)->str_gen();
    }
    dd << ");" << endl;

    // generate merged and/or rdm
    dd << active()->generate_sources(indent, tag, index(), merged_->index(), merged_->label(), use_blas);

  }
  return dd.str();
}



string Tensor::generate_active(string indent, const string tag, const int ninptensors, const bool use_blas) const {
  assert(label_.find("Gamma") != string::npos);
  stringstream dd;
  if (!merged_) {
    dd << active()->generate(indent, tag, index());
  } else {

#ifdef debug_tasks
    dd << indent <<"// associated with merged" << endl;
#endif

    // add fdata
    list<shared_ptr<const Index>>& merged = merged_->index();
    // fdata tensor should be last to mirror gamma footer
    dd << indent << "std::unique_ptr<" << DataType << "[]> fdata = in("<< ninptensors-1 << ")->get_block(";
    for (auto j = merged.rbegin(); j != merged.rend(); ++j) {
      if (j != merged.rbegin()) dd << ", ";
      dd << (*j)->str_gen();
    }
    dd << ");" << endl;

    if (use_blas) {
      dd << indent << "std::unique_ptr<" << DataType << "[]> fdata_sorted(new " << DataType << "["<< merged_->label() << "->get_size(fhash)]);" << endl;

      // make sort_indices for merged op
      vector<int> done;
      dd << indent << "sort_indices<";
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
        dd << i << ",";

      dd << "0,1,1,1";

      // add source data dimensions
      dd << ">(fdata, fdata_sorted, " ;
      for (auto iter = merged.rbegin(); iter != merged.rend(); ++iter) {
        if (iter != merged.rbegin()) dd << ", ";
          dd << (*iter)->str_gen() << ".size()";
      }
      dd << ");" << endl;
    }

    // generate merged and/or rdm
    dd << active()->generate(indent, tag, index(), merged_->index(), merged_->label(), use_blas);

  }
  return dd.str();
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

OutStream Tensor::generate_gamma_sources(const int ic, const bool use_blas, const bool der, const shared_ptr<Tensor> source, const list<shared_ptr<const Index>> di) const {
  assert(label_.find("Gamma") != string::npos);
  OutStream out;

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
  vector<string> rdmn = active()->required_rdm();
  int ninptensors;
  if (merged_) {
    ninptensors = rdmn.size()+1;
  } else {
    ninptensors = rdmn.size();
  }

  out << generate_gamma_header_sources(ic, use_blas, der, nindex);
  out << generate_gamma_body_sources(ic, use_blas, der, nindex, ninptensors, merged, source, di);
  out << generate_gamma_footer_sources(ic, use_blas, der, nindex, ninptensors, merged);

  return out;
}


OutStream Tensor::generate_gamma(const int ic, const bool use_blas, const bool der) const {
  assert(label_.find("Gamma") != string::npos);
  OutStream out;

  // if it is derivative code, do not generate anything
  if (der) return out;

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
  vector<string> rdmn = active()->required_rdm();
  int ninptensors;
  if (merged_) {
    ninptensors = rdmn.size()+1;
  } else {
    ninptensors = rdmn.size();
  }

  out << generate_gamma_header(ic, use_blas, der, nindex, ninptensors);
  out << generate_gamma_body(ic, use_blas, der, nindex, ninptensors, merged);
  out << generate_gamma_footer(ic, use_blas, der, nindex, ninptensors, merged);

  return out;
}

OutStream Tensor::generate_gamma_header_sources(const int ic, const bool use_blas, const bool der, const int nindex) const {
  OutStream out;

#ifdef debug_tasks
  out.tt << "class Task" << ic << " : public Task {" <<  "  // merged with gamma" << endl;
#else
  out.tt << "class Task" << ic << " : public Task {" << endl;
#endif
  out.tt << "  protected:" << endl;
  out.tt << "    std::array<std::shared_ptr<Tensor>,5> out_;" << endl;
  out.tt << "    std::array<std::shared_ptr<const Tensor>," << (merged_? "2":"1") << "> in_;" << endl;
  out.tt << "    class Task_local : public SubTask_Merged<" << nindex << "," << (merged_? "2":"1") << ",5> {" << endl;    // five here is output tensors
  out.tt << "      protected:" << endl;
  out.tt << "        const std::array<std::shared_ptr<const IndexRange>,3> range_;" << endl;
  out.tt << endl;

  out.tt << "        const Index& b(const size_t& i) const { return this->block(i); }" << endl;
  out.tt << "        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }" << endl;
  out.tt << "        std::shared_ptr<Tensor> out(const size_t& i) { return this->out_tensor(i); }" << endl;
  out.tt << endl;
  out.tt << "      public:" << endl;
  out.tt << "        Task_local(const std::array<const Index," << nindex << ">& block, const std::array<std::shared_ptr<const Tensor>," << (merged_? "2":"1") <<  ">& in," << endl;
  out.tt << "                   const std::array<std::shared_ptr<Tensor>,5>& out,";
  out.tt << " std::array<std::shared_ptr<const IndexRange>,3>& ran)" << endl;
  out.tt << "          : SubTask_Merged<" << nindex << "," << (merged_? "2":"1") << ",5>(block, in, out), range_(ran) { }" << endl;
  out.tt << endl;
  out.tt << endl;
  out.tt << "        void compute() override;" << endl;

  out.dd << "void Task" << ic << "::Task_local::compute() {" << endl;

  return out;
}


OutStream Tensor::generate_gamma_body_sources(const int ic, const bool use_blas, const bool der, const int nindex, const int ninptensors, const list<shared_ptr<const Index>>& merged, const shared_ptr<Tensor> source, const list<shared_ptr<const Index>> di) const {
  OutStream out;

  string indent ="  ";
  int bcnt = 0;
  for (auto i = index_.rbegin(); i != index_.rend(); ++i, bcnt++) {
    out.dd << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  }
  if (merged_) {
    for (auto i = merged.rbegin(); i != merged.rend(); ++i, bcnt++)
      out.dd << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  }

  // generate gamma get block, true does a move_block

  out.dd << source->generate_get_block(indent, "i0", "in(0)");
  out.dd << source->generate_sort_indices(indent, "i0", "in(0)", di) << endl;

#ifdef debug_tasks // if needed, eg debug
  out.dd << indent << "// tensor label (calculated on-the-fly): " << label() << endl;
#endif
  // TODO currently not considering blas
  if (merged_) {
    if (use_blas && !index_.empty()) out.dd << generate_scratch_area(indent, "o", "out", true);
  }
  // now generate codes for rdm
  out.dd << generate_active_sources(indent, "o", ninptensors, use_blas, source);
  out.dd << "}" << endl << endl << endl;

  return out;
}


OutStream Tensor::generate_gamma_footer_sources(const int ic, const bool use_blas, const bool der, const int nindex, const int ninptensors, const list<shared_ptr<const Index>>& merged) const {
  OutStream out;

  out.tt << "    };" << endl;
  out.tt << "" << endl;
  out.tt << "    std::vector<std::shared_ptr<Task_local>> subtasks_;" << endl;
  out.tt << "" << endl;

  out.tt << "    void compute_() override {" << endl;
  out.tt << "      if (!out_[0]->allocated())" << endl;
  out.tt << "        out_[0]->allocate();" << endl;
  out.tt << "      if (!out_[1]->allocated())" << endl;
  out.tt << "        out_[1]->allocate();" << endl;
  out.tt << "      if (!out_[2]->allocated())" << endl;
  out.tt << "        out_[2]->allocate();" << endl;
  out.tt << "      if (!out_[3]->allocated())" << endl;
  out.tt << "        out_[3]->allocate();" << endl;
  out.tt << "      if (!out_[4]->allocated())" << endl;
  out.tt << "        out_[4]->allocate();" << endl;
  out.tt << "      for (auto& i : subtasks_) i->compute();" << endl;
  out.tt << "    }" << endl << endl;

  out.tt << "  public:" << endl;
  out.tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);" << endl;

  out.cc << "Task" << ic << "::Task" << ic << "(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {" << endl;
  out.cc << "  array<shared_ptr<const Tensor>," << (merged_? "2" : "1") << "> in = {{";

  // write out tensors in increasing order
  for (auto i = 1;  i < (merged_? 2 : 1) + 1; ++i)
    out.cc << "t[" << i+4 << "]" << (i == (merged_? 2 : 1) ? "" : ", ");
  out.cc << "}};" << endl;

  out.cc << "  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};" << endl;
  out.cc << "  in_ = in;" << endl;
  out.cc << "  out_ = out;" << endl << endl;

  // over original outermost indices
  if (!index_.empty()) {
    out.cc << "  subtasks_.reserve(";
    for (auto i =index_.rbegin(); i != index_.rend(); ++i) {
      if (i != index_.rbegin()) out.cc << "*";
      out.cc << (*i)->generate_range() << "->nblock()";
    }
    if (merged_){
      for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
        out.cc << "*" << (*i)->generate_range() << "->nblock()";
      }
    }
    out.cc << ");" << endl;
  }
  // loops
  string cindent = "  ";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i, cindent += "  ")
    out.cc << cindent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  if (merged_)
    for (auto i = merged.rbegin(); i != merged.rend(); ++i, cindent += "  ")
      out.cc << cindent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  // parallel if
  string listind = "";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    if (i != index_.rbegin()) listind += ", ";
    listind += (*i)->str_gen();
  }
  // the input t decides whether it is local or not --- wrong when merged
  out.cc << cindent << "if (t[5]->is_local("<< listind << "))" << endl;
  cindent += "  ";
  // add subtasks
  out.cc << cindent  << "subtasks_.push_back(make_shared<Task_local>(array<const Index," << nindex << ">{{" << listind;
  if (merged_) {
    out.cc << (index_.empty() ? "" : ", ");
    for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
      out.cc << (*i)->str_gen();
      if (i != --merged.rend()) out.cc << ", ";
    }
  }
  out.cc << "}}, in, out, range));" << endl;
  out.cc << "}" << endl << endl << endl;

  out.tt << "    ~Task" << ic << "() {}" << endl;
  out.tt << "};" << endl << endl;


  return out;
}


OutStream Tensor::generate_gamma_header(const int ic, const bool use_blas, const bool der, const int nindex, const int ninptensors) const {
  OutStream out;

#ifdef debug_tasks
  out.tt << "class Task" << ic << " : public Task {" <<  "  // associated with gamma" << endl;
#else
  out.tt << "class Task" << ic << " : public Task {" << endl;
#endif
  out.tt << "  protected:" << endl;
  out.tt << "    std::shared_ptr<Tensor> out_;" << endl;
  out.tt << "    std::array<std::shared_ptr<const Tensor>," << ninptensors << "> in_;" << endl;
  out.tt << "    class Task_local : public SubTask<" << nindex << "," << ninptensors << "> {" << endl;
  out.tt << "      protected:" << endl;
  out.tt << "        const std::array<std::shared_ptr<const IndexRange>,3> range_;" << endl;
  out.tt << endl;

  out.tt << "        const Index& b(const size_t& i) const { return this->block(i); }" << endl;
  out.tt << "        std::shared_ptr<const Tensor> in(const size_t& i) const { return this->in_tensor(i); }" << endl;
  out.tt << "        std::shared_ptr<Tensor> out() { return this->out_tensor(); }" << endl;
  out.tt << endl;
  out.tt << "      public:" << endl;
  out.tt << "        Task_local(const std::array<const Index," << nindex << ">& block, const std::array<std::shared_ptr<const Tensor>," << ninptensors <<  ">& in, std::shared_ptr<Tensor>& out," << endl;
  out.tt << "                   std::array<std::shared_ptr<const IndexRange>,3>& ran)" << endl;
  out.tt << "          : SubTask<" << nindex << "," << ninptensors << ">(block, in, out), range_(ran) { }" << endl;
  out.tt << endl;
  out.tt << endl;
  out.tt << "        void compute() override;" << endl;

  out.dd << "void Task" << ic << "::Task_local::compute() {" << endl;

  return out;
}

OutStream Tensor::generate_gamma_body(const int ic, const bool use_blas, const bool der, const int nindex, const int ninptensors, const list<shared_ptr<const Index>>& merged) const {
  OutStream out;

  string indent ="  ";
  // map indices
  int bcnt = 0;
  for (auto i = index_.rbegin(); i != index_.rend(); ++i, bcnt++)
    out.dd << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  if (merged_) {
    for (auto i = merged.rbegin(); i != merged.rend(); ++i, bcnt++)
      out.dd << indent << "const Index " << (*i)->str_gen() << " = b(" << bcnt << ");" << endl;
  }

  // generate gamma get block, true does a move_block
  out.dd << generate_get_block(indent, "o", "out()", /*move=*/true, /*noscale=*/true);
  if (merged_) {
    if (use_blas && !index_.empty()) out.dd << generate_scratch_area(indent, "o", "out", true);
  }
  // now generate codes for rdm
  out.dd << generate_active(indent, "o", ninptensors, use_blas);

  // generate gamma put block
  out.dd << indent << "out()->add_block(odata";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i)
    out.dd << ", " << (*i)->str_gen();
  out.dd << ");" << endl;
  out.dd << "}" << endl << endl << endl;

  return out;
}

OutStream Tensor::generate_gamma_footer(const int ic, const bool use_blas, const bool der, const int nindex, const int ninptensors, const list<shared_ptr<const Index>>& merged) const {
  OutStream out;

  out.tt << "    };" << endl;
  out.tt << "" << endl;
  out.tt << "    std::vector<std::shared_ptr<Task_local>> subtasks_;" << endl;
  out.tt << "" << endl;

  out.tt << "    void compute_() override {" << endl;
  out.tt << "      if (!out_->allocated())" << endl;
  out.tt << "        out_->allocate();" << endl;
  out.tt << "      for (auto& i : subtasks_) i->compute();" << endl;
  out.tt << "    }" << endl << endl;

  out.tt << "  public:" << endl;
  out.tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,3> range);" << endl;

  out.cc << "Task" << ic << "::Task" << ic << "(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {" << endl;
  out.cc << "  array<shared_ptr<const Tensor>," << ninptensors << "> in = {{";

  // write out tensors in increasing order
  for (auto i = 1;  i < ninptensors + 1; ++i)
    out.cc << "t[" << i << "]" << (i == ninptensors ? "" : ", ");
  out.cc << "}};" << endl;

  out.cc << "  out_ = t[0];" << endl;
  out.cc << "  in_ = in;" << endl << endl;

  // over original outermost indices
  if (!index_.empty()) {
    out.cc << "  subtasks_.reserve(";
    for (auto i =index_.rbegin(); i != index_.rend(); ++i) {
      if (i != index_.rbegin()) out.cc << "*";
      out.cc << (*i)->generate_range() << "->nblock()";
    }
    if (merged_){
      for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
        out.cc << "*" << (*i)->generate_range() << "->nblock()";
      }
    }
    out.cc << ");" << endl;
  }
  // loops
  string cindent = "  ";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i, cindent += "  ")
    out.cc << cindent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  if (merged_) {
    for (auto i = merged.rbegin(); i != merged.rend(); ++i, cindent += "  ")
      out.cc << cindent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  }
  // parallel if
  string listind = "";
  for (auto i = index_.rbegin(); i != index_.rend(); ++i) {
    if (i != index_.rbegin()) listind += ", ";
    listind += (*i)->str_gen();
  }
  out.cc << cindent << "if (t[0]->is_local("<< listind << "))" << endl;
  cindent += "  ";
  // add subtasks
  out.cc << cindent  << "subtasks_.push_back(make_shared<Task_local>(array<const Index," << nindex << ">{{" << listind;
  if (merged_) {
    out.cc << (index_.empty() ? "" : ", ");
    for (auto i = merged.rbegin(); i != merged.rend(); ++i) {
      out.cc << (*i)->str_gen();
      if (i != --merged.rend()) out.cc << ", ";
    }
  }
  out.cc << "}}, in, t[0], range));" << endl;
  out.cc << "}" << endl << endl << endl;

  out.tt << "    ~Task" << ic << "() {}" << endl;
  out.tt << "};" << endl << endl;


  return out;
}
