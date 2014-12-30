//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: dedci.cc
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


#include "dedci.h"
#include <algorithm>

using namespace std;
using namespace smith;




// local functions... (not a good practice...) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
static string merge__(vector<string> array) {
  stringstream ss;
  vector<string> done;
  for (auto& label : array) {
    size_t found = label.find("dagger");
    if (found != string::npos) {
      string tmp(label.begin(), label.begin() + found);
      label = tmp;
    }
    if (find(done.begin(), done.end(), label) != done.end()) continue;
    done.push_back(label);
    if (label == "f1" || label == "v2" || label == "h1") label = label + "_";
    if (label == "dci") label = "rdm0deriv_";
    ss << (label != array.front() ? ", " : "") << ((label == "proj") ? "deci" : label);
  }
  return ss.str();
}
static string merge__(list<string> array) { return merge__(vector<string>(array.begin(), array.end())); }
// local functions... (not a good practice...) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


OutStream Dedci::create_target(const int i) const {
  OutStream out;

  out.tt << "class Task" << i << " : public Task {" << endl;
  out.tt << "  protected:" << endl;
  out.tt << "    std::shared_ptr<Tensor> dec_;" << endl;
  out.tt << "    IndexRange closed_;" << endl;
  out.tt << "    IndexRange active_;" << endl;
  out.tt << "    IndexRange virt_;" << endl;
  out.tt << "    IndexRange ci_;" << endl;
  out.tt << "" << endl;
  out.tt << "    void compute_() {" << endl;
  out.tt << "      dec_->zero();" << endl;
  out.tt << "    }" << endl;
  out.tt << "" << endl;
  out.tt << "  public:" << endl;
  out.tt << "    Task" << i << "(std::vector<std::shared_ptr<Tensor>> t);" << endl;

  out.cc << "Task" << i << "::Task" << i << "(vector<shared_ptr<Tensor>> t) {" << endl;
  out.cc << "  dec_ =  t[0];" << endl;
  out.cc << "}" << endl << endl << endl;

  out.tt << "    ~Task" << i << "() {}" << endl;
  out.tt << "};" << endl << endl;

  out.ee << "  auto dedci_ = make_shared<Queue>();" << endl;
  out.ee << "  vector<shared_ptr<Tensor>> tensor" << i << " = {deci};" << endl;
  out.ee << "  auto task" << i << " = make_shared<Task" << i << ">(tensor" << i << ");" << endl;
  out.ee << "  dedci_->add_task(task" << i << ");" << endl << endl;

  return out;
}


shared_ptr<Tensor> Dedci::create_tensor(list<shared_ptr<const Index>> dm) const {
 return make_shared<Tensor>(1.0, "deci", dm);
}


OutStream Dedci::generate_task(const int ip, const int ic, const vector<string> op, const string scalar, const int iz, bool der) const {
  OutStream out;
  out.ee << "  vector<shared_ptr<Tensor>> tensor" << ic << " = {" << merge__(op) << "};" << endl;
  out.ee << "  auto task" << ic << " = make_shared<Task" << ic << ">(tensor" << ic << ", cindex" << (scalar.empty() ? "" : ", this->e0_") << ");" << endl;
  if (parent_) {
    assert(parent_->parent());
    out.ee << "  task" << ip << "->add_dep(task" << ic << ");" << endl;
    out.ee << "  task" << ic << "->add_dep(task" << iz << ");" << endl;
  } else {
    assert(depth() == 0);
    out.ee << "  task" << ic << "->add_dep(task" << iz << ");" << endl;
  }
  out.ee << "  dedci_->add_task(task" << ic << ");" << endl;
  out.ee << endl;
  return out;
}


OutStream Dedci::generate_compute_header(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors, const bool no_outside) const {
  const int ninptensors = tensors.size()-1;

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  const int nindex = ti.size();
  OutStream out;
  out.tt << "class Task" << ic << " : public Task {" << endl;
  out.tt << "  protected:" << endl;
  // if index is empty give dummy arg
  out.tt << "    class Task_local : public SubTask<" << (ti.empty() ? 1 : nindex) << "," << ninptensors << "> {" << endl;
  out.tt << "      protected:" << endl;
  out.tt << "        const std::array<std::shared_ptr<const IndexRange>,4> range_;" << endl << endl;

  out.tt << "        const Index& b(const size_t& i) const { return this->block(i); }" << endl;
  out.tt << "        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }" << endl;
  out.tt << "        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }" << endl;
  if (need_e0) out.tt << endl << "        double e0_;" << endl;
  out.tt << endl;
  out.tt << "      public:" << endl;
  // if index is empty use dummy index 1 to subtask
  if (ti.empty()) {
    out.tt << "        Task_local(const std::array<std::shared_ptr<const Tensor>," << ninptensors <<  ">& in, std::shared_ptr<Tensor>& out," << endl;
    out.tt << "                   std::array<std::shared_ptr<const IndexRange>,4>& ran" << (need_e0 ? ", const double e" : "") << ")" << endl;
    out.tt << "          : SubTask<1," << ninptensors << ">(std::array<const Index, 1>(), in, out), range_(ran)" << (need_e0 ? ", e0_(e)" : "") << " { }" << endl;
  } else {
    out.tt << "        Task_local(const std::array<const Index," << nindex << ">& block, const std::array<std::shared_ptr<const Tensor>," << ninptensors <<  ">& in, std::shared_ptr<Tensor>& out," << endl;
    out.tt << "                   std::array<std::shared_ptr<const IndexRange>,4>& ran" << (need_e0 ? ", const double e" : "") << ")" << endl;
    out.tt << "          : SubTask<" << nindex << "," << ninptensors << ">(block, in, out), range_(ran)" << (need_e0 ? ", e0_(e)" : "") << " { }" << endl;
  }
  out.tt << endl;
  out.tt << endl;
  out.tt << "        void compute() override;" << endl;

  out.dd << "void Task" << ic << "::Task_local::compute() {" << endl;

  if (!no_outside) {
    list<shared_ptr<const Index>> ti_copy = ti;
    if (depth() == 0) {
      if (ti.size() > 1) {
        for (auto i = ti_copy.begin(), j = ++ti_copy.begin(); i != ti_copy.end(); ++i, ++i, ++j, ++j)
          swap(*i, *j);
      }
    }

    int cnt = 0;
    for (auto i = ti_copy.rbegin(); i != ti_copy.rend(); ++i)
      out.dd << "  const Index " << (*i)->str_gen() << " = b(" << cnt++ << ");" << endl;
    out.dd << endl;
  }

  return out;
}


OutStream Dedci::generate_compute_footer(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors) const {
  const int ninptensors = tensors.size()-1;
  assert(ninptensors > 0);
  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  OutStream out;
  out.dd << "}" << endl << endl << endl;

  out.tt << "    };" << endl;
  out.tt << "" << endl;
  out.tt << "    std::vector<std::shared_ptr<Task_local>> subtasks_;" << endl;
  out.tt << "" << endl;

  out.tt << "    void compute_() override {" << endl;
  out.tt << "      for (auto& i : subtasks_) {" << endl;
  out.tt << "        i->compute();" << endl;
  out.tt << "      }" << endl;
  out.tt << "    }" << endl << endl;

  out.tt << "  public:" << endl;
  out.tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor>> t, std::array<std::shared_ptr<const IndexRange>,4> range" << (need_e0 ? ", double e" : "" ) <<  ");" << endl;

  out.cc << "Task" << ic << "::Task" << ic << "(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range" << (need_e0 ? ", double e" : "" ) <<  ") {" << endl;
  out.cc << "  array<shared_ptr<const Tensor>," << ninptensors << "> in = {{";
  for (auto i = 1; i < ninptensors + 1; ++i)
    out.cc << "t[" << i << "]" << (i < ninptensors ? ", " : "");
  out.cc << "}};" << endl << endl;

  // over original outermost indices
  if (!ti.empty()) {
    out.cc << "  subtasks_.reserve(";
    for (auto i = ti.begin(); i != ti.end(); ++i) {
      if (i != ti.begin()) out.cc << "*";
      out.cc << (*i)->generate_range() << "->nblock()";
    }
    out.cc << ");" << endl;
  }
  // loops
  string indent = "  ";
  for (auto i = ti.begin(); i != ti.end(); ++i, indent += "  ")
    out.cc << indent << "for (auto& " << (*i)->str_gen() << " : *" << (*i)->generate_range() << ")" << endl;
  // add subtasks
  if (!ti.empty()) {
    out.cc << indent  << "subtasks_.push_back(make_shared<Task_local>(std::array<const Index," << ti.size() << ">{{";
    for (auto i = ti.rbegin(); i != ti.rend(); ++i) {
      if (i != ti.rbegin()) out.cc << ", ";
      out.cc << (*i)->str_gen();
    }
    out.cc << "}}, in, t[0], range" << (need_e0 ? ", e" : "") << "));" << endl;
  } else {
    out.cc << indent  << "subtasks_.push_back(make_shared<Task_local>(in, t[0], range" << (need_e0 ? ", e" : "") << "));" << endl;
  }
  out.cc << "}" << endl << endl << endl;

  out.tt << "    ~Task" << ic << "() {}" << endl;
  out.tt << "};" << endl << endl;
  return out;
}


OutStream Dedci::generate_bc(const shared_ptr<BinaryContraction> i) const {
  OutStream out;

  if (depth() != 0) {
    const string bindent = "  ";
    string dindent = bindent;

    out.dd << target_->generate_get_block(dindent, "o", "out()", true);
    out.dd << target_->generate_scratch_area(dindent, "o", "out()", true); // true means zero-out

    list<shared_ptr<const Index>> ti = depth() != 0 ? i->target_indices() : i->tensor()->index();

    // inner loop (where similar indices in dgemm tensors are summed over) will show up here
    // but only if outer loop is not empty
    list<shared_ptr<const Index>> di = i->loop_indices();
//  di.reverse();

    vector<string> close2;
    if (ti.size() != 0) {
      out.dd << endl;
      for (auto iter = di.rbegin(); iter != di.rend(); ++iter, dindent += "  ") {
        string index = (*iter)->str_gen();
        out.dd << dindent << "for (auto& " << index << " : *" << (*iter)->generate_range("_") << ") {" << endl;
        close2.push_back(dindent + "}");
      }
    } else {
      int cnt = 0;
      for (auto k = di.begin(); k != di.end(); ++k, cnt++) out.dd << dindent << "const Index " <<  (*k)->str_gen() << " = b(" << cnt << ");" << endl;
      out.dd << endl;
    }

    // retrieving tensor_
    out.dd << i->tensor()->generate_get_block(dindent, "i0", "in(0)");
    out.dd << i->tensor()->generate_sort_indices(dindent, "i0", "in(0)", di) << endl;
    // retrieving subtree_
    out.dd << i->next_target()->generate_get_block(dindent, "i1", "in(1)");
    out.dd << i->next_target()->generate_sort_indices(dindent, "i1", "in(1)", di) << endl;

    // call dgemm
    {
      pair<string, string> t0 = i->tensor()->generate_dim(di);
      pair<string, string> t1 = i->next_target()->generate_dim(di);
      if (t0.first != "" || t1.first != "") {
        out.dd << dindent << "dgemm_(\"T\", \"N\", ";
        string tt0 = t0.first == "" ? "1" : t0.first;
        string tt1 = t1.first == "" ? "1" : t1.first;
        string ss0 = t1.second== "" ? "1" : t1.second;
        out.dd << tt0 << ", " << tt1 << ", " << ss0 << "," << endl;
        out.dd << dindent << "       1.0, i0data_sorted, " << ss0 << ", i1data_sorted, " << ss0 << "," << endl
           << dindent << "       1.0, odata_sorted, " << tt0;
        out.dd << ");" << endl;
      } else {
        string ss0 = t1.second== "" ? "1" : t1.second;
        out.dd << dindent << "odata_sorted[0] += ddot_(" << ss0 << ", i0data_sorted, 1, i1data_sorted, 1);" << endl;
      }
    }

    if (ti.size() != 0) {
      for (auto iter = close2.rbegin(); iter != close2.rend(); ++iter)
        out.dd << *iter << endl;
      out.dd << endl;
    }
    // Inner loop ends here

    // sort buffer
    {
      out.dd << i->target()->generate_sort_indices_target(bindent, "o", di, i->tensor(), i->next_target());
    }
    // put buffer
    {
      string label = target_->label();
      // new interface requires indices for put_block
      out.dd << bindent << "out()->put_block(odata";
      list<shared_ptr<const Index>> ti = depth() != 0 ? i->target_indices() : i->tensor()->index();
      for (auto i = ti.rbegin(); i != ti.rend(); ++i)
        out.dd << ", " << (*i)->str_gen();
      out.dd << ");" << endl;
    }

  } else { // bc depth = 0
    // depth should only be zero and here in residual tree
    throw logic_error("shouldn't happen in Dedci::generate_bc");
  }

  return out;
}









