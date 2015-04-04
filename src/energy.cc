//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: energy.cc
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


#include "energy.h"
#include "constants.h"

using namespace std;
using namespace smith;


OutStream Energy::generate_task(const int ip, const int ic, const vector<string> op, const string scalar, const int, bool, bool) const {
  OutStream out;
  out.ee << "  auto tensor" << ic << " = vector<shared_ptr<Tensor>>{" << merge__(op) << "};" << endl;
  out.ee << "  auto task" << ic << " = make_shared<Task" << ic << ">(tensor" << ic << ", pindex" << (scalar.empty() ? "" : ", this->e0_") << ");" << endl;

  if (parent_) {
    if (ip != ic)
      out.ee << "  task" << ip << "->add_dep(task" << ic << ");" << endl;
  } else {
    assert(depth() == 0);
  }
  out.ee << "  " << label_ << "q->add_task(task" << ic << ");" << endl;
  out.ee << endl;
  return out;
}


OutStream Energy::generate_compute_header(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors, const bool no_outside) const {
  vector<string> labels;
  for (auto i = ++tensors.begin(); i != tensors.end(); ++i)
    labels.push_back((*i)->label());
  const int ninptensors = count_distinct_tensors__(labels);

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  const int nindex = ti.size();
  OutStream out;
  out.tt << "class Task" << ic << " : public AccTask {" << endl;
  out.tt << "  protected:" << endl;
  // if index is empty give dummy arg
  out.tt << "    class Task_local : public SubTask<" << (ti.empty() ? 1 : nindex) << "," << ninptensors << "> {" << endl;
  out.tt << "      protected:" << endl;
  out.tt << "        const std::array<std::shared_ptr<const IndexRange>,3> range_;" << endl << endl;

  out.tt << "        const Index& b(const size_t& i) const { return this->block(i); }" << endl;
  out.tt << "        const std::shared_ptr<const Tensor>& in(const size_t& i) const { return this->in_tensor(i); }" << endl;
  out.tt << "        const std::shared_ptr<Tensor>& out() const { return this->out_tensor(); }" << endl;
  out.tt << "        double target_;" << endl;
  if (need_e0)  out.tt << "        double e0_;" << endl;
  out.tt << endl;
  out.tt << "      public:" << endl;
  // if index is empty use dummy index 1 to subtask
  if (ti.empty()) {
    out.tt << "        Task_local(const std::array<std::shared_ptr<const Tensor>," << ninptensors <<  ">& in, std::shared_ptr<Tensor>& out," << endl;
    out.tt << "                   std::array<std::shared_ptr<const IndexRange>,3>& ran" << (need_e0 ? ", const double e" : "") << ")" << endl;
    out.tt << "          : SubTask<1," << ninptensors << ">(std::array<const Index, 1>(), in, out), range_(ran)" << (need_e0 ? ", e0_(e)" : "") << " { }" << endl;
  } else {
    out.tt << "        Task_local(const std::array<const Index," << nindex << ">& block, const std::array<std::shared_ptr<const Tensor>," << ninptensors <<  ">& in, std::shared_ptr<Tensor>& out," << endl;
    out.tt << "                   std::array<std::shared_ptr<const IndexRange>,3>& ran" << (need_e0 ? ", const double e" : "") << ")" << endl;
    out.tt << "          : SubTask<" << nindex << "," << ninptensors << ">(block, in, out), range_(ran)" << (need_e0 ? ", e0_(e)" : "") << " { }" << endl;
  }
  out.tt << endl;
  out.tt << "        double target() const { return target_; }" << endl;
  out.tt << endl;
  out.tt << "        void compute() override;" << endl;

  out.dd << "void Task" << ic << "::Task_local::compute() {" << endl;
  out.dd << "  target_ = 0.0;" << endl;

  if (!no_outside) {
    list<shared_ptr<const Index>> ti_copy = ti;
    if (depth() == 0) {
      for (auto i = ti_copy.begin(), j = ++ti_copy.begin(); i != ti_copy.end(); ++i, ++i, ++j, ++j)
        swap(*i, *j);
    }

    int cnt = 0;
    for (auto i = ti_copy.rbegin(); i != ti_copy.rend(); ++i)
      out.dd << "  const Index " << (*i)->str_gen() << " = b(" << cnt++ << ");" << endl;
    out.dd << endl;
  }

  return out;
}


OutStream Energy::generate_compute_footer(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors) const {
  vector<string> labels;
  for (auto i = ++tensors.begin(); i != tensors.end(); ++i)
    labels.push_back((*i)->label());
  const int ninptensors = count_distinct_tensors__(labels);
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
  out.tt << "      this->target_ = 0.0;" << endl;
  out.tt << "      for (auto& i : subtasks_) {" << endl;
  out.tt << "        i->compute();" << endl;
  out.tt << "        this->target_ += i->target();" << endl;
  out.tt << "      }" << endl;
  out.tt << "    }" << endl << endl;

  out.tt << "  public:" << endl;
  out.tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor>> t,  std::array<std::shared_ptr<const IndexRange>,3> range" << (need_e0 ? ", const double e" : "") << ");" << endl;

  out.cc << "Task" << ic << "::Task" << ic << "(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range" << (need_e0 ? ", const double e" : "") << ") {" << endl;
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
    out.cc << indent  << "subtasks_.push_back(make_shared<Task_local>(array<const Index," << ti.size() << ">{{";
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


OutStream Energy::generate_bc(const shared_ptr<BinaryContraction> i) const {
  OutStream out;

  if (depth() != 0) {
    const string bindent = "  ";
    string dindent = bindent;
    // skip if energy tree depth is 1
    if (depth() != 1) {
      out.dd << target_->generate_get_block(dindent, "o", "out()", true);
      out.dd << target_->generate_scratch_area(dindent, "o", "out()", true); // true means zero-out
    }

    list<shared_ptr<const Index>> ti = depth() != 0 ? i->target_indices() : i->tensor()->index();

    // inner loop will show up here
    // but only if outer loop is not empty
    list<shared_ptr<const Index>> di = i->loop_indices();
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
    string inlabel("in("); inlabel += (same_tensor__(i->tensor()->label(), i->next_target()->label()) ? "0)" : "1)");
    out.dd << i->next_target()->generate_get_block(dindent, "i1", inlabel);
    out.dd << i->next_target()->generate_sort_indices(dindent, "i1", inlabel, di) << endl;

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
        if (depth() != 1) {
          string ss0 = t1.second== "" ? "1" : t1.second;
          out.dd << dindent << "odata_sorted[0] += ddot_(" << ss0 << ", i0data_sorted, 1, i1data_sorted, 1);" << endl;
        } else {
          string ss0 = t1.second== "" ? "1" : t1.second;
          out.dd << dindent << "target_ += ddot_(" << ss0 << ", i0data_sorted, 1, i1data_sorted, 1);" << endl;
        }
      }
    }

    if (ti.size() != 0) {
      for (auto iter = close2.rbegin(); iter != close2.rend(); ++iter)
        out.dd << *iter << endl;
      out.dd << endl;
    }
    // Inner loop ends here

    // skip if energy tree depth is 1
    if (depth() != 1) {
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
    }
  } else {
    // depth should not equal 0 in energy tree
    throw logic_error("shouldn't happen in Energy::generate_bc");
  }


  return out;
}

