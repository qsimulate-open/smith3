//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: residual.cc
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
#include "residual.h"

using namespace std;
using namespace smith;


OutStream Residual::create_target(const int i) const {
  OutStream out;

  out.tt << "class Task" << i << " : public Task {" << endl;
  out.tt << "  protected:" << endl;
  out.tt << "    std::array<std::shared_ptr<Tensor>,1> tensor_;" << endl;
  out.tt << "    const bool reset_;" << endl;
  out.tt << "" << endl;
  out.tt << "    void compute_() {" << endl;
  out.tt << "      if (reset_) " << target_name__(label_) << "_->zero();" << endl;
  out.tt << "    }" << endl;
  out.tt << "" << endl;
  out.tt << "  public:" << endl;
  out.tt << "    Task" << i << "(std::array<std::shared_ptr<Tensor>,1> t, const bool reset) : tensor_(t), reset_(reset) { }" << endl;
  out.tt << "};" << endl << endl;

  out.ee << "  auto " << label_ << "q = make_shared<Queue>();" << endl;
  out.ee << "  auto tensor" << i << " = array<shared_ptr<Tensor>,1>{{" << target_name__(label_) << "}};" << endl;
  out.ee << "  auto task" << i << " = make_shared<Task" << i << ">(tensor" << i << ", reset);" << endl;
  out.ee << "  " << label_ << "q->add_task(task" << i << ");" << endl << endl;

  return out;
}


OutStream Residual::generate_task(const int ip, const int ic, const vector<string> op, const string scalar, const int i0, bool der, bool diagonal) const {
  stringstream tmp;

  // when there is no gamma under this, we must skip for off-digonal
  string indent = "";

  if (diagonal) {
    tmp << "  shared_ptr<Task" << ic << "> task" << ic << ";" << endl;
    tmp << "  if (diagonal) {" << endl;
    indent += "  ";
  }
  tmp << indent << "  auto tensor" << ic << " = array<shared_ptr<Tensor>," << count_distinct_tensors__(op) << ">{{" << merge__(op, label_) << "}};" << endl;
  tmp << indent << "  " << (diagonal ? "" : "auto ") << "task" << ic
                << " = make_shared<Task" << ic << ">(tensor" << ic << (der || label_=="deci" ? ", cindex" : ", pindex") << (scalar.empty() ? "" : ", this->e0_") << ");" << endl;

  const bool is_gamma = op.front().find("Gamma") != string::npos;
  if (!is_gamma) {
    if (parent_) {
      assert(parent_->parent());
      tmp << indent << "  task" << ip << "->add_dep(task" << ic << ");" << endl;
      tmp << indent << "  task" << ic << "->add_dep(task" << i0 << ");" << endl;
    } else {
      assert(depth() == 0);
      tmp << indent << "  task" << ic << "->add_dep(task" << i0 << ");" << endl;
    }
    tmp << indent << "  " << label_ << "q->add_task(task" << ic << ");" << endl;
  }
  if (diagonal)
    tmp << "  }" << endl;

  if (!is_gamma)
    tmp << endl;

  OutStream out;
  if (!is_gamma) {
    out.ee << tmp.str();
  } else {
    out.gg << tmp.str();
  }
  return out;
}


OutStream Residual::generate_compute_header(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors, const bool no_outside) const {
  vector<string> labels;
  for (auto& i : tensors)
    labels.push_back(i->label());
  const int ninptensors = count_distinct_tensors__(labels);

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  const int nindex = ti.size();
  OutStream out;
  out.tt << "class Task" << ic << " : public Task {" << endl;
  out.tt << "  protected:" << endl;
  // if index is empty give dummy arg
  out.tt << "    std::array<std::shared_ptr<Tensor>," << ninptensors << "> tensor_;" << endl << endl;
  if (need_e0)
    out.tt << "    const double e0_;" << endl;
  out.tt << "    void compute_() override;" << endl;
  out.tt << endl;

  out.dd << "void Task" << ic << "::compute_() {" << endl;
  return out;
}


OutStream Residual::generate_compute_footer(const int ic, const list<shared_ptr<const Index>> ti, const vector<shared_ptr<Tensor>> tensors) const {
  vector<string> labels;
  for (auto& i : tensors)
    labels.push_back(i->label());
  const int ninptensors = count_distinct_tensors__(labels);
  assert(ninptensors > 1);

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  OutStream out;
  out.dd << "}" << endl << endl << endl;

  out.tt << "  public:" << endl;
  out.tt << "    Task" << ic << "(std::array<std::shared_ptr<Tensor>," << ninptensors << "> t" << (need_e0 ? ", const double e" : "") << ") : tensor_(t)"
                       << (need_e0 ? ", e0_(e)" : "") << " { }" << endl;
  out.tt << "};" << endl << endl;
  return out;
}


OutStream Residual::generate_bc(const shared_ptr<BinaryContraction> i) const {
  OutStream out;
  if (depth() != 0) {
    if (i->tensor()->label().find("Gamma") != string::npos)
      out.dd << "  tensor_[1]->init();" << endl;
    if (i->next_target()->label().find("Gamma") != string::npos)
      out.dd << "  tensor_[2]->init();" << endl;

    out.dd << "  auto ta0 = tensor_[0]->tiledarray<" << target_->index().size() << ">();" << endl;
    out.dd << "  auto ta1 = tensor_[1]->tiledarray<" << i->tensor()->index().size() << ">();" << endl;
    if (!same_tensor__(i->tensor()->label(), i->next_target()->label()))
      out.dd << "  auto ta2 = tensor_[2]->tiledarray<" << i->next_target()->index().size() << ">();" << endl;
    out.dd << "  madness::World::get_default().gop.fence();" << endl;
    out.dd << "  " << target_->generate_ta("ta0", true) << " += ";
    // retrieving tensor_
    out.dd << i->tensor()->generate_ta("ta1") << (target_->index().size() == 0 ? ".dot(" : " * ");
    // retrieving subtree_
    string inlabel("ta"); inlabel += (same_tensor__(i->tensor()->label(), i->next_target()->label()) ? "1" : "2");
    out.dd << i->next_target()->generate_ta(inlabel) << (target_->index().size() == 0 ? ").get()" : "") << ";" << endl;
    out.dd << "  madness::World::get_default().gop.fence();" << endl;
    out.dd << "  *tensor_[0] = make_shared<Tensor>(*ta0);" << endl;
  } else {  // now at bc depth 0
    // making residual vector...
    list<shared_ptr<const Index>> proj = i->target_index();
    list<shared_ptr<const Index>> res;
    assert(!(proj.size() & 1));
    for (auto i = proj.begin(); i != proj.end(); ++i, ++i) {
      auto j = i; ++j;
      res.push_back(*j);
      res.push_back(*i);
    }
    auto residual = make_shared<Tensor>(1.0, target_name__(label_), res);
    vector<shared_ptr<Tensor>> op2 = { i->next_target() };
    out << generate_compute_operators(residual, op2, i->dagger());
  }


  return out;
}


shared_ptr<Tensor> Residual::create_tensor(list<shared_ptr<const Index>> dm) const {
 return make_shared<Tensor>(1.0, "r", dm);
}







