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


// TODO one can merge. only difference is the second argument of merge__
OutStream Energy::generate_task(const int ip, const int ic, const vector<string> op, const string scalar, const int, bool, bool) const {
  OutStream out;
  out.ee << "  auto tensor" << ic << " = array<shared_ptr<Tensor>," << count_distinct_tensors__(op) << ">{{" << merge__(op) << "}};" << endl;
  out.ee << "  auto task" << ic << " = make_shared<Task" << ic << ">(tensor" << ic << (scalar.empty() ? "" : ", this->e0_") << ");" << endl;

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


// TODO one can merge. only difference is "Acc"
OutStream Energy::generate_compute_header(const int ic, const vector<shared_ptr<Tensor>> tensors) const {
  vector<string> labels;
  for (auto& i : tensors)
    labels.push_back(i->label());
  const int ninptensors = count_distinct_tensors__(labels);
  assert(ninptensors > 1);

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  OutStream out;
  out.tt << "class Task" << ic << " : public AccTask {" << endl;
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


OutStream Energy::generate_bc(const shared_ptr<BinaryContraction> i) const {
  OutStream out;

  if (depth() != 0) {
    if (i->tensor()->label().find("Gamma") != string::npos)
      out.dd << "  tensor_[1]->init();" << endl;
    if (i->next_target()->label().find("Gamma") != string::npos)
      out.dd << "  tensor_[2]->init();" << endl;

    out.dd << "  auto ta0 = tensor_[0]->tiledarray<" << target_->index().size() << ">(true);" << endl;
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
    if (depth() > 1) {
      out.dd << "  *tensor_[0] = *make_shared<Tensor>(*ta0);" << endl;
    } else {
      out.dd << "  target_ += (*ta0)(\"\");" << endl;
    }
  } else {
    // depth should not equal 0 in energy tree
    throw logic_error("shouldn't happen in Energy::generate_bc");
  }
  return out;
}

