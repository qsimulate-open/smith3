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
#include "tree.h"

using namespace std;
using namespace smith;


// local functions... (not a good practice...) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
static string label__(const string& lbl) {
  string label = lbl;
  if (label == "proj") label = "r";
  size_t found = label.find("dagger");
  if (found != string::npos) {
    string tmp(label.begin(), label.begin() + found);
    label = tmp;
  }
  return label;
}
// local functions... (not a good practice...) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


OutStream Tree::create_target(const int i) const {
  OutStream out;
  if (!is_energy_tree()) {
    out.tt << "class Task" << i << " : public Task {" << endl;
    out.tt << "  protected:" << endl;
    out.tt << "    std::shared_ptr<TATensor<" << DataType << ",4>> tensor_;" << endl;
    out.tt << "    const bool reset_;" << endl;
    out.tt << "" << endl;
    out.tt << "    void compute_() {" << endl;
    out.tt << "      if (reset_) tensor_[0]->zero();" << endl;
    out.tt << "    }" << endl;
    out.tt << "" << endl;
    out.tt << "  public:" << endl;
    out.tt << "    Task" << i << "(std::shared_ptr<TATensor<" << DataType << ",4>> t, const bool reset) : tensor_(t), reset_(reset) { }" << endl;
    out.tt << "};" << endl << endl;

    out.ee << "  auto " << label_ << "q = make_shared<Queue>();" << endl;
    out.ee << "  auto task" << i << " = make_shared<Task" << i << ">(" << target_name__(label_) << ", reset);" << endl;
    out.ee << "  " << label_ << "q->add_task(task" << i << ");" << endl << endl;
  }
  return out;
}


OutStream Tree::generate_task(const int ip, const int ic, const vector<string> op, const string scalar, const int i0, bool der, bool diagonal) const {
  stringstream tmp;

  // when there is no gamma under this, we must skip for off-digonal
  string indent = "";

  const bool is_gamma = op.front().find("Gamma") != string::npos;

  if (diagonal) {
    tmp << "  shared_ptr<Task" << ic << "> task" << ic << ";" << endl;
    tmp << "  if (diagonal) {" << endl;
    indent += "  ";
  }
  tmp << indent << "  " << (diagonal ? "" : "auto ") << "task" << ic
                << " = make_shared<Task" << ic << ">(" << merge__(op, label_) << (is_gamma ? (der ? ", cindex" : ", pindex") : "")
                << (scalar.empty() ? "" : ", this->e0_") << ");" << endl;

  if (!is_gamma) {
    if (parent_) {
      assert(parent_->parent());
      if (ip != ic)
        tmp << indent << "  task" << ip << "->add_dep(task" << ic << ");" << endl;
    } else
      assert(depth() == 0);

    if (!is_energy_tree())
      tmp << indent << "  task" << ic << "->add_dep(task" << i0 << ");" << endl;
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


OutStream Tree::generate_compute_header(const int ic, const vector<shared_ptr<Tensor>> tensors) const {
  vector<string> labels;
  for (auto& i : tensors)
    labels.push_back(i->label());
  const int ninptensors = count_distinct_tensors__(labels);

  bool need_e0 = false;
  for (auto& s : tensors)
    if (!s->scalar().empty()) need_e0 = true;

  OutStream out;
  out.tt << "class Task" << ic << " : public " << (is_energy_tree() ? "Acc" : "") << "Task {" << endl;
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


OutStream Tree::generate_bc(const shared_ptr<BinaryContraction> i) const {
  OutStream out;
  if (depth() != 0) {
    out.dd << "  if (!ta0_->initialized())" << endl;
    out.dd << "    ta0_->fill_local(0.0);" << endl;
    if (i->tensor()->label().find("Gamma") != string::npos)
      out.dd << "  ta1_->init();" << endl;
    if (i->next_target()->label().find("Gamma") != string::npos)
      out.dd << "  ta2_->init();" << endl;

    out.dd << "  madness::World::get_default().gop.fence();" << endl;
    out.dd << "  " << target_->generate_ta("ta0_", true) << " += ";
    // retrieving tensor_
    out.dd << i->tensor()->generate_ta("ta1_") << (target_->index().size() == 0 ? ".dot(" : " * ");
    // retrieving subtree_
    string inlabel("ta"); inlabel += (same_tensor__(i->tensor()->label(), i->next_target()->label()) ? "1_" : "2_");
    out.dd << i->next_target()->generate_ta(inlabel) << (target_->index().size() == 0 ? ").get()" : "") << ";" << endl;
    out.dd << "  madness::World::get_default().gop.fence();" << endl;
    if (is_energy_tree() && depth() == 1)
      out.dd << "  target_ += (*ta0)(\"\");" << endl;
  } else {  // now at bc depth 0
    assert(!is_energy_tree());
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


OutStream Tree::generate_compute_operators(shared_ptr<Tensor> target, const vector<shared_ptr<Tensor>> op, const bool dagger) const {
  OutStream out;
  // needed in case the tensor labels are repeated..
  vector<shared_ptr<Tensor>> uniq_tensors;
  vector<string> tensor_labels;
  // map redundant tensor label to unique tensor label
  map<string, int> op_tensor_lab;

  int uniq_cnt = 1;
  for (auto& i : op) {
    string label = label__(i->label());
    // if we already have the label in the list
    if (find(tensor_labels.begin(), tensor_labels.end(), label) != tensor_labels.end()) {
      continue;
    }
    // otherwise not yet found, so add to list
    tensor_labels.push_back(label);
    uniq_tensors.push_back(i);
    op_tensor_lab[label] = uniq_cnt;
    ++uniq_cnt;
  }

  if (depth() != 0) {
    out.dd << "  if (!ta0_->initialized())" << endl;
    out.dd << "    ta0_->fill_local(0.0);" << endl;
  }
  out.dd << "  madness::World::get_default().gop.fence();" << endl;
  out.dd << "  " << target->generate_ta("ta0_", true) << " += ";

  // add the source data to the target
  int j = 0;
  for (auto s = op.begin(); s != op.end(); ++s, ++j) {
    // uses map to give label number consistent with operator, needed in case label is repeated (eg ccaa)
    string label = label__((*s)->label());
    stringstream instr; instr << "ta" << op_tensor_lab[label] << "_";

    out.dd << (s == op.begin() ? "" : "\n     + ") << (*s)->generate_ta(instr.str());

    // if this is at the top-level and needs to be daggered:
    if (depth() == 0 && dagger) {
      shared_ptr<Tensor> top = make_shared<Tensor>(**s);
      // swap operators so that tensor is daggered
      if (top->index().size() != 4) {
         throw logic_error("Daggered object is only supported for 4-index tensors");
      } else {
        list<shared_ptr<const Index>> di = target->index();
        auto k0 = di.begin(); auto k1 = k0; ++k1; auto k2 = k1; ++k2; auto k3 = k2; ++k3;
        using Iter = list<shared_ptr<const Index>>::iterator;
        const list<pair<Iter, Iter>> map {{k0, k2}, {k2, k0}, {k1, k3}, {k3, k1}};
        list<shared_ptr<const Index>> tmp;
        for (auto& k : top->index()) {
          for (auto l = map.begin(); l != map.end(); ++l) {
            if (k->identical(*l->first)) {
              tmp.push_back(*l->second);
              break;
            }
            auto ll = l;
            if (++ll == map.end()) throw logic_error("should not happen: dagger");
          }
        }
        top->index() = tmp;
        out.dd << " + " << top->generate_ta(instr.str());
      }
    }
  }
  out.dd << ";" << endl;
  out.dd << "  madness::World::get_default().gop.fence();" << endl;
  return out;
}


OutStream Tree::generate_task(const int ic, const vector<shared_ptr<Tensor>> op, const list<shared_ptr<Tensor>> g, const int iz, const bool diagonal) const {
  OutStream out;

  vector<string> ops;
  for (auto& i : op) ops.push_back(i->label());
  int ip = -1;
  if (parent_) ip = parent_->parent()->num();

  string scalar;
  for (auto& i : op) {
    if (!i->scalar().empty()) {
      if (i->scalar() != "e0") throw logic_error("unknown scalar. Tree::generate_task");
      if (!scalar.empty() && i->scalar() != scalar) throw logic_error("multiple scalars. Tree::generate_task");
      scalar = i->scalar();
    }
  }

  // if gamma, we need to add dependency.
  // this one is virtual, ie tree specific
  out << generate_task(ip, ic, ops, scalar, iz, /*der*/false, /*diagonal*/diagonal);

// TODO at this moment all gammas are recomputed.
#if 0
  for (auto& i : op) {
    if (i->label().find("Gamma") != string::npos)
      out.ee << "  task" << ic << "->" << add_depend(i, g) << endl;
  }
  out.ee << endl;
#endif

  return out;
}


string Tree::add_depend(const shared_ptr<const Tensor> o, const list<shared_ptr<Tensor>> gamma) const {
  stringstream out;
  if (!parent_) {
    assert(!gamma.empty());
    for (auto& i : gamma) {
      if (o->label() == i->label()) out << "add_dep(task" << i->num() << ");";
    }
  } else {
    out << parent_->parent()->add_depend(o, gamma);
  }
  return out.str();
}


tuple<OutStream, int, int, vector<shared_ptr<Tensor>>>
  BinaryContraction::generate_task_list(int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  OutStream out, tmp;
  string depends, tasks, specials;

  for (auto& i : subtree_) {
    tie(tmp, tcnt, t0, itensors) = i->generate_task_list(tcnt, t0, gamma, itensors);
    out << tmp;
  }
  return make_tuple(out, tcnt, t0, itensors);
}


tuple<OutStream, int, int, vector<shared_ptr<Tensor>>>
  Tree::generate_task_list(int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  // here ss is dependencies, tt is tasks
  OutStream out, tmp;

  if (depth() == 0) { //////////////////// zero depth /////////////////////////////
    if (root_targets()) {
      // process tree with target indices eg, ci derivative, density matrix
      num_ = tcnt;
      // save density task zero
      t0 = tcnt;

      // virtual
      out << create_target(tcnt);
      ++tcnt;

      for (auto& j : bc_) {
        // if at top bc, add a task to for top level contraction (proj)
        if (depth() == 0 ) {
          vector<shared_ptr<Tensor>> source_tensors = j->tensors_vec();
          const bool diagonal = j->diagonal_only();
          num_ = tcnt;
          for (auto& s : source_tensors) {
            // if it contains a new intermediate tensor, dump a constructor
            if (find(itensors.begin(), itensors.end(), s) == itensors.end() && s->label().find("I") != string::npos) {
              itensors.push_back(s);
              out.ee << s->constructor_str(diagonal) << endl;
            }
          }
          out << generate_task(num_, source_tensors, gamma, t0, diagonal);

          list<shared_ptr<const Index>> proj = j->target_index();
          // write out headers
          out << generate_compute_header(num_, source_tensors);

          list<shared_ptr<const Index>> dm;
          if (proj.size() > 1) {
            for (auto i = proj.begin(); i != proj.end(); ++i, ++i) {
              auto k = i; ++k;
              dm.push_back(*k);
              dm.push_back(*i);
            }
          } else if (proj.size() == 1) {
            for (auto& i : proj) dm.push_back(i);
          } else {
            throw logic_error("Tree::generate_task_list, should not have empty target index here.");
          }

          // virtual
          shared_ptr<Tensor> proj_tensor = create_tensor(dm);

          vector<shared_ptr<Tensor>> op2 = { j->next_target() };
          out << generate_compute_operators(proj_tensor, op2, j->dagger());
          out << generate_compute_footer(num_, source_tensors);
          // increment task counter
          ++tcnt;
        }

        tie(tmp, tcnt, t0, itensors) = j->generate_task_list(tcnt, t0, gamma, itensors);
        out << tmp;

      }

    } else {  // trees without root target indices
      out.ee << "  auto " << label() << "q = make_shared<Queue>();" << endl;
      num_ = tcnt;
      for (auto& j : bc_) {
        tie(tmp, tcnt, t0, itensors) = j->generate_task_list(tcnt, t0, gamma, itensors);
        out << tmp;
      }
    }
  } else { //////////////////// non-zero depth /////////////////////////////
    tie(tmp, tcnt, t0, itensors) = generate_steps(tcnt, t0, gamma, itensors);
    out << tmp;
  }

  return make_tuple(out, tcnt, t0, itensors);
}


tuple<OutStream, int, int, vector<shared_ptr<Tensor>>>
    Tree::generate_steps(int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  OutStream out, tmp;
  /////////////////////////////////////////////////////////////////
  // if op_ is not empty, we add a task that adds up op_.
  /////////////////////////////////////////////////////////////////
  if (!op_.empty()) {

    // step through operators and if they are new, construct them.
    if (find(itensors.begin(), itensors.end(), target_) == itensors.end()) {
      itensors.push_back(target_);
      out.ee << target_->constructor_str(diagonal_only()) << endl;
    }

    vector<shared_ptr<Tensor>> op = {target_};
    op.insert(op.end(), op_.begin(), op_.end());
    // check if op_ has gamma tensor
    bool diagonal = nogamma_upstream();
    for (auto& i : op_) {
      string label = i->label();
      diagonal &= label.find("Gamma") == string::npos;
    }
    out << generate_task(tcnt, op, gamma, t0, diagonal);

    // make sure no duplicates in tensor list for compute header & footer
    vector<shared_ptr<Tensor>> uniq_tensors;
    vector<string> tensor_labels;
    for (auto& i : op) {
      string label = label__(i->label());
      if (find(tensor_labels.begin(), tensor_labels.end(), label) != tensor_labels.end()) continue;
      tensor_labels.push_back(label);
      uniq_tensors.push_back(i);
    }

    out << generate_compute_header(tcnt, uniq_tensors);
    out << generate_compute_operators(target_, op_);
    out << generate_compute_footer(tcnt, uniq_tensors);

    ++tcnt;
  }

  /////////////////////////////////////////////////////////////////
  // step through BinaryContraction
  /////////////////////////////////////////////////////////////////
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    vector<shared_ptr<Tensor>> source_tensors = (*i)->tensors_vec();

    const bool diagonal = (*i)->diagonal_only();
    for (auto& s : source_tensors) {
      // if it contains a new intermediate tensor, dump a constructor
      if (find(itensors.begin(), itensors.end(), s) == itensors.end() && s->label().find("I") != string::npos) {
        itensors.push_back(s);
        out.ee << s->constructor_str(diagonal) << endl;
      }
    }
    // saving a counter to a protected member for dependency checks
    num_ = tcnt;
    out << generate_task(num_, source_tensors, gamma, t0, diagonal);

    // write out headers
    out << generate_compute_header(num_, source_tensors);

    // use virtual function to generate a task for this binary contraction
    out << generate_bc(*i);
    out << generate_compute_footer(num_, source_tensors);

    // increment tcnt before going to subtrees
    ++tcnt;
    // triggers a recursive call
    tie(tmp, tcnt, t0, itensors) = (*i)->generate_task_list(tcnt, t0, gamma, itensors);
    out << tmp;
  } // end bc

  return make_tuple(out, tcnt, t0, itensors);
}


OutStream Tree::generate_compute_footer(const int ic, const vector<shared_ptr<Tensor>> tensors) const {
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

