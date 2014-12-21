//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: tree.cc
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


#include "tree.h"
#include "energy.h"
#include "residual.h"
#include "density.h"
#include "density1.h"
#include "density2.h"
#include "correction.h"
#include "dedci.h"
#include "constants.h"
#include <algorithm>
#include <stdexcept>
#include <map>

using namespace std;
using namespace smith;


Tree::Tree(shared_ptr<Equation> eq, string lab) : parent_(NULL), tree_name_(eq->name()), num_(-1), label_(lab), root_targets_(eq->targets()) {
  // First make ListTensor for all the diagrams
  list<shared_ptr<Diagram>> d = eq->diagram();

  const bool rt_targets = eq->targets();

  for (auto i = d.begin(); i != d.end(); ++i) {
    shared_ptr<ListTensor> tmp = make_shared<ListTensor>(*i);
    // All internal tensor should be included in the active part
    tmp->absorb_all_internal();

    // rearrange brakets and reindex associated tensors, ok if not complex
    tmp->absorb_ket();

    shared_ptr<Tensor> first = tmp->front();
    shared_ptr<ListTensor> rest = tmp->rest();

    // convert to tree and then bc
    shared_ptr<Tree> tr;
    if (label_ == "residual") {
      tr = make_shared<Residual>(rest, lab, rt_targets);
    } else if (label_ == "energy") {
      tr = make_shared<Energy>(rest, lab, rt_targets);
    } else if (label_ == "dedci") {
      tr = make_shared<Dedci>(rest, lab, rt_targets);
    } else if (label_ == "correction") {
      tr = make_shared<Correction>(rest, lab, rt_targets);
    } else if (label_ == "density") {
      tr = make_shared<Density>(rest, lab, rt_targets);
    } else if (label_ == "density1") {
      tr = make_shared<Density1>(rest, lab, rt_targets);
    } else if (label_ == "density2") {
      tr = make_shared<Density2>(rest, lab, rt_targets);
    } else {
      throw logic_error("Error Tree::Tree, code generation for this tree type not implemented");
    }
    list<shared_ptr<Tree>> lt; lt.push_back(tr);
    shared_ptr<BinaryContraction> b = make_shared<BinaryContraction>(lt, first, (*i)->target_index());
    bc_.push_back(b);
  }

  factorize();
  set_parent_sub();
  set_target_rec();
}


Tree::Tree(const shared_ptr<ListTensor> l, string lab, const bool t) : num_(-1), label_(lab), root_targets_(t) {
  target_ = l->target();
  dagger_ = l->dagger();
  if (l->length() > 1) {
    shared_ptr<BinaryContraction> bc = make_shared<BinaryContraction>(target_, l, lab, root_targets_);
    bc_.push_back(bc);
  } else {
    shared_ptr<Tensor> t = l->front();
    t->set_factor(l->fac());
    t->set_scalar(l->scalar());
    op_.push_back(t);
  }
}


BinaryContraction::BinaryContraction(shared_ptr<Tensor> o, shared_ptr<ListTensor> l, string lab, const bool rt) : label_(lab) {
  // o is a target tensor NOT included in l
  target_ = o;
  tensor_ = l->front();
  shared_ptr<ListTensor> rest = l->rest();

  shared_ptr<Tree> tr;
  if (label_ == "residual") {
    tr = make_shared<Residual>(rest, lab, rt);
  } else if (label_ == "energy") {
    tr = make_shared<Energy>(rest, lab, rt);
  } else if (label_ == "dedci") {
    tr = make_shared<Dedci>(rest, lab, rt);
  } else if (label_ == "correction") {
    tr = make_shared<Correction>(rest, lab, rt);
  } else if (label_ == "density") {
    tr = make_shared<Density>(rest, lab, rt);
  } else if (label_ == "density1") {
    tr = make_shared<Density1>(rest, lab, rt);
  } else if (label_ == "density2") {
    tr = make_shared<Density2>(rest, lab, rt);
  } else {
    throw logic_error("Error BinaryContraction::BinaryContraction, code generation for this tree type not implemented");
  }
  subtree_.push_back(tr);

}


void Tree::factorize() {
  list<list<shared_ptr<BinaryContraction>>::iterator> done;
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    auto j = i; ++j;
    for ( ; j != bc_.end(); ++j) {
      // bc is factorized according to excitation target indices
      if (*(*i)->tensor() == *(*j)->tensor() && (*i)->dagger() == (*j)->dagger() && (*i)->target_index_str() == (*j)->target_index_str()) {
        done.push_back(j);
        (*i)->subtree().insert((*i)->subtree().end(), (*j)->subtree().begin(), (*j)->subtree().end());
      }
    }
  }
  for (auto& i : done) bc_.erase(i);
  for (auto& i : bc_) i->factorize();
}


void BinaryContraction::factorize() {
  list<list<shared_ptr<Tree>>::iterator> done;
  // and then merge subtrees
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    auto j = i; ++j;
    for ( ; j != subtree_.end(); ++j) {
      if (find(done.begin(), done.end(), j) != done.end()) continue;
      if ((*i)->merge(*j)) done.push_back(j);
    }
  }
  for (auto i = done.begin(); i != done.end(); ++i) subtree_.erase(*i);
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) (*i)->factorize();
}


bool Tree::merge(shared_ptr<Tree> o) {
  bool out = false;
  if (o->bc_.size() > 0) {
    shared_ptr<Tensor> a = bc_.front()->tensor();
    shared_ptr<Tensor> b = o->bc_.front()->tensor();
    if (*a == *b) {
      out = true;
      bc_.insert(bc_.end(), o->bc_.begin(), o->bc_.end());
      for (auto i = bc_.begin(); i != bc_.end(); ++i) (*i)->set_target(target_);
    }
  } else if (o->op_.size() > 0) {
    out = true;
    op_.insert(op_.end(), o->op_.begin(), o->op_.end());
  }
  return out;
}


void BinaryContraction::set_target_rec() {
  if (!subtree_.empty()) {
    auto i = subtree_.begin();
    shared_ptr<Tensor> first = (*i)->target();
    for (++i ; i != subtree_.end(); ++i) {
      (*i)->set_target(first);
    }
    for (i = subtree_.begin() ; i != subtree_.end(); ++i)
      (*i)->set_target_rec();
  }
}


void Tree::set_target_rec() {
  if (!bc_.empty()) {
    for (auto& i : bc_) i->set_target_rec();
  }
}


void BinaryContraction::set_parent_subtree() {
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    (*i)->set_parent(this);
    (*i)->set_parent_sub();
  }
}


void Tree::set_parent_sub() {
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    (*i)->set_parent(this);
    (*i)->set_parent_subtree();
  }
}


int BinaryContraction::depth() const { return parent_->depth(); }

int Tree::depth() const { return parent_ ? parent_->parent()->depth()+1 : 0; }


bool BinaryContraction::dagger() const {
  return subtree_.front()->dagger();
}


shared_ptr<Tensor> BinaryContraction::next_target() {
  return subtree().front()->target();
}


void Tree::sort_gamma(std::list<std::shared_ptr<Tensor>> o) {
  gamma_ = o;
  list<shared_ptr<Tensor>> g = gather_gamma();
  for (auto& i : g) find_gamma(i);
}


void Tree::find_gamma(shared_ptr<Tensor> o) {
  bool found = false;
  for (auto& i : gamma_) {
    if ((*i) == (*o)) {
      found = true;
      o->set_alias(i);
      break;
    }
  }
  if (!found) gamma_.push_back(o);
}


// return gamma below this
list<shared_ptr<Tensor>> Tree::gather_gamma() const {
  list<shared_ptr<Tensor>> out;
  for (auto& i : bc_) {
    for (auto& j : i->subtree()) {
      // recursive call
      list<shared_ptr<Tensor>> tmp = j->gather_gamma();
      out.insert(out.end(), tmp.begin(), tmp.end());
    }
  }
  for (auto& i : op_)
    if (i->label().find("Gamma") != string::npos) out.push_back(i);
  return out;
}


string BinaryContraction::target_index_str() const {
  stringstream zz;
  if (!target_index_.empty()) {
    zz << "(";
    for (auto i = target_index_.begin(); i != target_index_.end(); ++i) {
      if (i != target_index_.begin()) zz << ", ";
      zz << (*i)->str(false);
    }
    zz << ")";
  }
  return zz.str();
}


void BinaryContraction::print() const {
  string indent = "";
  for (int i = 0; i != depth(); ++i) indent += "  ";
  cout << indent << (target_ ? (target_->str() + " = ") : "") << tensor_->str() << (depth() == 0 ? target_index_str() : "") << " * " << subtree_.front()->target()->str() << endl;
  for (auto& i : subtree_) i->print();
}


void Tree::print() const {
  string indent = "";
  for (int i = 0; i != depth(); ++i) indent += "  ";
  if (target_) {
    for (auto i = op_.begin(); i != op_.end(); ++i)
      cout << indent << target_->str() << " += " << (*i)->str() << (dagger_ ? " *" : "") << endl;
  }
  for (auto i = bc_.begin(); i != bc_.end(); ++i)
    (*i)->print();
}


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



list<shared_ptr<const Index>> BinaryContraction::target_indices() {
  // returns a list of target indices
  return target_->index();
}


list<shared_ptr<const Index>> BinaryContraction::loop_indices() {
  // returns a list of inner loop indices.
  list<shared_ptr<const Index>> out;
  list<shared_ptr<const Index>> ti = target_->index();
  for (auto iter = tensor_->index().begin(); iter != tensor_->index().end(); ++iter) {
    bool found = false;
    for (auto i = ti.begin(); i != ti.end(); ++i) {
      if ((*i)->identical(*iter)) {
        found = true;
        break;
      }
    }
    if (!found) {
      out.push_back(*iter);
    }
  }
  return out;
}


string Tree::generate_compute_operators(const string indent, const shared_ptr<Tensor> target, const vector<shared_ptr<Tensor>> op,
                                        const bool dagger) const {
  stringstream tt;

  vector<string> close;
  string cindent = indent + "    ";
  // note that I am using reverse_iterator
  //tt << target->generate_loop(cindent, close);
  // get data
  tt << target->generate_get_block(cindent, "o", "out()", true);

  // needed in case the tensor labels are repeated..
  vector<shared_ptr<Tensor>> uniq_tensors;
  vector<string> tensor_labels;
  // map redundant tensor label to unique tensor label
  map<string, int> op_tensor_lab;

  int uniq_cnt = 0;
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

  // add the source data to the target
  int j = 0;
  for (auto s = op.begin(); s != op.end(); ++s, ++j) {
    stringstream uu; uu << "i" << j;
    tt << cindent << "{" << endl;

    // uses map to give label number consistent with operator, needed in case label is repeated (eg ccaa)
    string label = label__((*s)->label());
    stringstream instr; instr << "in(" << op_tensor_lab[label] << ")";

    tt << (*s)->generate_get_block(cindent+"  ", uu.str(), instr.str());
    list<shared_ptr<const Index>> di = target->index();
    tt << (*s)->generate_sort_indices(cindent+"  ", uu.str(), instr.str(), di, true);
    tt << cindent << "}" << endl;

    // if this is at the top-level and needs to be daggered:
    if (depth() == 0 && dagger) {
      shared_ptr<Tensor> top = make_shared<Tensor>(**s);
      // swap operators so that tensor is daggered
      if (top->index().size() != 4) {
         throw logic_error("Daggered object is only supported for 4-index tensors");
      } else {
        auto k0 = di.begin(); auto k1 = k0; ++k1; auto k2 = k1; ++k2; auto k3 = k2; ++k3;
        list<pair<list<shared_ptr<const Index>>::iterator, list<shared_ptr<const Index>>::iterator>> map;
        map.push_back(make_pair(k0, k2));
        map.push_back(make_pair(k2, k0));
        map.push_back(make_pair(k1, k3));
        map.push_back(make_pair(k3, k1));
        list<shared_ptr<const Index>> tmp;
        for (auto& k : top->index()) {
          for (auto l = map.begin(); l != map.end(); ++l) {
            if (k->identical(*l->first)) {
              tmp.push_back(*l->second);
              break;
            }
            auto ll = l;
            if (++ll == map.end()) throw logic_error("should not happen: dagger stuffs");
          }
        }
        top->index() = tmp;
        tt << cindent << "{" << endl;
        tt << top->generate_get_block(cindent+"  ", uu.str(), instr.str());
        list<shared_ptr<const Index>> di = target->index();
        tt << top->generate_sort_indices(cindent+"  ", uu.str(), instr.str(), di, true);
        tt << cindent << "}" << endl;
      }
    }
  }

  // put data
  {
    string label = target->label();
    list<shared_ptr<const Index>> ti = target->index();
    tt << cindent << "out()->put_block(odata";
    for (auto i = ti.rbegin(); i != ti.rend(); ++i)
      tt << ", " << (*i)->str_gen();
    tt << ");" << endl;
  }
  for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
        tt << *iter << endl;
  return tt.str();
}


string Tree::generate_task(const string indent, const int ic, const vector<shared_ptr<Tensor>> op, const list<shared_ptr<Tensor>> g, const int iz) const {
  stringstream ss;

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
  ss << generate_task(indent, ip, ic, ops, scalar, iz);
  for (auto& i : op) {
    if (i->label().find("Gamma") != string::npos)
      ss << indent << "task" << ic << "->" << add_depend(i, g) << endl;
  }
  ss << endl;

  return ss.str();
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


tuple<string, string, string, int, int, vector<shared_ptr<Tensor>>>
  BinaryContraction::generate_task_list(int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  stringstream ss, tt, cc;
  string depends, tasks, specials;

  for (auto& i : subtree_) {
    tie(depends, tasks, specials, tcnt, t0, itensors) = i->generate_task_list(tcnt, t0, gamma, itensors);
    ss << depends;
    tt << tasks;
    cc << specials;
  }
  return make_tuple(ss.str(), tt.str(), cc.str(), tcnt, t0, itensors);
}


tuple<string, string, string, int, int, vector<shared_ptr<Tensor>>>
  Tree::generate_task_list(int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  // here ss is dependencies, tt is tasks
  stringstream ss, tt, cc;
  string depends, tasks, specials;
  string indent = "      ";

  if (depth() == 0) { //////////////////// zero depth /////////////////////////////
    if (root_targets()) {
      // process tree with target indices eg, ci derivative, density matrix (residual zero level task is created above (forest::generate_header))
      if (label() != "residual") {
        num_ = tcnt;
        // save density task zero
        t0 = tcnt;

        // virtual
        auto rtmp = create_target(indent, tcnt);
        ss << get<0>(rtmp);
        tt << get<1>(rtmp);
        cc << get<2>(rtmp);
        ++tcnt;

        for (auto& j : bc_) {
          // if at top bc, add a task to for top level contraction (proj)
          if (depth() == 0 ) {
            vector<shared_ptr<Tensor>> source_tensors = j->tensors_vec();
            num_ = tcnt;
            for (auto& s : source_tensors) {
              // if it contains a new intermediate tensor, dump a constructor
              if (find(itensors.begin(), itensors.end(), s) == itensors.end() && s->label().find("I") != string::npos) {
                itensors.push_back(s);
                ss << s->constructor_str(indent) << endl;
              }
            }
            ss << generate_task(indent, num_, source_tensors, gamma, t0);

            list<shared_ptr<const Index>> proj = j->target_index();
            // write out headers
            {
              list<shared_ptr<const Index>> ti = depth() != 0 ? j->target_indices() : proj;
              // if outer loop is empty, send inner loop indices to header
              if (ti.size() == 0) {
                assert(depth() != 0);
                list<shared_ptr<const Index>> di = j->loop_indices();
                tt << generate_compute_header(num_, di, source_tensors, true);
              } else {
                tt << generate_compute_header(num_, ti, source_tensors);
              }
            }

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
            tt << generate_compute_operators(indent, proj_tensor, op2, j->dagger());


            {
              // send outer loop indices if outer loop indices exist, otherwise send inner indices
              list<shared_ptr<const Index>> ti = depth() != 0 ? j->target_indices() : proj;
              if (depth() == 0)
                if (ti.size() > 1) {
                  for (auto m = ti.begin(), n = ++ti.begin(); m != ti.end(); ++m, ++m, ++n, ++n)
                    swap(*m, *n);
                }
              if (ti.size() == 0) {
                assert(depth() != 0);
                // sending inner indices
                list<shared_ptr<const Index>> di = j->loop_indices();
                di.reverse();
                string stringtt, stringcc;
                tie(stringtt, stringcc) = generate_compute_footer(num_, di, source_tensors);
                tt << stringtt;
                cc << stringcc;
              } else {
                // sending outer indices
                string stringtt, stringcc;
                tie(stringtt, stringcc) = generate_compute_footer(num_, ti, source_tensors);
                tt << stringtt;
                cc << stringcc;
              }
            }
            // increment task counter
            ++tcnt;
          }

          tie(depends, tasks, specials, tcnt, t0, itensors) = j->generate_task_list(tcnt, t0, gamma, itensors);
          ss << depends;
          tt << tasks;
          cc << specials;

        }
      } else { // residual

        // step through op and bc. Careful, triggers recursive call
        tie(depends, tasks, specials, tcnt, t0, itensors) = generate_steps(indent, tcnt, t0, gamma, itensors);
        ss << depends;
        tt << tasks;
        cc << specials;
      }

    } else {  // trees without root target indices
      ss << "      auto " << label() << "_ = std::make_shared<Queue>();" << endl;
      num_ = tcnt;
      for (auto& j : bc_) {
        tie(depends, tasks, specials, tcnt, t0, itensors) = j->generate_task_list(tcnt, t0, gamma, itensors);
        ss << depends;
        tt << tasks;
        cc << specials;
      }
    }
  } else { //////////////////// non-zero depth /////////////////////////////
    tie(depends, tasks, specials, tcnt, t0, itensors) = generate_steps(indent, tcnt, t0, gamma, itensors);
    ss << depends;
    tt << tasks;
    cc << specials;
  }

  return make_tuple(ss.str(), tt.str(), cc.str(), tcnt, t0, itensors);
}


tuple<string, string, string, int, int, vector<shared_ptr<Tensor>>>
    Tree::generate_steps(string indent, int tcnt, int t0, const list<shared_ptr<Tensor>> gamma, vector<shared_ptr<Tensor>> itensors) const {
  stringstream ss, tt, cc;
  string depends, tasks, specials;
  /////////////////////////////////////////////////////////////////
  // if op_ is not empty, we add a task that adds up op_.
  /////////////////////////////////////////////////////////////////
  if (!op_.empty()) {

    // step through operators and if they are new, construct them.
    if (find(itensors.begin(), itensors.end(), target_) == itensors.end()) {
      itensors.push_back(target_);
      ss << target_->constructor_str(indent) << endl;
    }

    vector<shared_ptr<Tensor>> op = {target_};
    op.insert(op.end(), op_.begin(), op_.end());
    ss << generate_task(indent, tcnt, op, gamma, t0);

    list<shared_ptr<const Index>> ti = target_->index();

    // make sure no duplicates in tensor list for compute header & footer
    vector<shared_ptr<Tensor>> uniq_tensors;
    vector<string> tensor_labels;
    for (auto& i : op) {
      string label = label__(i->label());
      if (find(tensor_labels.begin(), tensor_labels.end(), label) != tensor_labels.end()) continue;
      tensor_labels.push_back(label);
      uniq_tensors.push_back(i);
    }

    tt << generate_compute_header(tcnt, ti, uniq_tensors);
    tt << generate_compute_operators(indent, target_, op_);

    string stringtt, stringcc;
    tie(stringtt, stringcc) = generate_compute_footer(tcnt, ti, uniq_tensors);
    tt << stringtt;
    cc << stringcc;

    ++tcnt;
  }

  /////////////////////////////////////////////////////////////////
  // step through BinaryContraction
  /////////////////////////////////////////////////////////////////
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    vector<shared_ptr<Tensor>> source_tensors = (*i)->tensors_vec();

    for (auto& s : source_tensors) {
      // if it contains a new intermediate tensor, dump a constructor
      if (find(itensors.begin(), itensors.end(), s) == itensors.end() && s->label().find("I") != string::npos) {
        itensors.push_back(s);
        ss << s->constructor_str(indent) << endl;
      }
    }
    // saving a counter to a protected member for dependency checks
    num_ = tcnt;
    ss << generate_task(indent, num_, source_tensors, gamma, t0);

    // write out headers
    {
      list<shared_ptr<const Index>> ti = depth() != 0 ? (*i)->target_indices() : (*i)->target_index();
      // if outer loop is empty, send inner loop indices to header
      if (ti.size() == 0) {
        assert(depth() != 0);
        list<shared_ptr<const Index>> di = (*i)->loop_indices();
        tt << generate_compute_header(num_, di, source_tensors, true);
      } else {
        tt << generate_compute_header(num_, ti, source_tensors);
      }
    }

    // use virtual function to generate a task for this binary contraction
    pair<string, string> btmp = generate_bc(indent, (*i));
    ss << btmp.first;
    tt << btmp.second;

    {
      // send outer loop indices if outer loop indices exist, otherwise send inner indices
      list<shared_ptr<const Index>> ti = depth() != 0 ? (*i)->target_indices() : (*i)->target_index();
      if (depth() == 0)
        for (auto i = ti.begin(), j = ++ti.begin(); i != ti.end(); ++i, ++i, ++j, ++j)
          swap(*i, *j);
      if (ti.size() == 0) {
        assert(depth() != 0);
        // sending inner indices
        list<shared_ptr<const Index>> di = (*i)->loop_indices();
        di.reverse();
        string stringtt, stringcc;
        tie(stringtt, stringcc) = generate_compute_footer(num_, di, source_tensors);
        tt << stringtt;
        cc << stringcc;
      } else {
        // sending outer indices
        string stringtt, stringcc;
        tie(stringtt,stringcc) = generate_compute_footer(num_, ti, source_tensors);
        tt << stringtt;
        cc << stringcc;
      }
    }

    ///////////////////////////////////////////////////////////////////////

    // increment tcnt before going to subtrees
    ++tcnt;
    // triggers a recursive call
    tie(depends, tasks, specials, tcnt, t0, itensors) = (*i)->generate_task_list(tcnt, t0, gamma, itensors);
    ss << depends;
    tt << tasks;
    cc << specials;
  } // end bc

  return make_tuple(ss.str(), tt.str(), cc.str(), tcnt, t0, itensors);
}


vector<shared_ptr<Tensor>> BinaryContraction::tensors_vec() {
  vector<shared_ptr<Tensor>> out;
  if (target_) out.push_back(target_);
  out.push_back(tensor_);
  if (!subtree_.empty())
    out.push_back(subtree_.front()->target());
  return out;
}


vector<int> BinaryContraction::required_rdm(vector<int> orig) const {
  vector<int> out = orig;
  for (auto& i : subtree_)
    out = i->required_rdm(out);

  sort(out.begin(), out.end());
  return out;
}

vector<int> Tree::required_rdm(vector<int> orig) const {
  vector<int> out = orig;
  for (auto& i : bc_)
    out = i->required_rdm(out);

  for (auto& i : op_) {
    if (i->label().find("Gamma") == string::npos) continue;
    vector<int> rdmn = i->active()->required_rdm();
    for (auto& j: rdmn)
      if (find(out.begin(), out.end(), j) == out.end()) out.push_back(j);
  }
  sort(out.begin(), out.end());
  return out;
}



