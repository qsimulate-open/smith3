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
#include "constants.h"
#include <algorithm>
#include <stdexcept>

using namespace std;
using namespace smith;


shared_ptr<Tensor> BinaryContraction::next_target() {
  return subtree().front()->target();
}

BinaryContraction::BinaryContraction(shared_ptr<Tensor> o, shared_ptr<ListTensor> l) {
  // o is a target tensor NOT included in l
  target_ = o;
  tensor_ = l->front();
  shared_ptr<ListTensor> rest = l->rest();
  shared_ptr<Tree> tr(new Tree(rest));
  subtree_.push_back(tr);
}


void BinaryContraction::print() const {
  string indent = "";
  for (int i = 0; i != depth(); ++i) indent += "  ";
  cout << indent << (target_ ? (target_->str()+" = ") : "") << tensor_->str() << " * " << subtree_.front()->target()->str() << endl;
  for (auto& i : subtree_) i->print();
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


bool BinaryContraction::dagger() const {
  return subtree_.front()->dagger();
}


Tree::Tree(const shared_ptr<ListTensor> l) : num_(-1) {
  target_ = l->target();
  dagger_ = l->dagger();
  if (l->length() > 1) {
    shared_ptr<BinaryContraction> bc(new BinaryContraction(target_, l));
    bc_.push_back(bc);
  } else {
    shared_ptr<Tensor> t = l->front();
    t->set_factor(l->fac());
    t->set_scalar(l->scalar());
    op_.push_back(t);
  }
}


Tree::Tree(shared_ptr<Equation> eq) : parent_(NULL), tree_name_(eq->name()), num_(-1) {

  // First make ListTensor for all the diagrams
  list<shared_ptr<Diagram> > d = eq->diagram();
  for (auto i = d.begin(); i != d.end(); ++i) {
    shared_ptr<ListTensor> tmp(new ListTensor(*i));
    // All internal tensor should be included in the active part
    tmp->absorb_all_internal();   

    shared_ptr<Tensor> first = tmp->front();
    shared_ptr<ListTensor> rest = tmp->rest();

    // convert it to tree
    shared_ptr<Tree> tr(new Tree(rest));
    list<shared_ptr<Tree> > lt; lt.push_back(tr);
    shared_ptr<BinaryContraction> b(new BinaryContraction(lt, first));
    bc_.push_back(b);
  }

  factorize();
  set_parent_sub();
  set_target_rec();
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


bool Tree::done() const {
  // vectensor should be really a tensor
#if 0
  bool out = true;
  for (auto i = tensor_.begin(); i != tensor_.end(); ++i) out &= ((*i)->length() == 1);
  return out;
#endif
// todo check this
 cout << "Check Tree::done" << endl;
 return true;
}


void BinaryContraction::factorize() {
  list<list<shared_ptr<Tree> >::iterator> done;
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


void Tree::sort_gamma() {
  cout << "debugging sort_gamma " << endl; 

  list<shared_ptr<Tensor> > g = gammas();
  for (auto& i : g) {
    i->print();
    // find like tensors in list using overloaded ==
    find_gamma(i);
  }
  cout << "\nUnique Gamma Tensors: " << endl;
  for (auto j = gamma_.begin(); j != gamma_.end(); ++j)
    cout << "gamma_ list: " << (*j)->label() << (j == --gamma_.end() ? "\n" : "") << endl;
}


void Tree::find_gamma(shared_ptr<Tensor> o) {
  bool found = false;
  for (auto& i : gamma_) {
    if ((*i) == (*o)) {
      cout << "Rename: " << o->label() <<  " to " << i->label() <<  endl;
      found = true;
      break;
    }
  }
  if (!found) gamma_.push_back(o);
}


// return gamma below this
list<shared_ptr<Tensor> > Tree::gammas() const {
  list<shared_ptr<Tensor> > out;
  for (auto& i : bc_) {
    for (auto& j : i->subtree()) {
      // recursive call 
      list<shared_ptr<Tensor> > tmp = j->gammas();
      out.insert(out.end(), tmp.begin(), tmp.end());
    }
  }
  for (auto& i : op_)
    if (i->label().find("Gamma") != string::npos) out.push_back(i); 
  return out;
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


void Tree::factorize() {
  list<list<shared_ptr<BinaryContraction> >::iterator> done;
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    auto j = i; ++j;
    for ( ; j != bc_.end(); ++j) {
      if (*(*i)->tensor() == *(*j)->tensor() && (*i)->dagger() == (*j)->dagger()) {
        done.push_back(j);
        (*i)->subtree().insert((*i)->subtree().end(), (*j)->subtree().begin(), (*j)->subtree().end());
      }
    }
  }
  for (auto& i : done) bc_.erase(i);
  for (auto& i : bc_) i->factorize();
}


int BinaryContraction::depth() const { return parent_->depth(); }
int Tree::depth() const { return parent_ ? parent_->parent()->depth()+1 : 0; }


void Tree::set_target_rec() {
  if (!bc_.empty()) {
    for (auto& i : bc_) i->set_target_rec();
  }
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


static int icnt;
// not a good implementation...
vector<shared_ptr<Tensor> > ii;

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
    if (label == "f1" || label == "v2") label = "this->" + label + "_";
    ss << (label != array.front() ? ", " : "") << ((label == "proj") ? "r" : label);
  }
  return ss.str();
}
static string merge__(list<string> array) { return merge__(vector<string>(array.begin(), array.end())); }
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




list<shared_ptr<Index> > BinaryContraction::target_indices() {
  // returns a list of target indices
  return target_->index();
}


list<shared_ptr<Index> > BinaryContraction::loop_indices() {
  // returns a list of inner loop indices.
  list<shared_ptr<Index> > out;
  list<shared_ptr<Index> > ti = target_->index();
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


string Tree::generate_compute_header(const int ic, const vector<shared_ptr<Tensor> > tensors, const bool enlist) const {
  stringstream tt;
  tt << "template <typename T>" << endl;
  if (!enlist) {
    tt << "class Task" << ic << " : public Task<T> {" << endl;
  } else {
    tt << "class Task" << ic << " : public EnergyTask<T> {" << endl;
  }
  tt << "  protected:" << endl;
  tt << "    IndexRange closed_;" << endl;
  tt << "    IndexRange active_;" << endl;
  tt << "    IndexRange virt_;" << endl;

  vector<string> done;
  for (auto& s : tensors) {
    string label = label__(s->label());
    if (find(done.begin(), done.end(), label) != done.end()) continue;
    done.push_back(label);
    tt << "    std::shared_ptr<Tensor<T> > " << label << ";" << endl;
    // check if scalar needs to be added  
    if (!s->scalar().empty()) {
      if (s->scalar() == "e0" ) {
        tt << "    double " << s->scalar() << "_;" << endl;
      } else { 
        throw logic_error ("Implimentation needed");
      }
    }
  }

  tt << "" << endl;
  tt << "    void compute_() {" << endl;
  if (enlist)
    tt << "      this->energy_ = 0.0;" << endl;
  return tt.str();
}


string Tree::generate_compute_footer(const int ic, const vector<shared_ptr<Tensor> > tensors, const bool enlist) const {
  stringstream tt;
  tt << "    };" << endl;
  tt << "" << endl;
  tt << "  public:" << endl;
  if (!enlist) {
    bool found = false;
    for (auto& s : tensors)
      if (!s->scalar().empty()) found = true;
    tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i" << (found ? ", const double e" : "") << ") : Task<T>() {" << endl;
  } else {
    tt << "    Task" << ic << "(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {" << endl;
  }
  tt << "      closed_ = i[0];" << endl;
  tt << "      active_ = i[1];" << endl;
  tt << "      virt_   = i[2];" << endl;

  int i = 0;
  vector<string> done;
  for (auto s = tensors.begin(); s != tensors.end(); ++s) {
    string label = label__((*s)->label());

    if (find(done.begin(), done.end(), label) != done.end()) continue;
    done.push_back(label);
    tt << "      " << (label == "proj" ? "r" : label) << " = t[" << i << "];" << endl;
    ++i;

    // check if scalar needs to be added  
    if (!(*s)->scalar().empty()) {
      if ((*s)->scalar() == "e0") {
        tt << "      e0_ = e; " << endl;
      } else { 
        throw logic_error("Tree::generate_compute_footer, only \"e0\" is accepted so far");
      }
    }
  }


  tt << "    };" << endl;
  tt << "    ~Task" << icnt << "() {};" << endl;
  tt << "};" << endl << endl;
  return tt.str();
}


string Tree::generate_gamma(const int ic, const shared_ptr<Tensor> gamma, const bool enlist) const {
  stringstream tt;
  vector<int> rdmn = gamma->active()->required_rdm();

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
  tt << "    std::shared_ptr<Tensor<T> > " << gamma->label() << ";" << endl;
  for (auto& i: rdmn)
      tt << "    std::shared_ptr<Tensor<T> > rdm" << i << ";" << endl;
  if (gamma->merged()) 
    tt << "    std::shared_ptr<Tensor<T> > " << gamma->merged()->label() << ";" << endl;
  tt << endl;
  // loops
  tt << "    void compute_() {" << endl;
  string indent ="      ";
  vector<string> close;
  tt << gamma->generate_gamma(indent, close, "o", true); 

  // generate gamma put block
  tt << indent << gamma->label() << "->put_block(ohash, odata);" << endl;
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
  tt << "      " << gamma->label() << "  = t[0];" << endl;
  //this is related to what rdms we have
  {
    int itmp = 1;
    for (auto i = rdmn.begin(); i != rdmn.end(); ++i, ++itmp) {
      tt << "      rdm" << (*i) << "    = t[" << itmp << "];" << endl;
    }
    // related to merged tensor
    if (gamma->merged()) tt << "      "<< gamma->merged()->label() << "      = t[" <<  itmp << "];" << endl;
  }
  tt << "    };" << endl;
  tt << "    ~Task" << icnt << "() {};" << endl;
  tt << "};" << endl << endl;
  return tt.str();
}


string Tree::generate_compute_operators(const string indent, const shared_ptr<Tensor> target, const list<shared_ptr<Tensor> > op,
                                        const bool dagger) const {
  stringstream tt;

  vector<string> close;
  string cindent = indent;
  // note that I am using reverse_iterator
  tt << target->generate_loop(cindent, close);
  // get data
  tt << target->generate_get_block(cindent, "o", true);

  // add the source data to the target
  int j = 0;
  for (auto s = op.begin(); s != op.end(); ++s, ++j) {
    stringstream uu; uu << "i" << j;
    tt << cindent << "{" << endl;
    tt << (*s)->generate_get_block(cindent+"  ", uu.str());
    list<shared_ptr<Index> > di = target->index();
    tt << (*s)->generate_sort_indices(cindent+"  ", uu.str(), di, true);
    tt << cindent << "}" << endl;

    // if this is at the top-level and nees to be daggered:
    if (depth() == 0 && dagger) {
      shared_ptr<Tensor> top(new Tensor(**s));
      // swap operators so that tensor is daggered
      if (top->index().size() != 4) {
         throw logic_error("Daggered object is only supported for 4-index tensors");
      } else {
        auto k0 = di.begin(); auto k1 = k0; ++k1; auto k2 = k1; ++k2; auto k3 = k2; ++k3;
        list<pair<list<shared_ptr<Index> >::iterator, list<shared_ptr<Index> >::iterator> > map;
        map.push_back(make_pair(k0, k2));
        map.push_back(make_pair(k2, k0));
        map.push_back(make_pair(k1, k3));
        map.push_back(make_pair(k3, k1));
        list<shared_ptr<Index> > tmp;
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
        tt << top->generate_get_block(cindent+"  ", uu.str());
        list<shared_ptr<Index> > di = target->index();
        tt << top->generate_sort_indices(cindent+"  ", uu.str(), di, true);
        tt << cindent << "}" << endl;
      }
    }
  }

  // put data
  {
    string label = target->label();
    tt << cindent << (label == "proj" ? "r" : label) << "->put_block(ohash, odata);" << endl;
  }
  for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
        tt << *iter << endl;
  return tt.str();
}


string Tree::generate_task(const string indent, const int ic, const vector<shared_ptr<Tensor> > op, const bool enlist) const {
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

  return generate_task(indent, ip, ic, ops, enlist, scalar);
}


string Tree::generate_task(const string indent, const int ip, const int ic, const vector<string> op, const bool enlist, const string scalar) const {
  stringstream ss;
  ss << indent << "std::vector<std::shared_ptr<Tensor<T> > > tensor" << ic << " = {" << merge__(op) << "};" << endl;
  ss << indent << "std::shared_ptr<Task" << ic << "<T> > task"
               << ic << "(new Task" << ic << "<T>(tensor" << ic << ", index" << (scalar.empty() ? "" : ", this->e0_") << "));" << endl;

  if (!enlist) {
    if (parent_) {
      assert(parent_->parent());
      ss << indent << "task" << ip << "->add_dep(task" << ic << ");" << endl;
      ss << indent << "task" << ic << "->add_dep(task0);" << endl;
    } else {
      assert(depth() == 0);
      ss << indent << "task" << ic << "->add_dep(task0);" << endl;
    }
    ss << indent << "queue_->add_task(task" << ic << ");" << endl;
  } else {
    if (parent_) {
      if (ip != ic)
        ss << indent << "task" << ip << "->add_dep(task" << ic << ");" << endl;
    } else {
      assert(depth() == 0);
    }
    ss << indent << "energy_->add_task(task" << ic << ");" << endl;
  }
  ss << endl;
  return ss.str();
}


pair<string, string> Tree::generate_task_list(const bool enlist, const shared_ptr<Tree> energy) const {
  // ss is the X.h driver routine
  stringstream ss;
  // tt is the X_tasks.h driver routine
  stringstream tt;

  // just to make sure
  if (energy && depth() != 0) throw logic_error("Energy tree passed in a wrong place.");

  string indent = "      ";
  string vectensor = "std::vector<std::shared_ptr<Tensor<T> > >";
  // if this is the top tree, we want to initialize index, as well as to create a task that zeros out the residual vector
  if (depth() == 0) {
    assert(icnt == 0 && tree_name_ != "");
    ss << header(tree_name_);
    tt << header(tree_name_ + "_tasks");

    ss << "#ifndef __SRC_SMITH_" << tree_name_ << "_H " << endl;
    ss << "#define __SRC_SMITH_" << tree_name_ << "_H " << endl;
    ss << "" << endl;
    ss << "#include <src/smith/spinfreebase.h>" << endl;
    ss << "#include <src/scf/fock.h>" << endl;
    ss << "#include <src/util/f77.h>" << endl;
    ss << "#include <iostream>" << endl;
    ss << "#include <iomanip>" << endl;
    ss << "#include <src/smith/queue.h>" << endl;
    ss << "#include <src/smith/" << tree_name_ << "_tasks.h>" << endl;
    ss << "#include <src/smith/smith_info.h>" << endl;
    ss << "" << endl;
    ss << "namespace bagel {" << endl;
    ss << "namespace SMITH {" << endl;
    ss << "namespace " << tree_name_ << "{" << endl;
    ss << "" << endl;
    ss << "template <typename T>" << endl;
    ss << "class " << tree_name_ << " : public SpinFreeMethod<T>, SMITH_info {" << endl;
    ss << "  protected:" << endl;
    ss << "    std::shared_ptr<Tensor<T> > t2;" << endl;
    ss << "    std::shared_ptr<Tensor<T> > r;" << endl;
    ss << "    double e0_;" << endl;
    ss << "" << endl;
    ss << "    std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > make_queue_() {" << endl;
    ss << "      std::shared_ptr<Queue<T> > queue_(new Queue<T>());" << endl;

    ss << indent << "std::vector<IndexRange> index = {this->closed_, this->active_, this->virt_};" << endl << endl;
    ss << indent << vectensor << " tensor0 = {r};" << endl;
    ss << indent << "std::shared_ptr<Task0<T> > task0(new Task0<T>(tensor0, index));" << endl;
    ss << indent << "queue_->add_task(task0);" << endl << endl;

    tt << "#ifndef __SRC_SMITH_" << tree_name_ << "_TASKS_H " << endl;
    tt << "#define __SRC_SMITH_" << tree_name_ << "_TASKS_H " << endl;
    tt << "" << endl;
    tt << "#include <memory>" << endl;
    tt << "#include <algorithm>" << endl;
    tt << "#include <src/smith/indexrange.h>" << endl;
    tt << "#include <src/smith/tensor.h>" << endl;
    tt << "#include <src/smith/task.h>" << endl;
    tt << "#include <vector>" << endl;
    tt << "" << endl;
    tt << "namespace bagel {" << endl;
    tt << "namespace SMITH {" << endl;
    tt << "namespace " << tree_name_ << "{" << endl;
    tt << "" << endl;

    // here generate Task0 that zeros out the residual
    tt << "template <typename T>" << endl;
    tt << "class Task0 : public Task<T> {" << endl;
    tt << "  protected:" << endl;
    tt << "    std::shared_ptr<Tensor<T> > r_;" << endl;
    tt << "    IndexRange closed_;" << endl;
    tt << "    IndexRange active_;" << endl;
    tt << "    IndexRange virt_;" << endl;
    tt << "" << endl;
    tt << "    void compute_() {" << endl;
    tt << "      r_->zero();" << endl;
    tt << "    };  " << endl;
    tt << "" << endl;
    tt << "  public:" << endl;
    tt << "    Task0(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {" << endl;
    tt << "      r_ =  t[0];" << endl;
    tt << "      closed_ = i[0];" << endl;
    tt << "      active_ = i[1];" << endl;
    tt << "      virt_   = i[2];" << endl;
    tt << "    };  " << endl;
    tt << "    ~Task0() {}; " << endl;
    tt << "};" << endl << endl;

    ++icnt;
  }

  /////////////////////////////////////////////////////////////////
  // if op_ is not empty, we add a task that adds up op_.
  /////////////////////////////////////////////////////////////////
  if (!op_.empty()) {

    // Gamma should be constructed here
    for (auto& i : op_)
      if (i->label().find("Gamma") != string::npos) {
        // reconstruct gamma label if necessary
        ss << i->constructor_str(indent) << endl; 
      }
    
    // step through operators and if they are new, construct them.
    if (find(ii.begin(), ii.end(), target_) == ii.end()) {
      ii.push_back(target_);
      ss << target_->constructor_str(indent) << endl;
    }

    vector<shared_ptr<Tensor> > op = {{target_}};
    op.insert(op.end(), op_.begin(), op_.end());
    ss << generate_task(indent, icnt, op, enlist);

    vector<shared_ptr<Tensor> > op2(op.begin(), op.end());
    tt << generate_compute_header(icnt, op2, enlist);
    tt << generate_compute_operators(indent, target_, op_);
    tt << generate_compute_footer(icnt, op2, enlist);

    ++icnt;

    for (auto& i : op_) {
      if (i->label().find("Gamma") != string::npos) {
        tt << generate_gamma(icnt, i, enlist);
        vector<string> tmp = {i->label()};
        vector<int> rdms = i->active()->required_rdm(); 
        for (auto& j : rdms) {
          stringstream zz; 
          zz << "this->rdm" << j << "_";
          tmp.push_back(zz.str());
        }
        if (i->merged()) {
          stringstream mm;
          mm << "this->" << i->merged()->label() << "_";
          tmp.push_back(mm.str());
        }
        ss << generate_task(indent, icnt-1, icnt, tmp, enlist);
        ++icnt;
      }
    }

  }

  /////////////////////////////////////////////////////////////////
  // step through BinaryContraction
  /////////////////////////////////////////////////////////////////
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    vector<shared_ptr<Tensor> > strs = (*i)->tensors_str();

    for (auto& s : strs) {
      // if it contains a new intermediate tensor, dump a constructor
      // (while registering in ii - note that ii is a static variable here).
      if (find(ii.begin(), ii.end(), s) == ii.end() && s->label().find("I") != string::npos) {
        ii.push_back(s);
        ss << s->constructor_str(indent) << endl;
      }
    }
    // saving a counter to a protected member for dependency checks
    num_ = icnt;
    ss << generate_task(indent, icnt, strs, enlist);

    // here we generate a task for this binary contraction
    tt << generate_compute_header(num_, strs, enlist);

    if (depth() != 0) {
      vector<string> close;
      string cindent = indent;
      list<shared_ptr<Index> > ti = (*i)->target_indices();
      // note that I am using reverse_iterator
      // skip if this is the last step in energy contribution
      if (!enlist || depth() != 1) {
        for (auto iter = ti.begin(); iter != ti.end(); ++iter, cindent += "  ") {
          string cindex = (*iter)->str_gen();
          tt << cindent << "for (auto& " << cindex << " : " << (*iter)->generate() << ") {" << endl;
          close.push_back(cindent + "}");
        }

        tt << target_->generate_get_block(cindent, "o", true);
        tt << target_->generate_scratch_area(cindent, "o", true); // true means zero-out
      }

      // inner loop will show up here
      tt << endl;
      vector<string> close2;
      string dindent = cindent;
      list<shared_ptr<Index> > di = (*i)->loop_indices();
      for (auto iter = di.rbegin(); iter != di.rend(); ++iter, dindent += "  ") {
        string cindex = (*iter)->str_gen();
        tt << dindent << "for (auto& " << cindex << " : " << (*iter)->generate() << ") {" << endl;
        close2.push_back(dindent + "}");
      }

      // retrieving tensor_
      tt << (*i)->tensor()->generate_get_block(dindent, "i0");
      tt << (*i)->tensor()->generate_sort_indices(dindent, "i0", di) << endl;
      // retrieving subtree_
      tt << (*i)->next_target()->generate_get_block(dindent, "i1");
      tt << (*i)->next_target()->generate_sort_indices(dindent, "i1", di) << endl;

      // call dgemm
      {
        pair<string, string> t0 = (*i)->tensor()->generate_dim(di);
        pair<string, string> t1 = (*i)->next_target()->generate_dim(di);
        if (t0.first != "" || t1.first != "") {
          tt << dindent << "dgemm_(\"T\", \"N\", ";
          string tt0 = t0.first == "" ? "1" : t0.first;
          string tt1 = t1.first == "" ? "1" : t1.first;
          string ss0 = t1.second== "" ? "1" : t1.second;
          tt << tt0 << ", " << tt1 << ", " << ss0 << "," << endl;
          tt << dindent << "       1.0, i0data_sorted, " << ss0 << ", i1data_sorted, " << ss0 << "," << endl
             << dindent << "       1.0, odata_sorted, " << tt0;
          tt << ");" << endl;
        } else {
          // so far I am expecting the case of energy contribution
          if (!enlist || depth() != 1) throw logic_error("so far I am expecting the case of energy contribution");
          string ss0 = t1.second== "" ? "1" : t1.second;
          tt << dindent << "this->energy_ += ddot_(" << ss0 << ", i0data_sorted, 1, i1data_sorted, 1);" << endl;
        }
      }

      for (auto iter = close2.rbegin(); iter != close2.rend(); ++iter)
        tt << *iter << endl;
      tt << endl;
      // Inner loop ends here

      // skip if this is the last step in energy contribution
      if (!enlist || depth() != 1) {
        // sort buffer
        {
          tt << (*i)->target()->generate_sort_indices_target(cindent, "o", di, (*i)->tensor(), (*i)->next_target());
        }
        // put buffer
        {
          string label = target_->label();
          tt << cindent << (label == "proj" ? "r" : label) << "->put_block(ohash, odata);" << endl;
        }
        for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
          tt << *iter << endl;
      }

    } else {
      // making residual vector...
      list<shared_ptr<Index> > proj = bc_.front()->tensor()->index();
      list<shared_ptr<Index> > res;
      assert(!(proj.size() & 1));
      for (auto ii = proj.begin(); ii != proj.end(); ++ii, ++ii) {
        auto j = ii; ++j;
        res.push_back(*j);
        res.push_back(*ii);
      }
      shared_ptr<Tensor> residual(new Tensor(1.0, "r", res));
      list<shared_ptr<Tensor> > op2(1, (*i)->next_target());
      tt << generate_compute_operators(indent, residual, op2, (*i)->dagger());
    }

    tt << generate_compute_footer(num_, strs, enlist);


    // increment icnt before going to subtrees
    ++icnt;
    // triggers a recursive call
    pair<string, string> tmp = (*i)->generate_task_list(enlist);
    ss << tmp.first;
    tt << tmp.second;
  }

  if (depth() == 0) {
    // energy queue here
    ss << "      std::shared_ptr<Queue<T> > energy_(new Queue<T>());" << endl;
    if (energy) {
      energy->num_ = icnt;
      for (auto& i : energy->bc_) {
        pair<string, string> tmp = i->generate_task_list(true); // true triggers energy list
        ss << tmp.first;
        tt << tmp.second;
      }
      ++icnt;
    }
    ss << "      return make_pair(queue_, energy_);" << endl;
    ss << "    };" << endl;
    ss << endl;
    ss << "  public:" << endl;
    ss << "    " << tree_name_ << "(std::shared_ptr<const Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info() {" << endl;
    ss << "      this->eig_ = this->f1_->diag();" << endl;
    ss << "      t2 = this->v2_->clone();" << endl;
    ss << "      e0_ = this->compute_e0();" << endl;
    ss << "#if 1" << endl;
    ss << "      this->update_amplitude_start(t2, this->v2_);" << endl;
    ss << "      t2->scale(2.0);" << endl;
    ss << "#endif" << endl;
    ss << "      r = t2->clone();" << endl;
    ss << "    };" << endl;
    ss << "    ~" << tree_name_ << "() {}; " << endl;
    ss << "" << endl;
    ss << "    void solve() {" << endl;
    ss << "      this->print_iteration();" << endl;
    ss << "      int iter = 0;" << endl;
    ss << "      for ( ; iter != maxiter_; ++iter) {" << endl;
    ss << "        std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > q = make_queue_();" << endl;
    ss << "        std::shared_ptr<Queue<T> > queue = q.first;" << endl;
    ss << "        std::shared_ptr<Queue<T> > energ = q.second;" << endl;
    ss << "        while (!queue->done())" << endl;
    ss << "          queue->next_compute();" << endl;
    ss << "        r->scale(0.25); // FIXME" << endl;
    ss << "//      *r = *(r->add_dagger());" << endl;
    ss << "        this->update_amplitude(t2, r);" << endl;
    ss << "        const double err = r->rms();" << endl;
    ss << "        r->zero();" << endl;
    ss << "        const double en = energy(energ);" << endl;
    ss << "        this->print_iteration(iter, en, err);" << endl;
    ss << "        if (err < thresh_residual()) break;" << endl;
    ss << "      }" << endl;
    ss << "      this->print_iteration(iter == maxiter_);" << endl;
    ss << "    };" << endl;
    ss << "" << endl;
    ss << "    double energy(std::shared_ptr<Queue<T> > energ) {" << endl;
    ss << "      double en = 0.0;" << endl;
    ss << "      while (!energ->done()) {" << endl;
    ss << "        std::shared_ptr<Task<T> > c = energ->next_compute();" << endl;
    ss << "        en += c->energy() * 0.25; // FIXME" << endl;
    ss << "      }   " << endl;
    ss << "      return en; " << endl;
    ss << "    };  " << endl;
    ss << "};" << endl;
    ss << endl;
    ss << "}" << endl;
    ss << "}" << endl;
    ss << "}" << endl;

    ss << "#endif" << endl << endl;

    tt << endl;
    tt << "}" << endl;
    tt << "}" << endl;
    tt << "}" << endl;
    tt << "#endif" << endl << endl;
  }
  return make_pair(ss.str(), tt.str());
}


pair<string, string> BinaryContraction::generate_task_list(const bool enlist) const {
  stringstream ss, tt;
  for (auto& i : subtree_) {
    pair<string, string> tmp = i->generate_task_list(enlist);
    ss << tmp.first;
    tt << tmp.second;
  }
  return make_pair(ss.str(), tt.str());
}


vector<shared_ptr<Tensor> > BinaryContraction::tensors_str() {
  vector<shared_ptr<Tensor> > out;
  if (target_) out.push_back(target_);
  out.push_back(tensor_);
  if (!subtree_.empty())
    out.push_back(subtree_.front()->target());
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


vector<int> BinaryContraction::required_rdm(vector<int> orig) const {
  vector<int> out = orig;
  for (auto& i : subtree_)
    out = i->required_rdm(out);

  sort(out.begin(), out.end());
  return out;
}
