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
  cout << indent
       << (target_ ? (target_->str()+" = ") : "") << tensor_->str() << " * " << subtree_.front()->target()->str() << endl;
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    (*i)->print();
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
    op_.push_back(t); 
  }
}


Tree::Tree(shared_ptr<Equation> eq) : parent_(NULL), num_(-1), tree_name_(eq->name()) {

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
  for (auto i = done.begin(); i != done.end(); ++i)
    bc_.erase(*i);
  for (auto i = bc_.begin(); i != bc_.end(); ++i)
    (*i)->factorize();
}


int BinaryContraction::depth() const { return parent_->depth(); }
int Tree::depth() const { return parent_ ? parent_->parent()->depth()+1 : 0; }


void Tree::set_target_rec() {
  if (!bc_.empty()) {
    for (auto i = bc_.begin() ; i != bc_.end(); ++i)
      (*i)->set_target_rec();
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

static string merge__(vector<shared_ptr<Tensor> > array) {
  stringstream ss;
  for (auto i = array.begin(); i != array.end(); ++i) {
    string label = (*i)->label(); 
    if (label == "f1" || label == "v2") label = "this->" + label + "_";
    ss << (i != array.begin() ? ", " : "") << ((label == "proj") ? "r" : label);
  }
  return ss.str();
}


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


pair<string, string> Tree::generate_task_list() const {
  // ss is the driver routine
  stringstream ss;
  // tt is the driver routine
  stringstream tt;

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
    ss << "#include <src/smith/smith.h>" << endl;
    ss << "" << endl;
    ss << "namespace SMITH {" << endl;
    ss << "namespace " << tree_name_ << "{" << endl;
    ss << "" << endl;
    ss << "template <typename T>" << endl;
    ss << "class " << tree_name_ << " : public SpinFreeMethod<T>, SMITH_info {" << endl;
    ss << "  protected:" << endl;
    ss << "    std::shared_ptr<Queue<T> > queue_;" << endl;
    ss << "    std::shared_ptr<Queue<T> > energy_;" << endl;
    ss << "    std::shared_ptr<Tensor<T> > t2;" << endl;
    ss << "    std::shared_ptr<Tensor<T> > r;" << endl;
    ss << endl;
    ss << "  public:" << endl;
    ss << "    " << tree_name_ << "(std::shared_ptr<Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info(), ";
    ss <<                                                             "queue_(new Queue<T>()), energy_(new Queue<T>()) {" << endl;
    ss << "      this->eig_ = this->f1_->diag();" << endl;
    ss << "      t2 = this->v2_->clone();" << endl;
    ss << "      r = t2->clone();" << endl;

    ss << indent << "std::vector<IndexRange> index = vec(this->closed_, this->act_, this->virt_);" << endl << endl;
    ss << indent << vectensor << " tensor0 = vec(r);" << endl;
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
    tt << "namespace SMITH {" << endl;
    tt << "namespace " << tree_name_ << "{" << endl;
    tt << "" << endl;

    // here generate Task0 that zeros out the residual
    tt << "template <typename T>" << endl;
    tt << "class Task0 : public Task<T> {" << endl;
    tt << "  protected:" << endl;
    tt << "    std::shared_ptr<Tensor<T> > r2_;" << endl;
    tt << "    IndexRange closed_;" << endl;
    tt << "    IndexRange act_;" << endl;
    tt << "    IndexRange virt_;" << endl;
    tt << "" << endl;
    tt << "    void compute_() {" << endl;
    tt << "      r2_->zero();" << endl;
    tt << "    };  " << endl;
    tt << "" << endl;
    tt << "  public:" << endl;
    tt << "    Task0(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {" << endl;
    tt << "      r2_ =  t[0];" << endl;
    tt << "      closed_ = i[0];" << endl;
    tt << "      act_    = i[1];" << endl;
    tt << "      virt_   = i[2];" << endl;
    tt << "    };  " << endl;
    tt << "    ~Task0() {}; " << endl;
    tt << "};" << endl << endl;

    ++icnt;
  }
  for (auto i = bc_.begin(); i != bc_.end(); ++i) {
    vector<shared_ptr<Tensor> > strs = (*i)->tensors_str();

    for (auto s = strs.begin(); s != strs.end(); ++s) {
      // if it contains a new intermediate tensor, dump a constructor
      // (while registering in ii - note that ii is a static variable here).
      if (find(ii.begin(), ii.end(), *s) == ii.end() && (*s)->label().find("I") != string::npos) {
        ii.push_back(*s);
        ss << (*s)->constructor_str(indent) << endl;
      }
    }
    ss << indent << vectensor << " tensor" << icnt << " = vec(" << merge__(strs) << ");" << endl;
    ss << indent << "std::shared_ptr<Task" << icnt << "<T> > task" << icnt << "(new Task" << icnt << "<T>(tensor" << icnt << ", index));" << endl;
    // saving a counter to a protected member for dependency checks
    num_ = icnt;
    if (parent_) {
      assert(parent_->parent());
      ss << indent << "task" << parent_->parent()->num() << "->add_dep(task" << num() << ");" << endl;
    } else {
      assert(depth() == 0);
      ss << indent << "task0->add_dep(task" << num() << ");" << endl;
    }
    ss << indent << "queue_->add_task(task" << num() << ");" << endl;
    ss << endl;

    // here we generate a task for this binary contraction
    tt << "template <typename T>" << endl;
    tt << "class Task" << num_ << " : public Task<T> {" << endl;
    tt << "  protected:" << endl;
    tt << "    IndexRange closed_;" << endl;
    tt << "    IndexRange act_;" << endl;
    tt << "    IndexRange virt_;" << endl;
    {
      for (auto s = strs.begin(); s != strs.end(); ++s) {
        string label = (*s)->label();
        tt << "    std::shared_ptr<Tensor<T> > " << (label == "proj" ? "r" : label) << ";" << endl;
      }
    }
    tt << "" << endl;
    tt << "    void compute_() {" << endl;
    if (depth() != 0) {
      vector<string> close;
      string cindent = indent;
      list<shared_ptr<Index> > ti = (*i)->target_indices();
      // note that I am using reverse_iterator
      for (auto iter = ti.begin(); iter != ti.end(); ++iter, cindent += "  ") {
        string cindex = (*iter)->str_gen();
        tt << cindent << "for (auto " << cindex << " = " << (*iter)->generate() << ".begin(); "
                                      << cindex << " != " << (*iter)->generate() << ".end(); "
                                      << "++" << cindex << ") {" << endl;
        close.push_back(cindent + "}");
      }

      tt << target_->generate_get_block(cindent, "o", true);
      tt << target_->generate_scratch_area(cindent, "o", true); // true means zero-out

      // inner loop will show up here
      tt << endl;
      vector<string> close2;
      string dindent = cindent;
      list<shared_ptr<Index> > di = (*i)->loop_indices();
      for (auto iter = di.rbegin(); iter != di.rend(); ++iter, dindent += "  ") { 
        string cindex = (*iter)->str_gen();
        tt << dindent << "for (auto " << cindex << " = " << (*iter)->generate() << ".begin(); "
                                      << cindex << " != " << (*iter)->generate() << ".end(); "
                                      << "++" << cindex << ") {" << endl;
        close2.push_back(dindent + "}");
      }

      // retrieving tensor_
      tt << (*i)->tensor()->generate_get_block(dindent, "i0");
      tt << (*i)->tensor()->generate_sort_indices(dindent, "i0", di) << endl;
      // retrieving subtree_
      tt << (*i)->subtree().front()->target()->generate_get_block(dindent, "i1");
      tt << (*i)->subtree().front()->target()->generate_sort_indices(dindent, "i1", di) << endl;

      // call dgemm
      tt << dindent << "dgemm_(\"T\", \"N\", ";
      {
        pair<string, string> t0 = (*i)->tensor()->generate_dim(di);
        pair<string, string> t1 = (*i)->subtree().front()->target()->generate_dim(di);
        assert(t0.second == t1.second);
        string tt0 = t0.first == "" ? "1" : t0.first;
        string tt1 = t1.first == "" ? "1" : t1.first;
        string ss0 = t1.second== "" ? "1" : t1.second;
        tt << tt0 << ", " << tt1 << ", " << ss0 << "," << endl;
        tt << dindent << "       1.0, i0data_sorted, " << ss0 << ", i1data_sorted, " << ss0 << "," << endl
           << dindent << "       1.0, odata_sorted, " << tt0;
      }
      tt << ");" << endl;

      for (auto iter = close2.rbegin(); iter != close2.rend(); ++iter)
        tt << *iter << endl;
      tt << endl;
      // Inner loop ends here

      {
        string label = target_->label();
        tt << cindent << (label == "proj" ? "r" : label) << "->put_block(ohash, odata);" << endl;
      }
      for (auto iter = close.rbegin(); iter != close.rend(); ++iter)
        tt << *iter << endl;
    }
    tt << "    };" << endl;
    tt << "" << endl;

    tt << "  public:" << endl;
    tt << "    Task" << num_ << "(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {" << endl;
    tt << "      closed_ = i[0];" << endl;
    tt << "      act_    = i[1];" << endl;
    tt << "      virt_   = i[2];" << endl;
    {
      int i = 0;
      for (auto s = strs.begin(); s != strs.end(); ++s, ++i) {
        string label = (*s)->label();
        tt << "      " << (label == "proj" ? "r" : label) << " = t[" << i << "];" << endl;
      }
    }

    tt << "    };" << endl;
    tt << "    ~Task" << num_ << "() {};" << endl;
    tt << "};" << endl << endl;


    // increment icnt before going to subtrees
    ++icnt;
    // triggers a recursive call
    pair<string, string> tmp = (*i)->generate_task_list();
    ss << tmp.first;
    tt << tmp.second;
  }
  if (depth() == 0) {
    ss << "    };" << endl;
    ss << "    ~" << tree_name_ << "() {}; " << endl;
    ss << "" << endl;
    ss << "    void solve() {" << endl;
    ss << "      t2->zero();" << endl;
    ss << "      this->print_iteration();" << endl;
    ss << "      int iter = 0;" << endl;
    ss << "      for ( ; iter != maxiter_; ++iter) {" << endl;
    ss << "        queue_->initialize();" << endl;
    ss << "        while (!queue_->done())" << endl;
    ss << "          queue_->next()->compute();" << endl;
    ss << "        update_amplitude(t2, r);" << endl;
    ss << "        const double en = energy();" << endl; 
    ss << "        const double err = r->rms();" << endl;
    ss << "        this->print_iteration(iter, en, err);" << endl;
    ss << "        if (err < thresh_residual()) break;" << endl;
    ss << "      }" << endl;
    ss << "      this->print_iteration(iter == maxiter_);" << endl;
    ss << "    };" << endl;
    ss << "" << endl;
    ss << "    double energy() {" << endl;
    ss << "      double en = 0.0;" << endl;
    ss << "      energy_->initialize();" << endl;
    ss << "      while (!energy_->done()) {" << endl;
    ss << "        std::shared_ptr<Task<T> > c = energy_->next();" << endl;
    ss << "        c->compute();" << endl;
    ss << "        en += c->energy();" << endl;
    ss << "      }   " << endl;
    ss << "      return en; " << endl;
    ss << "    };  " << endl;
    ss << "};" << endl;
    ss << endl;
    ss << "}" << endl;
    ss << "}" << endl;

    ss << "#endif" << endl << endl;

    tt << endl;
    tt << "}" << endl;
    tt << "}" << endl;
    tt << "#endif" << endl << endl;
  }
  return make_pair(ss.str(), tt.str());
}

pair<string, string> BinaryContraction::generate_task_list() const {
  stringstream ss, tt;
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    pair<string, string> tmp = (*i)->generate_task_list();
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
