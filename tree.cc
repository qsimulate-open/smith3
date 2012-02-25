//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "tree.h"
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


Tree::Tree(shared_ptr<Equation> eq) : parent_(NULL), num_(-1) {

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
  for (auto i = array.begin(); i != array.end(); ++i)
    ss << (i != array.begin() ? ", " : "") << (((*i)->label() == "proj") ? "r" : (*i)->label());
  return ss.str();
}

string Tree::generate_task_list() const {
  stringstream ss;
  string indent = "       ";
  string vectensor = "std::vector<std::shared_ptr<Tensor<T> > >";
  // if this is the top tree, we want to initialize index, as well as to create a task that zeros out the residual vector
  if (depth() == 0) {
    assert(icnt == 0);
    ss << indent << "std::vector<IndexRange> index = vec(this->closed_, this->act_, this->virt_);" << endl << endl;
    ss << indent << vectensor << " tensor0 = vec(r);" << endl;
    ss << indent << "std::shared_ptr<Task0<T> > t0(new Task0<T>(tensor0, index));" << endl << endl;
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
    ss << indent << "std::shared_ptr<Task" << icnt << "<T> > t" << icnt << "(new Task" << icnt << "<T>(tensor" << icnt << ", index));" << endl;
    // saving a counter to a protected member
    num_ = icnt;
    if (parent_) {
      assert(parent_->parent());
      ss << indent << "t" << parent_->parent()->num() << "->add_dep(t" << num() << ");" << endl;
    } else {
      assert(depth() == 0);
      ss << indent << "t0->add_dep(t" << num() << ");" << endl;
    }
    ss << endl;

    // increment icnt before going to subtrees
    ++icnt;
    // triggers a recursive call
    ss << (*i)->generate_task_list();
  }
  return ss.str();
}

string BinaryContraction::generate_task_list() const {
  stringstream ss;
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i)
    ss << (*i)->generate_task_list();
  return ss.str();
}


vector<shared_ptr<Tensor> > BinaryContraction::tensors_str() {
  vector<shared_ptr<Tensor> > out;
  if (target_) out.push_back(target_);
  out.push_back(tensor_);
  if (!subtree_.empty())
    out.push_back(subtree_.front()->target());
  return out;
}
