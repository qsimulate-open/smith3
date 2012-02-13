//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "tree.h"
#include <algorithm>

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
  cout << target_->str() << " = " << tensor_->str() << " * " << subtree_.front()->target()->str() << endl;
  for (auto i = subtree_.begin(); i != subtree_.end(); ++i) {
    (*i)->print();
  }
}


Tree::Tree(const shared_ptr<ListTensor> l) {
  target_ = l->target();
  if (l->length() > 1) {
    shared_ptr<BinaryContraction> bc(new BinaryContraction(target_, l)); 
    bc_.push_back(bc);
  } else {
    shared_ptr<Tensor> t = l->front();
    t->set_factor(l->fac());
    op_.push_back(t); 
  }
}


Tree::Tree(shared_ptr<Equation> eq) {

  // First make ListTensor for all the diagrams
  list<shared_ptr<Tree> > lt;
  list<shared_ptr<Diagram> > d = eq->diagram();
  for (auto i = d.begin(); i != d.end(); ++i) {
    shared_ptr<ListTensor> tmp(new ListTensor(*i));
    // All internal tensor should be included in the active part
    tmp->absorb_all_internal();

    // convert it to tree
    shared_ptr<Tree> tr(new Tree(tmp));
    lt.push_back(tr);
     
    tr->print();
  }

  // Second, make a nested tree
#if 0
  for (auto i = lt.begin(); i != lt.end(); ++i) {
    (*i)->print();
    done.push_back(i);

    // first make BinaryContraction object
    shared_ptr<Tensor> first = (*i)->front();
    shared_ptr<ListTensor> rest = (*i)->rest();
    shared_ptr<BinaryContraction> b(new BinaryContraction());

    auto j = i; ++j;
    for (; j != lt.end(); ++j) {
      // matching here; first tensor is the same + dagger_ is the same: then factorizable
      if (*first == *(*j)->front() && dag == (*j)->dagger()) {
        cout << "aaaa" << endl;
        done.push_back(j);
      }
    }
  }
#endif

}



void Tree::print() const {
  for (auto i = op_.begin(); i != op_.end(); ++i)
    cout << target_->str() << "+= " << (*i)->str() << endl;
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

