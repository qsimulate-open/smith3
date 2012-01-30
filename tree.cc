//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "tree.h"

using namespace std;

Tree::Tree(shared_ptr<Equation> eq) {

  list<shared_ptr<Diagram> > d = eq->diagram();
  for (auto i = d.begin(); i != d.end(); ++i) {
    shared_ptr<ListTensor> tmp(new ListTensor(*i));
    tensor_.push_back(tmp);
  }

}



void Tree::print() const {

  for (auto i = tensor_.begin(); i != tensor_.end(); ++i) {
    (*i)->print();
  }

}


bool Tree::done() const {
  // vectensor should be really a tensor 
  bool out = true;
  for (auto i = tensor_.begin(); i != tensor_.end(); ++i) out &= ((*i)->length() == 1);
  return out;
}

