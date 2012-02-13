//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "tree.h"
#include <algorithm>

using namespace std;

Tree::Tree(shared_ptr<Equation> eq) {

  list<shared_ptr<ListTensor> > lt;
  // First make ListTensor for all the diagrams
  list<shared_ptr<Diagram> > d = eq->diagram();
  for (auto i = d.begin(); i != d.end(); ++i) {
    shared_ptr<ListTensor> tmp(new ListTensor(*i));
    // All internal tensor should be included in the active part
    tmp->absorb_all_internal();
    lt.push_back(tmp);
  }

  // matching first tensors and make BinaryContraction
  list<list<shared_ptr<ListTensor> >::iterator> done;
  for (auto i = lt.begin(); i != lt.end(); ++i) {
    (*i)->print();
    if (find(done.begin(), done.end(), i) != done.end()) continue;
    done.push_back(i);

    // first make BinaryContraction object
    shared_ptr<Tensor> first = (*i)->front();
    shared_ptr<ListTensor> rest = (*i)->rest();

    bool dag = (*i)->dagger(); 

    auto j = i; ++j;
    for (; j != lt.end(); ++j) {
      // matching here; first tensor is the same + dagger_ is the same
      if (*first == *(*j)->front() && dag == (*j)->dagger()) {
        cout << "aaaa" << endl;
        done.push_back(j);
      }
    }
  }

}



void Tree::print() const {
#if 0
  for (auto i = tensor_.begin(); i != tensor_.end(); ++i)
    (*i)->print();
#endif
}


bool Tree::done() const {
  // vectensor should be really a tensor 
#if 0
  bool out = true;
  for (auto i = tensor_.begin(); i != tensor_.end(); ++i) out &= ((*i)->length() == 1);
  return out;
#endif
}

