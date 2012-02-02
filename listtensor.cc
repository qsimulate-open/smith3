//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "listtensor.h"

using namespace std;

ListTensor::ListTensor(shared_ptr<Diagram> d) {
  // factor
  fac_ = d->fac();
  // vector of tensors
  for (auto i = d->op().begin(); i != d->op().end(); ++i) {
    shared_ptr<Tensor> t(new Tensor(*i));
    list_.push_back(t);
  }
  if (d->rdm()) {
    shared_ptr<Tensor> t(new Tensor(d->rdm()));
    list_.push_back(t);
  }
  // dag
  dagger_ = d->dagger();
}


void ListTensor::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto i = list_.begin(); i != list_.end(); ++i) cout << (*i)->str();
  if (dagger_) cout << " ** Daggered object added **";
  cout << endl;
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    if ((*i)->active()) (*i)->active()->print("   ");
  }

  cout << endl;
  
}


void ListTensor::absorb_all_internal() {
  auto j = list_.begin();
  // first find active
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    if ((*i)->active()) j = i;
  }
  list<list<shared_ptr<Tensor> >::iterator> remove;
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    if ((*i)->all_active() && !(*i)->active()) {
      (*j)->merge(*i);
      remove.push_back(i);
    }
  }
  for (auto i = remove.begin(); i != remove.end(); ++i) list_.erase(*i);
}


