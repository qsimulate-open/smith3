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
  list_ = d->op();
  // active
  active_ = d->rdm();
  // dag
  dagger_ = d->dagger();
}


void ListTensor::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto i = list_.begin(); i != list_.end(); ++i) (*i)->print();

  // active operators
  cout << "[";
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    shared_ptr<Op> o = *i;
    if (o->num_active_nodagger() + o->num_active_dagger() != 0) {
      for (auto j = o->op().begin(); j != o->op().end(); ++j) {
        if (get<1>(*j) == -1 || get<1>(*j) == 0) continue;
        cout << (*get<0>(*j))->str();
      }
    }
  }
  cout << "]";
  if (dagger_) cout << " ** Daggered object added **";

  cout << endl;
  if (active_) active_->print("   ");
  cout << endl;
  
}

