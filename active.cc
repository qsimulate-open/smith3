//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "active.h"

using namespace std;

void RDM::print() const {

  cout << setw(5) << setprecision(2) << fac_ << " [";
  for (auto i = index_.begin(); i != index_.end(); ++i)
    cout << (*i)->str();
  for (auto i = delta_.begin(); i != delta_.end(); ++i)
    cout << " d(" << i->first->str(false) << i->second->str(false) << ")";
  cout << "]";

  if (done()) cout << "*";

  cout << endl;
}


// An application of "Wick's theorem"
list<shared_ptr<RDM> > RDM::reduce_one(list<shared_ptr<Index> >& done) const {
  // first find non-daggered operator which is not aligned
  list<shared_ptr<RDM> > out;

  for (auto i = index_.begin(); i != index_.end(); ++i) {

    // looking for a non-daggered index that is not registered in done
    if ((*i)->dagger() || find(done.begin(), done.end(), *i) != done.end())
      continue;

    for (auto j = i; j != index_.end(); ++j) {
      if (!(*j)->dagger() || j==i) continue;
      // if you find daggered object in the right hand side...
      list<shared_ptr<Index> > indx = index_;

      // find the indices to be deleted.
      vector<list<shared_ptr<Index> >::iterator> rml;
      for (auto k = indx.begin(); k != indx.end(); ++k)
        if (*k == *i || *k == *j) rml.push_back(k);
      assert(rml.size() == 2);

      list<pair<shared_ptr<Index>, shared_ptr<Index> > > dl = delta_;
      dl.push_back(make_pair(*rml[0],*rml[1]));
      double f = fac_;
      // Please note that this procedure does not change the sign (you can prove it in 30sec)
      if ((*i)->same_spin(*j)) {
        f *= 2.0;
      } else {
        // this case we need to replace a spin 
        const shared_ptr<Spin> s0 = (*rml[0])->spin();
        const shared_ptr<Spin> s1 = (*rml[1])->spin();
        for (auto k = indx.begin(); k != indx.end(); ++k) {
          if ((*k)->spin() == s0) {
//          (*k)->set_spin(s1);
          }
        }
      }

      // erasing indices which are push-backed in delta
      indx.erase(rml[0]);
      indx.erase(rml[1]);
      // create an RDM object
      shared_ptr<RDM> tmp(new RDM(indx,dl,f));
      out.push_back(tmp);
    }
    done.push_back(*i);
    break;
  }
  return out;
}


bool RDM::reduce_done(const list<shared_ptr<Index> >& done) const {
  // check if there is a annihilation operator which has creation operators in his right side
  bool out = true;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    // if non-dagger and not registered in done
    if (!(*i)->dagger() && find(done.begin(), done.end(), (*i)) == done.end()) {
      for (auto j = i; j != index_.end(); ++j) {
        if ((*j)->dagger()) out = false;
      } 
      break;
    }
  }
  return out;
}


bool RDM::done() const {
  // if operators are aligned as a0+ a1+ .. a1 a0

  bool out = true;
  // ann = false when we encounter annihilation operator:w
  bool ann = true;
  vector<shared_ptr<Spin> > dag;
  vector<shared_ptr<Spin> > nodag;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    bool idag = (*i)->dagger();
    ann &= idag; 
    if (idag) {
      if (!ann) {
        out = false;
        break;
      }
      dag.push_back((*i)->spin());
    } else if (!idag) {
      nodag.push_back((*i)->spin());
    }
  }
  out &= dag.size() == nodag.size();
  if (out) {
    auto j = dag.rbegin();
    for (auto i = dag.begin(); i != dag.end(); ++i, ++j) {
      out &= *i == *j;
    }
  }
  return out;
}


Active::Active(const list<shared_ptr<Index> >& in) {

  list<pair<shared_ptr<Index>, shared_ptr<Index> > > t;
  shared_ptr<RDM> tmp(new RDM(in, t, 1.0));

  // this sets list<RDM>
  cout << "====" << endl;
  tmp->print();
  reduce(tmp);

}


void Active::reduce(shared_ptr<RDM> in) {
  list<shared_ptr<Index> > d;
  list<std::pair<shared_ptr<RDM>, list<shared_ptr<Index> > > > buf(1, make_pair(in,d));

  while (buf.size() != 0) {
    list<std::pair<shared_ptr<RDM>, list<shared_ptr<Index> > > > buf2;

    for (auto iter = buf.begin(); iter != buf.end(); ++iter) { 
      shared_ptr<RDM> tmp = iter->first;
      list<shared_ptr<Index> > done = iter->second; 
      
      // taking delta
      list<shared_ptr<RDM> > out = tmp->reduce_one(done);
      // this is also needed!
      out.push_back(tmp);
      for (auto i = out.begin(); i != out.end(); ++i) {
        if ((*i)->reduce_done(done)) {
          rdm_.push_back(*i);
        } else {
          buf2.push_back(make_pair(*i,done));
        }
      }
    }
    buf = buf2;
  }
}


void Active::print() const {
  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) (*i)->print();
}
