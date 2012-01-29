//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "diagram.h"

using namespace std;


list<shared_ptr<Diagram> > Diagram::get_all() const {
  list<shared_ptr<Diagram> > out;

  const int max = 1 << num_general();
  for (int i = 0; i != max; ++i) {
    int j = i;
    shared_ptr<Diagram> d = copy();
    for (auto k = d->op().begin(); k != d->op().end(); ++k) (*k)->mutate_general(j);
    if (d->consistent_indices()) {
      out.push_back(d);
    }
  }
  return out;
}

// copying diagram with the same connectivity and so on.
shared_ptr<Diagram> Diagram::copy() const {
  // mapping of indices and spins. map<old, new>
  map<shared_ptr<Index>, shared_ptr<Index> > indexmap;
  map<shared_ptr<Spin>, shared_ptr<Spin> > spinmap;

  // creates Diagram without any info
  shared_ptr<Diagram> out(new Diagram());
  list<shared_ptr<Op> > outop;

  // loop over operators
  for (auto iter = op().begin(); iter != op().end(); ++iter) {
    // cloning...
    shared_ptr<Op> a = (*iter)->copy();

    auto j = a->op().begin();
    for (auto i = (*iter)->op().begin(); i != (*iter)->op().end(); ++i, ++j) {
      auto s = indexmap.find(*get<0>(*i));
      if (s == indexmap.end()) {
        indexmap.insert(make_pair(*get<0>(*i), *get<0>(*j)));
      } else {
        *get<0>(*j) = s->second; // one of the operators has s->second...
      }
      auto z = spinmap.find((*iter)->rho(get<2>(*i)));
      if (z == spinmap.end()) {
        spinmap.insert(make_pair((*iter)->rho(get<2>(*i)), a->rho(get<2>(*j))));
      } else {
        a->set_rho(get<2>(*j), z->second); // one of the operators has z->second...
      }
      get<1>(*j) = get<1>(*i);
    }
    outop.push_back(a);
  }
  out->set_op(outop);
  out->set_fac(fac_);
  return out;
}


void Diagram::refresh_indices() {
  map<shared_ptr<Index>, int> dict;
  map<shared_ptr<Index>, int> done;
  map<shared_ptr<Spin>, int> spin;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    (*i)->refresh_indices(dict, done, spin);
}


// this is not a const function because it refreshes the indices
void Diagram::print() {
  refresh_indices();
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto i = op_.begin(); i != op_.end(); ++i) (*i)->print();

  // active operators
  cout << "[";
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    shared_ptr<Op> o = *i;
    if (o->num_active_nodagger() + o->num_active_dagger() != 0) {
      for (auto j = o->op().begin(); j != o->op().end(); ++j) {
        if (get<1>(*j) == -1 || get<1>(*j) == 0) continue;
        cout << (*get<0>(*j))->str();
      }
    }
  }
  cout << "]";

  cout << endl;
  if (rdm_) rdm_->print("   ");
  cout << endl;
}


list<shared_ptr<Index> > Diagram::active_indices() const {
  list<shared_ptr<Index> > out;
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    shared_ptr<Op> o = *i;
    if (o->num_active_nodagger() + o->num_active_dagger() != 0) {
      for (auto j = o->op().begin(); j != o->op().end(); ++j) {
        if (get<1>(*j) == 2) out.push_back(*get<0>(*j));
      }
    }
  }
  return out;
}


void Diagram::print() const {
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto i = op_.begin(); i != op_.end(); ++i) (*i)->print();
  cout << endl;
}


bool Diagram::reduce_one_noactive(const int skip) {
  refresh_indices();

  bool found = false;
  // find the first dagger operator in list<Op>
  auto i = op_.begin();
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > data; // safe because they are held by tensors
  for (; i != op_.end(); ++i) {
    // this simultaneously eliminates one entry. op_ is therefore modified here
    data = (*i)->first_dagger_noactive();
    if (data.first) break;
  }
  if (!data.first) return false;

  // skip until it comes to skip
  int cnt = 0;
  shared_ptr<Spin> newspin, oldspin;
  for (auto j = op_.begin(); j != op_.end(); ++j) {
    // cannot contract with self
    if (i == j) continue;
    // all possible contraction pattern taken for *j (returned as a list).
    if (cnt + (*j)->num_nodagger() > skip) {
      tuple<double,shared_ptr<Spin>,shared_ptr<Spin> > tmp = (*j)->contract(data, skip-cnt);
      fac_ *= get<0>(tmp);
      newspin = get<1>(tmp);
      oldspin = get<2>(tmp);
      found = true;
      break;
    } else {
      cnt += (*j)->num_nodagger();
    }
  }
  for (auto j = op_.begin(); j != op_.end(); ++j) {
    for (auto k = (*j)->rho().begin(); k != (*j)->rho().end(); ++k) {
      if (*k == oldspin) *k = newspin;
    } 
  }
  return found;
}


bool Diagram::valid() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    if (!(*i)->contracted()) ++out;
  }
  return out > 1;
}


bool Diagram::done() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    if (!(*i)->contracted()) ++out;
  }
  return out == 0; 
}


bool Diagram::done_noactive() const {
  int out = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    out += (*i)->num_nodagger() + (*i)->num_dagger();
  }
  return out == 0; 
}


int Diagram::num_dagger() const {
  int out = 0;
  for (auto i = op_.begin(); i !=  op_.end(); ++i) out += (*i)->num_dagger();
  return out;
}


int Diagram::num_general() const {
  int cnt = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i) cnt += (*i)->num_general(); 
  return cnt;
}


bool Diagram::consistent_indices() const {
  int cnt1 = 0;
  int cnt2 = 0;
  for (auto i = op_.begin(); i != op_.end(); ++i) {
    cnt1 += (*i)->num_active_dagger(); 
    cnt2 += (*i)->num_active_nodagger(); 
  }
  return cnt1 == cnt2;
}

void Diagram::active() {
  // Index should be updated.
  refresh_indices();
  // Performs Wick in constructor of an Active object
  shared_ptr<Active> tmp(new Active(active_indices()));
  // Sets to member
  rdm_ = tmp; 
}


bool Diagram::permute() {
print();
  pair<bool, double> a = op_.back()->permute();
print();
  return false;
}




bool Diagram::identical(shared_ptr<Diagram> o) const {

  bool out = true;
  // first, they should be same size
  if (op_.size() != o->op().size()) out = false;
  // second, each indices should be the same (spin is not checked here)
  if (out) {
    for (auto i = op_.begin(), j = o->op().begin(); i != op_.end(); ++i, ++j) out &= (*i)->identical(*j);
  }
  // then, we check spins.
  if (out) {
    list<shared_ptr<Index> > act = active_indices();
    list<shared_ptr<Index> > oact = o->active_indices();
    map<shared_ptr<Spin>, shared_ptr<Spin> > myo;
    if (act.size() != oact.size()) {
      out = false;
    } else {
      for (auto i = act.begin(), j = oact.begin(); i != act.end(); ++i, ++j) {
        assert((*i)->identical(*j)); 
        shared_ptr<Spin> s = (*i)->spin();
        shared_ptr<Spin> os = (*j)->spin();
        auto iter = myo.find(s);
        if (myo.end() != iter) {
          // if s appears for the first time, register it
          myo.insert(make_pair(s,os));
        } else {
          if (os != iter->second) {
            out =false;
            break;
          }
        }
      }
    }
  } 

  return out;
}



