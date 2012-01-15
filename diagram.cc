//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "diagram.h"

using namespace std;

#if 0
Diagram::Diagram(const Diagram& o) : fac_(o.fac()) {
  // mapping of indices and spins. map<old, new>
  map<shared_ptr<Index>, shared_ptr<Index> > indexmap;
  map<shared_ptr<Spin>, shared_ptr<Spin> > spinmap;

  for (auto iter = o.op().begin(); iter != o.op().end(); ++iter) {
    // calling a copy costructor
    Op a(*iter);

    auto j = a.op().begin();
    for (auto i = iter->op().begin(); i != iter->op().end(); ++i, ++j) {
      auto s = indexmap.find(*get<0>(*i));
      if (s == indexmap.end()) {
        indexmap.insert(make_pair(*get<0>(*i), *get<0>(*j)));
      } else {
        get<0>(*j) = &s->second; // one of the operators has s->second...
      }
      auto z = spinmap.find(*get<2>(*i));
      if (z == spinmap.end()) {
        spinmap.insert(make_pair(*get<2>(*i), *get<2>(*j)));
      } else {
        get<2>(*j) = &z->second; // one of the operators has s->second...
      }
    }

    op_.push_back(a);
  }

}
#endif

void Diagram::refresh_indices() {
  map<shared_ptr<Index>, int> dict;
  map<shared_ptr<Index>, int> done;
  map<shared_ptr<Spin>, int> spin;
  for (auto i = op_.begin(); i != op_.end(); ++i)
    i->refresh_indices(dict, done, spin);
}


// this is not a const function because it refreshes the indices
void Diagram::print() {
  refresh_indices();
  cout << setw(4) << setprecision(1) << fixed <<  fac_ << " ";
  for (auto i = op_.begin(); i != op_.end(); ++i) i->print();
  cout << endl;
}


bool Diagram::reduce_one_noactive(const int skip) {
  bool found = false;
  // find the first dagger operator in list<Op>
  auto i = op_.begin();
  pair<shared_ptr<Index>*, shared_ptr<Spin>* > data; // safe because they are held by tensors
  for (; i != op_.end(); ++i) {
    // this simultaneously eliminates one entry. op_ is therefore modified here
    data = i->first_dagger_noactive();
    if (data.first) break;
  }
  if (!data.first) return false;

  // skip until it comes to skip
  int cnt = 0;
  for (auto j = op_.begin(); j != op_.end(); ++j) {
    // cannot contract with self
    if (i == j) continue;
    // all possible contraction pattern taken for *j (returned as a list).
    if (cnt + j->num_nodagger() > skip) {
      fac_ *= j->contract(data, skip-cnt);
      found = true;
      break;
    } else {
      cnt += j->num_nodagger();
    }
  }
  return found;
} 
