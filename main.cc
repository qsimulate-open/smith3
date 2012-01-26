//
// author : Toru Shiozaki
// date   : Jan 2012
//

// This program is supposed to perform Wick's theorem for multireference problems.
// Spin averaged quantities assumed.

#include <iostream>
#include <list>
#include "op.h"
#include "diagram.h"
#include "active.h"

using namespace std;

int main() {

  shared_ptr<Op> proj(new Op("proj", "c", "x", "a", "x"));
  shared_ptr<Op> f(new Op("f", "g", "g"));
  shared_ptr<Op> T(new Op("T", "a", "x", "c", "x"));

  list<shared_ptr<Op> > d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  shared_ptr<Diagram> di(new Diagram(d));
  list<shared_ptr<Diagram> > out = di->get_all();

  list<shared_ptr<Diagram> > final;

  while (out.front()->num_dagger()) {
    list<shared_ptr<Diagram> > out2;
    for (auto j = out.begin(); j != out.end(); ++j) {
      for (int i = 0; i != (*j)->num_dagger(); ++i) {
        shared_ptr<Diagram> n = (*j)->copy();
        bool done = n->reduce_one_noactive(i);
        if (!done) break;
        if (n->valid() || n->done()) {
          out2.push_back(n);
          if (n->done_noactive()) final.push_back(n);
        }
      }
    }
    out = out2;
  }

  for (auto iter = final.begin(); iter != final.end(); ++iter) (*iter)->print();

  for (auto iter = final.begin(); iter != final.end(); ++iter) {
    shared_ptr<Active> ac(new Active((*iter)->active_indices()));
    ac->print();
  }

}
