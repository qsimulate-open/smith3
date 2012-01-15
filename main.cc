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

using namespace std;

int main() {

  shared_ptr<Op> proj(new Op("proj", "c", "c", "a", "a"));
  shared_ptr<Op> f(new Op("f", "g", "g"));
  shared_ptr<Op> T(new Op("T", "a", "a", "c", "c"));

  list<shared_ptr<Op> > d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  shared_ptr<Diagram> di(new Diagram(d));
  di->print();

  cout << "first generation" << endl;

  list<shared_ptr<Diagram> > out;
  for (int i = 0; i != di->num_dagger(); ++i) {
    shared_ptr<Diagram> n = di->copy();
    bool done = n->reduce_one_noactive(i);
    if (!done) break;
    n->print();
    out.push_back(n);
  }

  cout << "second generation" << endl;

  for (auto j = out.begin(); j != out.end(); ++j) {
    list<shared_ptr<Diagram> > out2;
    for (int i = 0; i != (*j)->num_dagger(); ++i) {
      shared_ptr<Diagram> n = (*j)->copy();
      bool done = n->reduce_one_noactive(i);
      out2.push_back(n);
      if (!done) break;
      n->print();
    }
break;
  }

  return 0;
}

