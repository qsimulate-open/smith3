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

  shared_ptr<Op> proj(new Op("proj", "c", "x", "a", "a"));
  shared_ptr<Op> f(new Op("f", "g", "g"));
  shared_ptr<Op> T(new Op("T", "a", "a", "c", "x"));

  list<shared_ptr<Op> > d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  shared_ptr<Diagram> di(new Diagram(d));
  di->print();

  list<shared_ptr<Diagram> > out;
  out.push_back(di);

int cnt = 0;
  while (out.front()->num_dagger()) {
    list<shared_ptr<Diagram> > out2;
    for (auto j = out.begin(); j != out.end(); ++j) {
      for (int i = 0; i != (*j)->num_dagger(); ++i) {
        shared_ptr<Diagram> n = (*j)->copy();
        bool done = n->reduce_one_noactive(i);
        if (!done) break;
        if (n->valid() || n->done()) {
          out2.push_back(n);
        }
      }
    }
    out = out2;
if (cnt == 5) break;
++cnt;
for (auto iter = out.begin(); iter != out.end(); ++iter) (*iter)->print();
cout << "----------------" << endl;
  }

  for (auto iter = out.begin(); iter != out.end(); ++iter) (*iter)->print();

}
