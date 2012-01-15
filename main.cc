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

  Op proj("proj", "c", "c", "a", "a");
  Op f("f", "g", "g");
  Op T("T", "a", "a", "c", "c");

  list<Op> d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  Diagram di(d);
  di.print();


  cout << "first generation" << endl;

  list<Diagram> out;
  for (int i = 0; i != di.num_dagger(); ++i) {
    Diagram n(di);
    bool done = n.reduce_one_noactive(i);
    if (!done) break;
    n.print();
    out.push_back(n);
  }

  cout << "second generation" << endl;

  for (auto j = out.begin(); j != out.end(); ++j) {
    list<Diagram> out2;
j->print();
    for (int i = 0; i != j->num_dagger(); ++i) {
      Diagram n(*j);
      bool done = n.reduce_one_noactive(i);
      out2.push_back(n);
      if (!done) break;
      n.print();
    }
break;
  }

  return 0;
}

