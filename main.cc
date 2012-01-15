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

  TwoOp proj("proj", "c", "c", "a", "a");
  OneOp f("f", "g", "g");
  TwoOp T("T", "a", "a", "c", "c");

  list<Op> d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  Diagram di(d);
  di.print();

  for (int i = 0; i != 10; ++i) {
    Diagram n(di);
    cout << "i = " << i << endl;
    for (int j = 0; j != 10; ++i) n.reduce_one(j);
  }

  return 0;
}

