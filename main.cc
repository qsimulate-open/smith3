//
// author : Toru Shiozaki
// date   : Jan 2012
//

// This program is supposed to perform Wick's theorem for multireference problems.
// Spin averaged quantities assumed.

#include <iostream>
#include <list>
#include "twoop.h"
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

  return 0;
}

