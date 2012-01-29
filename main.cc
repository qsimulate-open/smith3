//
// author : Toru Shiozaki
// date   : Jan 2012
//

// This program is supposed to perform Wick's theorem for multireference problems.
// Spin averaged quantities assumed.

#include <iostream>
#include <list>
#include "equation.h"

using namespace std;

int main() {

  shared_ptr<Op> proj(new Op("proj", "c", "c", "x", "x"));
  shared_ptr<Op> f(new Op("f", "g", "g"));
  shared_ptr<Op> T(new Op("T", "x", "x", "c", "c"));

  list<shared_ptr<Op> > d;
  d.push_back(proj);
  d.push_back(f);
  d.push_back(T);

  shared_ptr<Diagram> di(new Diagram(d));
  shared_ptr<Equation> eq(new Equation(di));
  eq->factorize();
  eq->active();
  eq->print();

}
