//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#include "equation.h"
#include <cassert>

using namespace std;

Equation::Equation(shared_ptr<Diagram> in) {

  list<shared_ptr<Diagram> > out = in->get_all();
  assert(out.size() != 0);

  while (out.front()->num_dagger()) {
    list<shared_ptr<Diagram> > out2;
    for (auto j = out.begin(); j != out.end(); ++j) {
      for (int i = 0; i != (*j)->num_dagger(); ++i) {
        shared_ptr<Diagram> n = (*j)->copy();
        bool done = n->reduce_one_noactive(i);
        if (!done) break;
        if (n->valid() || n->done()) {
          out2.push_back(n);
          if (n->done_noactive()) diagram_.push_back(n);
        }
      }
    }
    out = out2;
  }
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) (*i)->refresh_indices();

}


// print. This triggers Diagram::refresh_indices().
void Equation::print() {
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) (*i)->print();
}

// processes active part
void Equation::active() {
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) (*i)->active();
}

// find identical terms
void Equation::factorize() {
  list<list<shared_ptr<Diagram> >::iterator> rm;
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) {
    bool found = false;
    // all possible permutations generated here
    do {
      // find identical
      auto j = i;
      for (++j ; j != diagram_.end(); ++j) {
        if ((*i)->identical(*j)) {
          found = true;
          (*j)->fac() += (*i)->fac();
          rm.push_back(i);
          if ((*j)->fac() == 0) throw logic_error("I don't think that this happens. Check!");
        }
      }
      if (found) break;
    } while ((*i)->permute());
  }
  for (auto iter = rm.begin(); iter != rm.end(); ++iter) diagram_.erase(*iter);
}

