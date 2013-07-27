//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: equation.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the SMITH3 package.
//
// The SMITH3 package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SMITH3 package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SMITH3 package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include "equation.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace smith;

Equation::Equation(shared_ptr<Diagram> in, std::string nam) : name_(nam) {

  list<shared_ptr<Diagram>> out = in->get_all();

  if (out.size() != 0) {
    while (out.front()->num_dagger()) {
      list<shared_ptr<Diagram>> out2;
      for (auto& j : out) {
        for (int i = 0; i != j->num_dagger(); ++i) {
          shared_ptr<Diagram> n = j->copy();
          bool found = n->reduce_one_noactive(i);
          if (!found) break;
          if (n->valid() || n->done()) {
            out2.push_back(n);
            if (n->done_noactive()) diagram_.push_back(n);
          }
        }
      }
      out = out2;
      if (out.size() == 0) break;
    }
    // collect target indices from excitation operators.
    for (auto& i : diagram_) i->refresh_indices();
  }
}


bool Equation::targets() const {
  bool out = false;
  for (auto& i : diagram_) {
    if (i->has_target_index()) out = true;
  }
  return out;
}


void Equation::term_select(list<string> t) {
  // go through diagrams and if do not contain desired target indices (t), remove.
  list<list<shared_ptr<Diagram>>::iterator> rm;
  bool keep = true;
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) {
    const list<shared_ptr<Operator>> ops = (*i)->op();
    for (auto& j : ops) {
      // find excitation operator
      if (j->label().empty()) {
        // compare op index label
        list<tuple<shared_ptr<Index>*,int, int>> q_ops = j->op();
        for (auto& k : q_ops) {
          bool found = false;
          for (auto& term : t)
            if ((*get<0>(k))->label() == term) found = true;
          if (!found) {
            rm.push_back(i);
            break;
          }
        }
      }
    }
  }
  for (auto& it : rm) diagram_.erase(it);
}


// print. This triggers Diagram::refresh_indices().
void Equation::print() {
  for (auto& i : diagram_) i->print();
}

// processes active part
void Equation::active() {
  for (auto& i : diagram_) i->active();
}

void Equation::refresh_indices() {
  for (auto& i : diagram_) i->refresh_indices();
}

// find identical terms
void Equation::duplicates() {
  duplicates_(false);
  refresh_indices();
  // TODO this is only valid for projection up to doubles
  // For any-order algorithm, we need to use a generic algorithm.
  duplicates_(true);
}


void Equation::absorb_ket() {
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) {
    // function designed to be used after contractions
    assert((*i)->done());
    if ((*i)->braket().second) {
      // if diagram does not have rdms
      if ((*i)->active_indices().empty()) {
          (*i)->set_ket(false);
          (*i)->set_bra(true);
      } else {
        // todo redo indices for diagrams with rdms
      }
    }
  }
}


void Equation::duplicates_(const bool proj) {
  list<list<shared_ptr<Diagram>>::iterator> rm;
  for (auto i = diagram_.begin(); i != diagram_.end(); ++i) {
    bool found = false;
    // all possible permutations generated here
    do {
      // find identical
      auto j = i;
      for (++j ; j != diagram_.end(); ++j) {
        if ((*i)->identical(*j)) {
          found = true;
          if (!proj) {
            (*j)->fac() += (*i)->fac();
            rm.push_back(i);
            if ((*j)->fac() == 0) throw logic_error("I don't think that this happens. Check! Equation::factorize1_");
          } else {
            (*j)->add_dagger();
            rm.push_back(i);
            if ((*j)->fac() != (*i)->fac()) throw logic_error("I don't think that this happens. Check! Equation::factorize2_");
          }
        }
      }
      if (found) break;
    } while ((*i)->permute(proj));
  }
  for (auto& it : rm) diagram_.erase(it);
}


