//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "op.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <map>

class Diagram {
  protected:
    std::list<std::shared_ptr<Op> > op_;
    double fac_;

  public:
    Diagram(std::list<std::shared_ptr<Op> > op) : op_(op), fac_(1.0) { };
    Diagram() : fac_(1.0) { };
    // copy constructor is complicated but preserves the same topology as this.
    ~Diagram() {};

    std::shared_ptr<Diagram> copy() const;

    double fac() const { return fac_; };

    const std::list<std::shared_ptr<Op> >& op() const { return op_; };
    void set_op(const std::list<std::shared_ptr<Op> >& o) { op_ = o; };

    // refresh the indices
    void refresh_indices();

    // printing function
    // CAUTION: it also refreshes the indices
    void print();

    int num_dagger() const {
      int out = 0;
      for (auto i = op_.begin(); i !=  op_.end(); ++i) out += (*i)->num_dagger();
      return out;
    };

    // this function performs one contraction ** IN PLACE **
    bool reduce_one_noactive(const int skip);

};

#endif
