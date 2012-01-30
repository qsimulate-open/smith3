//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "active.h"
#include "op.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <map>

class Diagram {
  protected:
    // a list of operators
    std::list<std::shared_ptr<Op> > op_;
    // a constant factor
    double fac_;
    // active part
    std::shared_ptr<Active> rdm_;

    // if this Diagram has a daggered counterpart (often the case for residual equations).
    bool dagger_;

  public:
    Diagram(std::list<std::shared_ptr<Op> > op) : op_(op), fac_(1.0), dagger_(false) { };
    Diagram() : fac_(1.0), dagger_(false) { };
    // copy constructor is complicated but preserves the same topology as this.
    ~Diagram() {};

    // returns a shared_ptr of a diagram that has the same topology as this.
    std::shared_ptr<Diagram> copy() const;

    // generate all combination of diagrams (related to general indices)
    std::list<std::shared_ptr<Diagram> > get_all() const;

    // get functions
    double& fac() { return fac_; };
    const double fac() const { return fac_; };

    // careful returns a const reference
    const std::list<std::shared_ptr<Op> >& op() const { return op_; };
    // set functions for private members
    void set_op(const std::list<std::shared_ptr<Op> >& o) { op_ = o; };
    void set_fac(const double a) { fac_ = a; };

    // refresh the indices
    void refresh_indices();

    // processes the active part
    void active();

    // daggered Diagram added to the sum
    void add_dagger() { dagger_ = true; };

    // permute indices in operators. return false when finished
    bool permute(const bool proj); 
    bool identical(std::shared_ptr<Diagram> o) const;

    // printing function
    // CAUTION: it also refreshes the indices
    void print();
    // this version does not refresh indices
    void print() const;

    // the number of daggered indices
    int num_dagger() const;
    // the number of general indices
    int num_general() const;
    // returns if this diagram has a consistent set of dagger and undaggered indices
    bool consistent_indices() const;

    // this function performs one contraction ** IN PLACE **
    bool reduce_one_noactive(const int skip);

    // returns if this diagram is still valid
    bool valid() const;
    // returns if this diagram is fully contracted and sorted
    bool done() const;
    // returns if this diagram is fully contracted (looking up only nonactive parts) 
    bool done_noactive() const;

    // gathers active indices
    std::list<std::shared_ptr<Index> > active_indices() const; 
};

#endif
