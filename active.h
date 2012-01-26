//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

//
// implements the active part
//


#ifndef __ACTIVE_H
#define __ACTIVE_H

#include <list>
#include <memory>
#include "op.h"

class RDM {
  protected:
    // prefactor
    double fac_;
    // operators that constitute RDM
    std::list<std::shared_ptr<Index> > index_;
    // kronecker's delta
    std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > > delta_;

  public:
    RDM(const std::list<std::shared_ptr<Index> >& in,
        const std::list<std::pair<std::shared_ptr<Index>, std::shared_ptr<Index> > >& in2,
        const double& f = 1.0)
      : index_(in), delta_(in2), fac_(f) { };
    ~RDM() {};

    void print() const;

    // returns private members
    double factor() const { return fac_; };

    // returns if this is in the final form
    bool done() const;
    bool reduce_done(const std::list<std::shared_ptr<Index> >& done) const;

    // One index is going to be annihilated. done is updated inside the function
    std::list<std::shared_ptr<RDM> > reduce_one(std::list<std::shared_ptr<Index> >& done) const;
};


class Active {
  protected:
    std::list<std::shared_ptr<RDM> > rdm_;
    void reduce(std::shared_ptr<RDM> in);

  public:
    Active(const std::list<std::shared_ptr<Index> >& in);
    ~Active() {};

    void print() const;

};

#endif
