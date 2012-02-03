//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#ifndef __TENSOR_H
#define __TENSOR_H

// In order to treat Active and Op on the same footing.
// To make it easy, this is again driven by Index::num_. 
// Perhaps it's better if it is done by pointer itself, but for the time being it's fine...
//
// Index::spin_ is not set. 

// this object should generate some codes.

#include <memory>
#include <string>
#include <list>
#include "op.h"
#include "active.h"

class Tensor {
  protected:
    // factor
    double factor_;
    // label of this tensor.
    std::string label_;
    // a list of indices
    std::list<std::shared_ptr<Index> > index_;

    // if this tensor is active, it has an internal structure
    std::shared_ptr<Active> active_;
    std::shared_ptr<Tensor> merged_;

  public:
    Tensor(const std::shared_ptr<Op> op);
    Tensor(const std::shared_ptr<Active> active);
    ~Tensor() {};

    std::list<std::shared_ptr<Index> >& index() { return index_; };
    const std::list<std::shared_ptr<Index> >& index() const { return index_; };

    std::string str() const;
    void set_factor(const double a) { factor_ = a; };

    std::shared_ptr<Active> active() { return active_; };

    bool all_active() const;

    void merge(std::shared_ptr<Tensor> o);

    std::string generate() const;
};

#endif

