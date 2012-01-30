//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __LISTTENSOR_H
#define __LISTTENSOR_H

#include "op.h"
#include "active.h"
#include "diagram.h"
#include <memory>
#include <list>

// class for a vector of tensors
class ListTensor {
  protected:
    double fac_;
    std::list<std::shared_ptr<Op> > list_;
    std::shared_ptr<Active> active_;

    bool dagger_;

  public:
    ListTensor(std::shared_ptr<Diagram> d);
    ListTensor(double f, std::list<std::shared_ptr<Op> > ve, std::shared_ptr<Active> ac, bool d)
      : fac_(f), list_(ve), active_(ac), dagger_(d) {};
    ~ListTensor() {};

    void print() const;

    int length() const { return list_.size(); };
};


#endif
