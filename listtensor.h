//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __LISTTENSOR_H
#define __LISTTENSOR_H

#include "tensor.h"
#include "diagram.h"
#include <memory>
#include <list>

// class for a vector of tensors
class ListTensor {
  protected:
    double fac_;
    std::list<std::shared_ptr<Tensor> > list_;

    bool dagger_;

  public:
    ListTensor(std::shared_ptr<Diagram> d);
    ListTensor(double f, std::list<std::shared_ptr<Tensor> > ve, bool d)
      : fac_(f), list_(ve), dagger_(d) {};
    ~ListTensor() {};

    void print() const;
    void absorb_all_internal();

    int length() const { return list_.size(); };
    std::shared_ptr<Tensor> front() const { return list_.front(); };
    std::shared_ptr<ListTensor> rest() const ;

    bool dagger() const { return dagger_; };

    std::string generate() const;
};


#endif
