//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: listtensor.h
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


#ifndef __LISTTENSOR_H
#define __LISTTENSOR_H

#include "tensor.h"
#include "diagram.h"
#include <memory>
#include <list>

namespace smith {

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
    std::shared_ptr<Tensor> target() const;

    double fac() const { return fac_; };
    bool dagger() const { return dagger_; };
};

}

#endif
