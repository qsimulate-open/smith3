//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: active.h
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


//
// implements the active part
//


#ifndef __ACTIVE_H
#define __ACTIVE_H

#include <string>
#include <list>
#include <memory>
#include "op.h"
#include "rdm.h"

class Active {
  protected:
    std::list<std::shared_ptr<RDM> > rdm_;
    void reduce(std::shared_ptr<RDM> in);

    mutable int count__;

  public:
    Active(const std::list<std::shared_ptr<Index> >& in);
    ~Active() {};

    void print(const std::string& indent = "") const;
    const std::list<std::shared_ptr<Index> > index() const;

    std::string generate(const std::string cindent, const std::string lab, const std::list<std::shared_ptr<Index> >& loop) const;
    std::vector<int> required_rdm() const;
};

#endif
