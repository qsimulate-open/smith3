//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: indexmap.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef _smith_indexmap_h
#define _smith_indexmap_h

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <list>
#include <stdexcept>

namespace smith {

// This defines the index classes. If you want to generalize this generator
// to more general cases (RASPT2, for instance), then just add some entry.
// Indices will be sorted using these numbers when tensors are canonicalized.

/// Defines index classes.
class IndexMap {
  protected:
    /// This is list of index classes.
    std::list<std::pair<std::string, std::pair<int,int> > > map_;
  public:
    /// Construct index classes.
    IndexMap() { 
      map_.push_back(std::make_pair("c", std::make_pair(0, 28)));
      map_.push_back(std::make_pair("x", std::make_pair(1, 6)));
      map_.push_back(std::make_pair("a", std::make_pair(2, 232)));
    };
    ~IndexMap() {};
    /// Returns map_ size.
    int num_orb_class() const { return map_.size(); };
    /// Also returns map_ size.
    int size() const { return num_orb_class(); };

    /// Returns class type based on map_.
    const int type(const std::string& type_) const {
      auto iter = map_.begin();
      for (; iter != map_.end(); ++iter) if (iter->first == type_) break;
      if (iter == map_.end()) throw std::runtime_error("key is no valid in Index::type()"); 
      return iter->second.first;
    };
    /// Returns index class beginning iterator.
    std::list<std::pair<std::string, std::pair<int,int> > >::const_iterator begin() const { return map_.begin(); };
    /// Returns index class end iterator.
    std::list<std::pair<std::string, std::pair<int,int> > >::const_iterator end() const { return map_.end(); };
};

}

#endif

