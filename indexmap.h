//
// Author : Toru Shiozaki
// Date   : Feb 2009
//

#ifndef _smith_indexmap_h
#define _smith_indexmap_h

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <list>
#include <stdexcept>

// This defines the index classes. If you want to generalize this generator
// to more general cases (RASPT2, for instance), then just add some entry.
// Indices will be sorted using these numbers when tensors are canonicalized.
class IndexMap {
  protected:
    std::list<std::pair<std::string, std::pair<int,int> > > map_;
  public:
    IndexMap() { 
      map_.push_back(std::make_pair("c", std::make_pair(0, 28)));
      map_.push_back(std::make_pair("x", std::make_pair(1, 6)));
      map_.push_back(std::make_pair("a", std::make_pair(2, 232)));
    };
    ~IndexMap() {};
    int num_orb_class() const { return map_.size(); };
    int size() const { return num_orb_class(); };

    const int type(const std::string& type_) const {
      auto iter = map_.begin();
      for (; iter != map_.end(); ++iter) if (iter->first == type_) break;
      if (iter == map_.end()) throw std::runtime_error("key is no valid in Index::type()"); 
      return iter->second.first;
    };
    std::list<std::pair<std::string, std::pair<int,int> > >::const_iterator begin() const { return map_.begin(); };
    std::list<std::pair<std::string, std::pair<int,int> > >::const_iterator end() const { return map_.end(); };
};

#endif

