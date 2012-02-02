//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __TREE_H
#define __TREE_H

#include <memory>
#include "equation.h"
#include "listtensor.h"

class Tree {
  protected:
    // recursive structure
    std::list<std::shared_ptr<Tree> > subtree_;
    std::list<std::shared_ptr<ListTensor> > tensor_;

  public:
    Tree(const std::shared_ptr<Equation> eq);
    ~Tree() {};

    bool done() const;
    void print() const;

};


#endif
