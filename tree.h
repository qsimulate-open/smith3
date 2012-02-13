//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __TREE_H
#define __TREE_H

#include <memory>
#include "equation.h"
#include "listtensor.h"

class Tree;

class BinaryContraction : public Tensor {
  protected:
    // *this = tensor_ * subtrees
    std::shared_ptr<Tensor> tensor_;
    std::list<std::shared_ptr<Tree> > subtree_;

  public:
    BinaryContraction(std::shared_ptr<Tensor> a, std::list<std::shared_ptr<Tree> > b) : Tensor(), tensor_(a), subtree_(b) {};
    ~BinaryContraction() {};
};

// Tree is a tensor
class Tree {
  protected:
    // recursive structure. BinaryContraction contains Tree.
    std::list<std::shared_ptr<Tensor> > tensors_;

  public:
    Tree(const std::shared_ptr<Equation> eq);
    ~Tree() {};

    bool done() const;
    void print() const;

    std::string generate() const;

};


#endif
