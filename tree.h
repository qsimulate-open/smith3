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

class BinaryContraction {
  protected:
    // target_ = tensor_ * subtrees
    std::shared_ptr<Tensor> target_;
    std::shared_ptr<Tensor> tensor_;
    std::list<std::shared_ptr<Tree> > subtree_;

  public:
    BinaryContraction(std::shared_ptr<Tensor> o, std::shared_ptr<ListTensor> l);
    ~BinaryContraction() {};

    void print() const;
};

// Tree is a tensor
class Tree {
  protected:
    // recursive structure. BinaryContraction contains Tree.
    std::list<std::shared_ptr<BinaryContraction> > bc_;
    // note that target_ can be NULL (at the very beginning) 
    std::shared_ptr<Tensor> target_;

    std::list<std::shared_ptr<Tensor> > op_;

  public:
    Tree(const std::shared_ptr<Equation> eq);
    Tree(const std::shared_ptr<ListTensor> l);
    ~Tree() {};

    bool done() const;
    void print() const;

    std::string generate() const;

};


#endif
