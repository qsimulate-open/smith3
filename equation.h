//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __EQUATION_H
#define __EQUATION_H

#include <list>
#include "diagram.h"

class Equation {
  protected:
    std::list<std::shared_ptr<Diagram> > diagram_;

  public:
    Equation(std::shared_ptr<Diagram>);
    ~Equation() {};

    void print();
    void active();
    void factorize();

};

#endif

