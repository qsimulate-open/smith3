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
    // list of diagrams
    std::list<std::shared_ptr<Diagram> > diagram_;
    // internal function used by factorize()
    void duplicates_(const bool);

    std::string name_;

  public:
    Equation(std::shared_ptr<Diagram>, std::string nam);
    ~Equation() {};

    // print function
    void print();
    // active parts are processed
    void active();
    // identifys the same terms
    void duplicates();
    // refresh indices in each diagram
    void refresh_indices();

    // returns the name of this equation
    std::string name() const { return name_; };

    std::list<std::shared_ptr<Diagram> > diagram() { return diagram_; };

};

#endif

