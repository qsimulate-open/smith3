//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "twoop.h"
#include <iostream>
#include <memory>
#include <list>
#include <map>

class Diagram {
  protected:
    std::list<Op> op_;

  public:
    Diagram(std::list<Op> op) : op_(op) {
    };
    ~Diagram() {};

    void print() const {
      std::map<std::shared_ptr<std::string>, int> dict;
      std::map<std::shared_ptr<std::string>, int> spin;
      for (auto i = op_.begin(); i != op_.end(); ++i)
        i->print(dict, spin);
      std::cout << std::endl;
    };

};

#endif
