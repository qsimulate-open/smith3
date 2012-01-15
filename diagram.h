//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "op.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <list>
#include <map>

class Diagram {
  protected:
    std::list<Op> op_;
    double fac_;

  public:
    Diagram(std::list<Op> op) : op_(op), fac_(1.0) { };
    ~Diagram() {};

    double fac() const { return fac_; };

    void print() const {
      std::cout << std::setw(4) << std::setprecision(1) << std::fixed <<  fac_ << " ";
      std::map<std::shared_ptr<std::string>, int> dict;
      std::map<std::shared_ptr<std::string>, int> spin;
      for (auto i = op_.begin(); i != op_.end(); ++i)
        i->print(dict, spin);
      std::cout << std::endl;
    };

    void reduce_one(const int skip) {
      // find the first dagger operator in list<Op>
      auto i = op_.begin();
      std::pair<std::shared_ptr<std::string>*, std::shared_ptr<std::string>* > data; // safe because they are held by tensors
      for (; i != op_.end(); ++i) {
        // this simultaneously eliminates one entry. op_ is therefore modified here
        data = i->first_dagger_noactive();
        if (data.first) break;
      }
      if (!data.first) return;

      // skip until it comes to skip
      int cnt = 0;
      for (auto j = op_.begin(); j != op_.end(); ++j) {
        // cannot contract with self
        if (i == j) continue;
        // all possible contraction pattern taken for *j (returned as a list).
        if (cnt + j->num_nodagger() > skip) {
          j->contract(data, skip-cnt);
          break;
        } else {
          cnt += j->num_nodagger();
        }
      }

      this->print();
    }; 

};

#endif
