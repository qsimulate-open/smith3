#ifndef __EQUATION_H
#define __EQUATION_H

#include "diagram.h"
#include <sstream>
#include <memory>
#include <list>
#include <vector>

class Equation {
  protected:
    std::list<std::shared_ptr<Diagram> > diagram_;
    std::string label_;
    std::string scalar_;

  public:
    Equation(const std::string l, const std::initializer_list<std::vector<std::shared_ptr<Tensor> > > in, const std::string scalar = "") : label_(l), scalar_(scalar) {
      std::list<int> max;
      for (auto& i : in) max.push_back(i.size());

      std::list<std::list<std::shared_ptr<Tensor> > > out;
      std::list<int> current(in.size(), 0);
      std::list<int> start = current;
      do {
        // set the current vector 
        std::list<std::shared_ptr<Tensor> > cc;
        auto inp = in.begin();
        for (auto i = current.begin(); i != current.end(); ++i, ++inp) cc.push_back((*inp)[*i]);
        out.push_back(cc);

        // get the next vector
        auto m = max.rbegin();
        for (auto i = current.rbegin(); i != current.rend(); ++i, ++m) {
          if (++*i == *m) {
            *i = 0; 
          } else {
            break;
          } 
        }
      } while (current != start);

      // construct Diagrams
      int cnt = 0;
      for (auto& i : out) {
        std::stringstream ss; ss << label_ << cnt; 
        diagram_.push_back(std::shared_ptr<Diagram>(new Diagram(i, ss.str())));
        ++cnt;
      }

    };


    std::string tree_label() const { return "t" + label_; }


    void merge(std::shared_ptr<Equation> o) {
      diagram_.insert(diagram_.end(), o->diagram_.begin(), o->diagram_.end());
    }

    std::string generate(std::initializer_list<std::shared_ptr<Equation> > o) const {
      std::stringstream ss;
      for (auto& i : diagram_)
        ss << i->construct_str();
      for (auto& i : diagram_)
        ss << "  shared_ptr<Diagram> " << i->diag_label() << "(new Diagram(" << i->label() << (scalar_.empty() ? "" : ", "+scalar_) << "));" << std::endl;
      for (auto& i : diagram_)
        ss << "  shared_ptr<Equation> " << i->eqn_label() << "(new Equation(" << i->diag_label() << ", theory));" << std::endl;
      for (auto i = diagram_.begin(); i != diagram_.end(); ++i)
        if (i != diagram_.begin())
          ss << "  " << diagram_.front()->eqn_label() << "->merge(" << (*i)->eqn_label() << ");" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->duplicates();" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->active();" << std::endl;
      ss << "  shared_ptr<Tree> " << tree_label() << "(new Tree(e" << diagram_.front()->label() << "));" << std::endl;
      for (auto i = o.begin(); i != o.end(); ++i) {
        ss << "  " << tree_label() << "->sort_gamma(" << (*i)->tree_label() << "->gamma());" << std::endl;
      }
      return ss.str();
    };


};


#endif