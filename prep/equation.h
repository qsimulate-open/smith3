#ifndef __EQUATION_H
#define __EQUATION_H

#include "diagram.h"
#include <sstream>
#include <memory>
#include <list>
#include <vector>
#include <string>

namespace SMITH3 {
namespace Prep {

class Equation {
  protected:
    std::list<std::shared_ptr<Diagram>> diagram_;
    std::string label_;
    double fac_;
    std::string tree_type_;

  public:
    Equation(const std::string l, const std::initializer_list<std::vector<std::shared_ptr<Tensor>> > in, const std::string scalar = "", const double d = 1.0) : label_(l), fac_(d), tree_type_("") {
      std::list<int> max;
      for (auto& i : in) max.push_back(i.size());

      std::list<std::list<std::shared_ptr<Tensor>>> out;
      std::list<int> current(in.size(), 0);
      std::list<int> start = current;
      do {
        // set the current vector 
        std::list<std::shared_ptr<Tensor>> cc;
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
        if (d == 1.0) {
          diagram_.push_back(std::shared_ptr<Diagram>(new Diagram(i, ss.str(), scalar)));
        } else {
          diagram_.push_back(std::shared_ptr<Diagram>(new Diagram(i, ss.str(), d)));
        }
        ++cnt;
      }

    };


    std::string tree_label() const { return "t" + label_; }

    void set_tree_type(std::string ttype) { tree_type_=ttype; };

    void merge(std::shared_ptr<Equation> o) {
      diagram_.insert(diagram_.end(), o->diagram_.begin(), o->diagram_.end());
    }

    std::string generate(std::initializer_list<std::shared_ptr<Equation>> o) const {
      std::stringstream ss;
      for (auto& i : diagram_) ss << i->construct_str();
      for (auto& i : diagram_) ss << i->diagram_str();
      for (auto& i : diagram_) ss << i->equation_str();

      for (auto i = diagram_.begin(); i != diagram_.end(); ++i)
        if (i != diagram_.begin())
          ss << "  " << diagram_.front()->eqn_label() << "->merge(" << (*i)->eqn_label() << ");" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->duplicates();" << std::endl;
      ss << "  " << diagram_.front()->eqn_label() << "->active();" << std::endl;

#if 1 // prune density matrix equation to keep following terms, for testing
      if (!tree_type_.empty() && tree_type_ == "density") {
        //ss << "  " << "list<string> terms = {\"x\"};" << std::endl;
        //ss << "  " << "list<string> terms = {\"c\"};" << std::endl;
          ss << "  " << "list<string> terms = {\"c\", \"x\"};" << std::endl;
        //ss << "  " << "list<string> terms = {\"c\", \"a\"};" << std::endl;
        ss << "  " <<  diagram_.front()->eqn_label() << "->term_select(terms);" << std::endl;
      }
#endif

      if (!tree_type_.empty() && tree_type_ == "residual") {
          ss << "  shared_ptr<Tree> " << tree_label() << "(new Residual(e" << diagram_.front()->label() << ", \"" << tree_type_ << "\"));" << std::endl;
      } else if (!tree_type_.empty() && tree_type_ == "energy") {
          ss << "  shared_ptr<Tree> " << tree_label() << "(new Energy(e" << diagram_.front()->label() << ", \"" << tree_type_ << "\"));" << std::endl;
      } else if (!tree_type_.empty() && tree_type_ == "correction") {
          ss << "  shared_ptr<Tree> " << tree_label() << "(new Correction(e" << diagram_.front()->label() << ", \"" << tree_type_ << "\"));" << std::endl;
      } else if (!tree_type_.empty() && tree_type_ == "density") {
          ss << "  shared_ptr<Tree> " << tree_label() << "(new Density(e" << diagram_.front()->label() << ", \"" << tree_type_ << "\"));" << std::endl;
      } else {
          throw std::logic_error("prep/equation.cc error, tree must be of derived type");
          // ss << "  shared_ptr<Tree> " << tree_label() << "(new Tree(e" << diagram_.front()->label() << "));" << std::endl;
      }

      if (o.size() == 0) {
        ss << "  " << tree_label() << "->sort_gamma();" << std::endl;
      } else {
        for (auto i = o.begin(); i != o.end(); ++i)
          ss << "  " << tree_label() << "->sort_gamma(" << (*i)->tree_label() << "->gamma());" << std::endl;
      }
      ss << std::endl;
      return ss.str();
    };


};

}}


#endif
