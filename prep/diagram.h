#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "tensor.h"
#include <initializer_list>

class Diagram {
  protected:
    std::list<std::shared_ptr<Tensor> > op_;
    std::string label_;
    std::string scalar_;

  public:
    Diagram(const std::list<std::shared_ptr<Tensor> > o, const std::string la, std::string s) : op_(o), label_(la), scalar_(s) {}; 

    const std::list<std::shared_ptr<Tensor> >& op() const { return op_; }; 
    std::string label() const { return label_; };
    std::string diag_label() const { return "d" + label_; };
    std::string eqn_label() const { return "e" + label_; };

    std::string index_str() const {
      std::stringstream ss;
      for (auto j = op_.begin(); j != op_.end(); ++j) {
        if (j != op_.begin()) ss << ", ";
        ss << (*j)->tag();
      }
      return ss.str();
    }

    std::string construct_str() const {
      std::stringstream ss;
      ss << "  list<shared_ptr<Op> > " << label_ << " = {" << index_str() << "};" << std::endl;
      return ss.str();
    }

    std::string scalar() const { return scalar_; };

    std::string diagram_str() const {
      std::stringstream ss;
      ss << "  shared_ptr<Diagram> " << diag_label() << "(new Diagram(" << label() << (scalar().empty() ? "" : ", \""+scalar()+"\"") << "));" << std::endl;
      return ss.str();
    };

    std::string equation_str() const {
      std::stringstream ss;
      ss << "  shared_ptr<Equation> " << eqn_label() << "(new Equation(" << diag_label() << ", theory));" << std::endl;
      return ss.str();
    }

};

#endif
