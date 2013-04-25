#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "tensor.h"
#include <initializer_list>

namespace SMITH3 {
namespace Prep {

class Diagram {
  protected:
    std::list<std::shared_ptr<Tensor>> op_;
    std::string label_;
    std::string scalar_;
    double fac_;

  public:
    Diagram(const std::list<std::shared_ptr<Tensor>> o, const std::string la, std::string s) : op_(o), label_(la), scalar_(s), fac_(1.0) {}; 
    Diagram(const std::list<std::shared_ptr<Tensor>> o, const std::string la, double d) : op_(o), label_(la), fac_(d) {}; 

    const std::list<std::shared_ptr<Tensor>>& op() const { return op_; }; 
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
      ss << "  list<shared_ptr<Operator>> " << label_ << " = {" << index_str() << "};" << std::endl;
      return ss.str();
    }

    std::string scalar() const { return scalar_; };

    std::string diagram_str() const {
      std::stringstream ss;
      if (fac_ == 1.0) {
        ss << "  shared_ptr<Diagram> " << diag_label() << "(new Diagram(" << label() << (scalar().empty() ? "" : ", \""+scalar()+"\"") << "));" << std::endl;
      } else {
        assert(scalar().empty());
        ss << "  shared_ptr<Diagram> " << diag_label() << "(new Diagram(" << label() << ", " << fac_ << "));" << std::endl;
      }
      return ss.str();
    };

    std::string equation_str() const {
      std::stringstream ss;
      ss << "  shared_ptr<Equation> " << eqn_label() << "(new Equation(" << diag_label() << ", theory));" << std::endl;
      return ss.str();
    }

};

}}

#endif
