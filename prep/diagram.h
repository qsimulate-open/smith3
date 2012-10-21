#ifndef __DIAGRAM_H
#define __DIAGRAM_H

#include "tensor.h"
#include <initializer_list>

class Diagram {
  protected:
    std::list<std::shared_ptr<Tensor> > op_;
    std::string label_;

  public:
    Diagram(const std::list<std::shared_ptr<Tensor> > o, const std::string la) : op_(o), label_(la) {}; 
#if 0
    Diagram(const std::initializer_list<std::shared_ptr<Tensor> > o) {
      for (auto& i : o) op_.push_back(i);
    }
#endif

    const std::list<std::shared_ptr<Tensor> >& op() const { return op_; }; 
    std::string label() const { return label_; };
    std::string diag_label() const { return "d" + label_; };
    std::string eqn_label() const { return "e" + label_; };

    std::string index_str() {
      std::stringstream ss;
      for (auto j = op_.begin(); j != op_.end(); ++j) {
        if (j != op_.begin()) ss << ", ";
        ss << (*j)->tag();
      }
      return ss.str();
    }

    std::string construct_str() {
      std::stringstream ss;
      ss << "  list<shared_ptr<Op> > " << label_ << " = {" << index_str() << "};" << std::endl;
      return ss.str();
    }

};

#endif
