//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "tensor.h"
#include <sstream>
#include <iomanip>

using namespace std;

Tensor::Tensor(const shared_ptr<Op> op) : factor_(1.0) {
  // label
  label_ = op->label();
  // op
  for (auto i = op->op().begin(); i != op->op().end(); ++i) {
    shared_ptr<Index> in = *get<0>(*i);
    index_.push_back(in);
  }

} 


Tensor::Tensor(const shared_ptr<Active> activ) : factor_(1.0) {
  // label
  label_ = "Gamma";
  // op
  index_ = activ->index();
  active_ = activ;
} 


std::string Tensor::str() const {
  stringstream ss;
  if (factor_ != 1.0) ss << " " << fixed << setw(4) << setprecision(2) << factor_ << " "; 
  ss << label_ << "(";
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    // we don't need the spin part here
    if (i != index_.begin()) ss << ", ";
    ss << (*i)->str(false);
  }
  ss << ")";

  if (merged_) ss << " << " << merged_->str();
  return ss.str();
}


// returns if all the indices are of active orbitals
bool Tensor::all_active() const {
  bool out = true;
  for (auto i = index_.begin(); i != index_.end(); ++i) {
    out &= (*i)->active(); 
  }
  return out;
}


// adds all-active tensor to Active_;
void Tensor::merge(shared_ptr<Tensor> a) {
  assert(active_);  
  merged_ = a; 
  list<list<shared_ptr<Index> >::iterator> remove;
  // remove Index that belongs to a
  for (auto i = a->index().begin(); i != a->index().end(); ++i) {
    const int n = (*i)->num();
    for (auto j = index_.begin(); j != index_.end(); ++j) {
      // careful, this code is driven by numbers
      if ((*j)->num() == n) {
        remove.push_back(j);
        break;
      }
    }
  }
  for (auto i = remove.begin(); i != remove.end(); ++i) {
    index_.erase(*i);
  }
}


bool Tensor::operator==(const Tensor& o) const {
  bool out = true;
  out &= label_ == o.label();
  out &= factor_ == o.factor();
  // TODO for the time being, active is not assumed to be factorized
  if (active() || o.active()) out = false;
  // index
  out &= index_.size() == o.index().size();
  if (out) {
    for (auto i = index_.begin(), j = o.index().begin(); i != index_.end(); ++i, ++j) {
      out &= (*i)->identical(*j);
    } 
  }
  return out;
}


string Tensor::constructor_str(std::string indent) const {
  stringstream ss;
  ss << indent << "std::vector<IndexRange> " << label_ << "_index";
  if (index_.empty()) {
    ss << ";" << endl;
  } else {
    ss << " = vec(";
    for (auto i = index_.begin(); i != index_.end(); ++i)
      ss << (i != index_.begin() ? ", " : "") << (*i)->generate();
    ss << ");" << endl;
  }
  ss << indent << "std::shared_ptr<Tensor<T> > " << label_ << "(new Tensor<T>(" << label_ << "_index, false));";
  return ss.str();
}
