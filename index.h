//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __INDEX_H
#define __INDEX_H

#include <string>
#include <sstream>
#include <cassert>

class Spin {
  protected:
    int num_;
  public:
    Spin() : num_(0) {};
    ~Spin() {};

    int num() const { return num_; };
    void set_num(const int i) { num_ = i; };

    std::string str() const {
      std::stringstream ss;
      ss << "(" << num_ << ")";
      return ss.str();
    };
};


class Index {
  protected:
    std::string label_;
    int num_;
    bool dagger_;
    std::shared_ptr<Spin> spin_;

  public:
    Index(std::string lab, bool dag) : label_(lab), num_(0), dagger_(dag) {};
    ~Index() {};

    int num() const { return num_; };
    bool dagger() const { return dagger_; };
    void set_num(const int i) { num_ = i; };
    const std::string label() const { return label_; };
    void set_label(const std::string& a) { label_ = a; };

    void set_spin(const std::shared_ptr<Spin> s) { spin_ = s; };
    std::shared_ptr<Spin> spin() { assert(spin_); return spin_; };
    const std::shared_ptr<Spin> spin() const { assert(spin_); return spin_; };

    bool same_spin(const std::shared_ptr<Index>& o) const { return o->spin() == spin(); };

    std::string str(const bool& opr = true) const {
      std::stringstream ss;
      ss << label_ << num_;
      if (dagger_ && opr) ss << "+";
      if (opr) {
        assert(spin_);
        ss << spin_->str();
      }
      return ss.str();
    };

    void print() const { std::cout << str() << std::endl; };
};


#endif
