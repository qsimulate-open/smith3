//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __INDEX_H
#define __INDEX_H

#include <string>
#include <sstream>

class Index {
  protected:
    std::string label_;
    int num_;
    bool dagger_;

  public:
    Index(std::string lab, bool dag) : label_(lab), num_(0), dagger_(dag) {};
    ~Index() {};

    int num() const { return num_; };
    bool dagger() const { return dagger_; };
    void set_num(const int i) { num_ = i; };
    const std::string label() const { return label_; };
    void set_label(const std::string& a) { label_ = a; };

    std::string str(const bool& opr = true) const {
      std::stringstream ss;
      ss << label_ << num_;
      if (dagger_ && opr) ss << "+";
      return ss.str();
    };

};


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

#endif
