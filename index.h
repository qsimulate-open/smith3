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

  public:
    Index(std::string lab) : label_(lab), num_(0) {};
    ~Index() {};

    const int num() const { return num_; };
    void set_num(const int i) { num_ = i; };
    const std::string label() const { return label_; };

    std::string str() const {
      std::stringstream ss;
      ss << label_ << num_; 
      return ss.str();
    };

};


class Spin {
  protected:
    int num_;
  public:
    Spin() : num_(0) {};
    ~Spin() {};

    void set_num(const int i) { num_ = i; };

    std::string str() const {
      std::stringstream ss;
      ss << "(" << num_ << ")";
      return ss.str();
    };
};

#endif
