//
// Author : Toru Shiozaki
// Date   : Jan 2012
//

#ifndef __TWOOP_H
#define __TWOOP_H

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>

// base class for spin-summed operators

class Op {
  protected:
    // tensor info
    std::string label_;
    std::list<std::tuple<std::shared_ptr<std::string>, bool, std::shared_ptr<std::string> > > tensor_;
    std::list<std::tuple<std::shared_ptr<std::string>, bool, std::shared_ptr<std::string> > > op_;
    
  public:
    Op(const std::string lab) : label_(lab) {};
    ~Op() {};

    std::string label() const { return label_; };

    void print(std::map<std::shared_ptr<std::string>, int>& dict, std::map<std::shared_ptr<std::string>, int>& spin) const {
      // tensor
      std::cout << label_ << "(";
      for (auto i = tensor_.begin(); i != tensor_.end(); ++i) {
        auto iter = dict.find(std::get<0>(*i));
        // if this pointer is not registered...
        if (iter == dict.end()) {
          const int c = dict.size();
          dict.insert(make_pair(std::get<0>(*i), c));
          std::cout << *std::get<0>(*i) << (std::get<1>(*i) ? "+_" : "_") << c << " "; 
        } else {
          std::cout << *std::get<0>(*i) << (std::get<1>(*i) ? "+_" : "_") << iter->second << " ";
        }
      }
      std::cout << ") {";
      for (auto i = op_.begin(); i != op_.end(); ++i) {
        auto iter = dict.find(std::get<0>(*i));
        if (iter == dict.end()) {
          const int c = dict.size();
          dict.insert(make_pair(std::get<0>(*i), c));
          std::cout << *std::get<0>(*i) << (std::get<1>(*i) ? "+_" : "_") << c << " "; 
        } else {
          std::cout << *std::get<0>(*i) << (std::get<1>(*i) ? "+_" : "_") << iter->second << " ";
        }

        auto siter = spin.find(std::get<2>(*i));
        if (siter == spin.end()) {
          const int c = spin.size();
          spin.insert(make_pair(std::get<2>(*i), c));
          std::cout << "(" << c << ") ";
        } else {
          std::cout << "(" << siter->second << ") ";
        }
      }
      std::cout << "}";
    };

};

// spin-summed two-body operator
// T^ab_dc a+_rho b+_sigma c_sigma d_rho

class TwoOp : public Op {

  protected:
    // operator info
    std::shared_ptr<std::string> rho_;
    std::shared_ptr<std::string> sigma_;

    std::shared_ptr<std::string> a_;
    std::shared_ptr<std::string> b_;
    std::shared_ptr<std::string> c_;
    std::shared_ptr<std::string> d_;

    void init() {
      // for convinience
      op_.push_back(std::make_tuple(a_, true,  rho_));
      op_.push_back(std::make_tuple(b_, true,  sigma_));
      op_.push_back(std::make_tuple(c_, false, sigma_));
      op_.push_back(std::make_tuple(d_, false, rho_));
      tensor_ = op_;
    };

  public:
    TwoOp(const std::string lab, const std::string al, const std::string bl, const std::string cl, const std::string dl)
      : Op(lab), a_(new std::string(al)), b_(new std::string(bl)), c_(new std::string(cl)), d_(new std::string(dl)),
        rho_(new std::string("rho")), sigma_(new std::string("sigma")) {
      init();
    };
    TwoOp(const TwoOp& o) : Op(o.label()), a_(new std::string(*o.a())), b_(new std::string(*o.b())), c_(new std::string(*o.c())), d_(new std::string(*o.d())),
        rho_(new std::string(*o.rho())), sigma_(new std::string(*o.sigma())) {
      init();
    };
    ~TwoOp() {};

    std::shared_ptr<std::string> rho() const { return rho_; };
    std::shared_ptr<std::string> sigma() const { return sigma_; };

    std::shared_ptr<std::string> a() const { return a_; };
    std::shared_ptr<std::string> b() const { return b_; };
    std::shared_ptr<std::string> c() const { return c_; };
    std::shared_ptr<std::string> d() const { return d_; };

};


class OneOp : public Op {

  protected:
    // operator info
    std::shared_ptr<std::string> rho_;
    std::shared_ptr<std::string> a_;
    std::shared_ptr<std::string> b_;

    void init() {
      // for convinience
      op_.push_back(std::make_tuple(a_, true,  rho_));
      op_.push_back(std::make_tuple(b_, false, rho_));
      tensor_ = op_;
    };

  public:
    OneOp(const std::string lab, const std::string al, const std::string bl)
      : Op(lab), a_(new std::string(al)), b_(new std::string(bl)), rho_(new std::string("rho")) {
      init();
    };
    OneOp(const OneOp& o) : Op(o.label()), a_(new std::string(*o.a())), b_(new std::string(*o.b())), rho_(new std::string(*o.rho())) {
      init();
    };
    ~OneOp() {};

    std::shared_ptr<std::string> rho() const { return rho_; };
    std::shared_ptr<std::string> a() const { return a_; };
    std::shared_ptr<std::string> b() const { return b_; };
};

#endif
