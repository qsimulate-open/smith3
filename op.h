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
#include <stdexcept>

// base class for spin-summed operators

class Op {
  protected:
    // tensor info
    std::string label_;
    std::list<std::tuple<std::shared_ptr<std::string>*, int, std::shared_ptr<std::string>* > > op_;
    
  public:
    Op(const std::string lab) : label_(lab) {};
    ~Op() {};

    std::string label() const { return label_; };

    int num_nodagger() const {
      int out = 0;
      for (auto i = op_.begin(); i != op_.end(); ++i) if (std::get<1>(*i)==1) ++out;
      return out;
    };
    int num_dagger() const {
      int out = 0;
      for (auto i = op_.begin(); i != op_.end(); ++i) if (std::get<1>(*i)==0) ++out;
      return out;
    };

    // caution, this changes op_ and returns daggered index
    std::pair<std::shared_ptr<std::string>*, std::shared_ptr<std::string>* > first_dagger_noactive() {
      std::pair<std::shared_ptr<std::string>*, std::shared_ptr<std::string>* > out;
      auto i = op_.begin(); 
      for (; i != op_.end(); ++i) {
        if (std::get<1>(*i)==0 && **std::get<0>(*i) != "x") { // "x" is active orbitals
          out = std::make_pair(std::get<0>(*i), std::get<2>(*i));
          break;
        }
      }
      if (out.first) std::get<1>(*i) = -1;
      return out;
    };

    // returns if you can contract two labels
    bool contractable(std::string a, std::string b) { return a == b || a == "g" || b == "g"; };

    std::shared_ptr<std::string>* survive(std::shared_ptr<std::string>* a, std::shared_ptr<std::string>* b) {
      if (**a == **b) return a;
      else if (**a == "g" && **b != "g") return b;
      else if (**a != "g" && **b == "g") return a;
      else throw std::runtime_error("strange in survive");
    };

    double contract(std::pair<std::shared_ptr<std::string>*, std::shared_ptr<std::string>* >& dat, const int skip) {
      int cnt = 0;
      auto i = op_.begin();
      double fac = 0.0;
      for (; i != op_.end(); ++i) {
        if (std::get<1>(*i)!=1) continue;
        if (contractable(**std::get<0>(*i), **dat.first)) {
          if (cnt == skip) {
            *std::get<0>(*i) = *survive(std::get<0>(*i), dat.first); 
            *dat.first = *std::get<0>(*i);
            if (*dat.second == *std::get<2>(*i)) {
              fac = 2.0;
            } else {
              fac = 1.0;
              *dat.second = *std::get<2>(*i);
            }
            break;
          } else {
            ++cnt;
          }
        }
      }
      std::get<1>(*i) = -1;
      return fac;
    }; 

    void print(std::map<std::shared_ptr<std::string>, int>& dict, std::map<std::shared_ptr<std::string>, int>& spin) const {
      // tensor
      std::cout << label_ << "(";
      for (auto i = op_.begin(); i != op_.end(); ++i) {
        auto iter = dict.find(*std::get<0>(*i));
        // if this pointer is not registered...
        if (iter == dict.end()) {
          const int c = dict.size();
          dict.insert(make_pair(*std::get<0>(*i), c));
          std::cout << **std::get<0>(*i) << c << " "; 
        } else {
          std::cout << **std::get<0>(*i) << iter->second << " ";
        }
      }
      std::cout << ") {";
      for (auto i = op_.begin(); i != op_.end(); ++i) {
        if (std::get<1>(*i) < 0) continue;
        auto iter = dict.find(*std::get<0>(*i));
        if (iter == dict.end()) {
          const int c = dict.size();
          dict.insert(make_pair(*std::get<0>(*i), c));
          std::cout << **std::get<0>(*i) << (std::get<1>(*i)==0 ? "+_" : "_") << c << " "; 
        } else {
          std::cout << **std::get<0>(*i) << (std::get<1>(*i)==0 ? "+_" : "_") << iter->second << " ";
        }

        auto siter = spin.find(*std::get<2>(*i));
        if (siter == spin.end()) {
          const int c = spin.size();
          spin.insert(make_pair(*std::get<2>(*i), c));
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
      op_.push_back(std::make_tuple(&a_, 0, &rho_));
      op_.push_back(std::make_tuple(&b_, 0, &sigma_));
      op_.push_back(std::make_tuple(&c_, 1, &sigma_));
      op_.push_back(std::make_tuple(&d_, 1, &rho_));
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
      op_.push_back(std::make_tuple(&a_, 0,  &rho_));
      op_.push_back(std::make_tuple(&b_, 1, &rho_));
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
