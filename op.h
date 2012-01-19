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
#include <list>
#include <tuple>
#include "index.h"

// base class for spin-summed operators

class Op {
  protected:
    // tensor info
    std::string label_;
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> > op_;

    // operator info
    std::vector<std::shared_ptr<Spin> > rho_;

    std::shared_ptr<Index> a_;
    std::shared_ptr<Index> b_;
    std::shared_ptr<Index> c_;
    std::shared_ptr<Index> d_;


  public:
    Op(const std::string lab, const std::string& ta, const std::string& tb, const std::string& tc, const std::string& td);
    Op(const std::string lab, const std::string& ta, const std::string& tb);
    Op() : label_("") { };
    virtual ~Op() {};

    // returns if this operator is completely contracted
    bool contracted() const;

    std::shared_ptr<Op> copy() const;

    std::string label() const { return label_; };

    void set_rho(const int i, std::shared_ptr<Spin> a) { rho_[i] = a; };
    std::vector<std::shared_ptr<Spin> >& rho() { return rho_; };
    std::shared_ptr<Spin> rho(const int i) const { return rho_.at(i); };
    const std::shared_ptr<Spin>* rho_ptr(const int i) const { return &rho_.at(i); };
    std::shared_ptr<Spin>* rho_ptr(const int i) { return &rho_.at(i); };

    std::shared_ptr<Index> a() const { return a_; };
    std::shared_ptr<Index> b() const { return b_; };
    std::shared_ptr<Index> c() const { return c_; };
    std::shared_ptr<Index> d() const { return d_; };

    const std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() const { return op_; };
    std::list<std::tuple<std::shared_ptr<Index>*, int, int> >& op() { return op_; };

    int num_nodagger() const;
    int num_dagger() const;

    // CAUTION:: this function returns the first daggered operator
    //           **AND** deletes the corresponding entry from this->op_.
    std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* > first_dagger_noactive();

    // returns if you can contract two labels
    bool contractable(std::string a, std::string b) { return a == b || a == "g" || b == "g"; };

    // returns which index to be kept when contraction is performed
    std::shared_ptr<Index>* survive(std::shared_ptr<Index>* a, std::shared_ptr<Index>* b);

    // perform a contraction (skipping first "skip" possibilities) and returns the factor to be multiplied.
    // Should be called from Diagram objects
    // fac, new, old
    std::tuple<double, std::shared_ptr<Spin>, std::shared_ptr<Spin> >
      contract(std::pair<std::shared_ptr<Index>*, std::shared_ptr<Spin>* >& dat, const int skip);

    // function to update num_ fields in Index and Spin. Should be called from Diagram objects
    void refresh_indices(std::map<std::shared_ptr<Index>, int>& dict,
                         std::map<std::shared_ptr<Index>, int>& done,
                         std::map<std::shared_ptr<Spin>, int>& spin);

    // printing this object out
    void print() const;
};



#endif
