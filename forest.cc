//
// SMITH3 - generates spin-free multireference electron correlation programs.
// Filename: forest.cc
// Copyright (C) 2013 Matthew MacLeod
//
// Author: Matthew MacLeod <matthew.macleod@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the SMITH3 package.
//
// The SMITH3 package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SMITH3 package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SMITH3 package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include "forest.h"
#include "constants.h"
#include <algorithm>
#include <stdexcept>
#include <map>
#include <tuple>

using namespace std;
using namespace smith;


void Forest::filter_gamma() {
  shared_ptr<Tree> res;

  for (auto& i : trees_) {
    list<shared_ptr<Tensor>> g;
    if (i->label() == "residual") {
      i->sort_gamma();
      res = i;
    } else {
      i->sort_gamma(res->gamma());
    }

    g = i->gamma();

    for (auto& j : g) {
      bool found = false;
      for (auto& k : gamma_) {
        if ((*j) == (*k)) {
          found = true;
          break;
        }
      }
      if (!found) gamma_.push_back(j);
    }
  }

}


pair<string, string> Forest::generate_code() const {
  stringstream ss, tt;
  string depends, tasks;
  pair<string, string> out;

  out = generate_headers();
  ss << out.first;
  tt << out.second;

  out = generate_gammas();
  ss << out.first;
  tt << out.second;

  for (auto& i : trees_) {
    tuple<string, string, int, int, vector<shared_ptr<Tensor>>> tmp = i->generate_task_list(icnt, i0, gamma_, itensors_);
    tie(depends, tasks, icnt, i0, itensors_) = tmp;
    ss << depends;
    tt << tasks;
  }

  out = generate_algorithm();
  ss << out.first;
  tt << out.second;

  return make_pair(ss.str(),tt.str());
}


pair<string, string> Forest::generate_headers() const {
  stringstream ss, tt;
  string indent = "      ";
    // save task zero
    icnt = 0;
    i0 = icnt;

    ss << header(forest_name_);
    tt << header(forest_name_ + "_tasks");

    ss << "#ifndef __SRC_SMITH_" << forest_name_ << "_H" << endl;
    ss << "#define __SRC_SMITH_" << forest_name_ << "_H" << endl;
    ss << "" << endl;
    ss << "#include <src/smith/spinfreebase.h>" << endl;
    ss << "#include <src/scf/fock.h>" << endl;
    ss << "#include <src/util/f77.h>" << endl;
    ss << "#include <iostream>" << endl;
    ss << "#include <tuple>" << endl;
    ss << "#include <iomanip>" << endl;
    ss << "#include <src/smith/queue.h>" << endl;
    ss << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl;
    ss << "#include <src/smith/smith_info.h>" << endl;
    ss << "" << endl;
    ss << "namespace bagel {" << endl;
    ss << "namespace SMITH {" << endl;
    ss << "namespace " << forest_name_ << "{" << endl;
    ss << "" << endl;
    ss << "template <typename T>" << endl;
    ss << "class " << forest_name_ << " : public SpinFreeMethod<T>{" << endl;
    ss << "  protected:" << endl;
    ss << "    using SpinFreeMethod<T>::ref_;" << endl;
    ss << "" << endl;
    ss << "    std::shared_ptr<Tensor<T>> t2;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> r;" << endl;
    ss << "    double e0_;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> sigma_;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> den1;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> den2;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> Den1;" << endl;
    ss << "    double correlated_norm;" << endl;
    ss << "    std::shared_ptr<Tensor<T>> deci;" << endl;
    ss << "" << endl;
    ss << "    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>,  std::shared_ptr<Queue<T>>,  std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {" << endl;
    ss << "      std::shared_ptr<Queue<T>> queue_(new Queue<T>());" << endl;
    ss << indent << "std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};" << endl;
    ss << indent << "std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};" << endl << endl;

    tt << "#ifndef __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
    tt << "#define __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
    tt << "" << endl;
    tt << "#include <memory>" << endl;
    tt << "#include <algorithm>" << endl;
    tt << "#include <src/smith/indexrange.h>" << endl;
    tt << "#include <src/smith/tensor.h>" << endl;
    tt << "#include <src/smith/task.h>" << endl;
    tt << "#include <src/smith/subtask.h>" << endl;
    tt << "#include <vector>" << endl;
    tt << "" << endl;
    tt << "namespace bagel {" << endl;
    tt << "namespace SMITH {" << endl;
    tt << "namespace " << forest_name_ << "{" << endl;
    tt << "" << endl;

    // virtual function, generate Task0 which zeros out the residual and starts zero level dependency queue.
    pair<string, string> rtmp = trees_.front()->create_target(indent, icnt);
    ss << rtmp.first;
    tt << rtmp.second;
    ++icnt;


  return make_pair(ss.str(), tt.str());
}


pair<string, string> Forest::generate_gammas() const {
  stringstream ss, tt;
  string indent = "      ";

  // All the gamma tensors (for all trees) should be defined here. Only distinct Gammas are computed.
  for (auto& i : gamma_) {

    i->set_num(icnt);
    assert(i->label().find("Gamma") != string::npos);

    ss << i->constructor_str(indent) << endl;

    // switch for blas, if true merged rdm*f1 tensor multiplication will use blas
    bool use_blas = false;
    tt << i->generate_gamma(icnt, use_blas, i->der());

    vector<string> tmp = {i->label()};
    vector<int> rdms = i->active()->required_rdm();
    if (i->der()) { // derivative rdm
      for (auto& j : rdms) {
        stringstream zz;
        zz << "this->rdm" << j << "deriv_";
        tmp.push_back(zz.str());
      }
    } else {  // normal rdms
      for (auto& j : rdms) {
        stringstream zz;
        zz << "this->rdm" << j << "_";
        tmp.push_back(zz.str());
      }
    }
    if (i->merged()) {
      stringstream mm;
      mm << "this->" << i->merged()->label() << "_";
      tmp.push_back(mm.str());
    }
    // virtual generate_task
    if (i->der()) {
      ss << trees_.front()->generate_task(indent, 0, icnt, tmp, "", 0, true);
    } else {
      ss << trees_.front()->generate_task(indent, 0, icnt, tmp);
    }
    ++icnt;
  }

  return make_pair(ss.str(),tt.str());
}


pair<string, string> Forest::generate_algorithm() const {
  stringstream ss, tt;
  string indent = "      ";

  // generate computational algorithm
  ss << "      return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);" << endl;
  ss << "    };" << endl;
  ss << endl;
  ss << "  public:" << endl;
  ss << "    " << forest_name_ << "(std::shared_ptr<const SMITH_Info> ref) : SpinFreeMethod<T>(ref) {" << endl;
  ss << "      this->eig_ = this->f1_->diag();" << endl;
  ss << "      t2 = this->v2_->clone();" << endl;
  ss << "      e0_ = this->e0();" << endl;
  ss << "      sigma_ = this->sigma();" << endl;
  ss << "      this->update_amplitude(t2, this->v2_, true);" << endl;
  ss << "      t2->scale(2.0);" << endl;
  ss << "      r = t2->clone();" << endl;
  ss << "      den1 = this->h1_->clone();" << endl;
  ss << "      den2 = this->h1_->clone();" << endl;
  ss << "      Den1 = this->v2_->clone();" << endl;
  ss << "      deci = this->rdm0deriv_->clone();" << endl;
  ss << "    };" << endl;
  ss << "    ~" << forest_name_ << "() {};" << endl;
  ss << "" << endl;
  ss << "    void solve() {" << endl;
  ss << "      this->print_iteration();" << endl;
  ss << "      int iter = 0;" << endl;
  ss << "      std::shared_ptr<Queue<T>> queue, energ, correct, dens2, dens1, Dens1, dec;" << endl;
  ss << "      for ( ; iter != ref_->maxiter(); ++iter) {" << endl;
  ss << "        std::tie(queue, energ, correct, dens2, dens1, Dens1, dec) = make_queue_();" << endl;
  ss << "        while (!queue->done())" << endl;
  ss << "          queue->next_compute();" << endl;
  ss << "        this->update_amplitude(t2, r);" << endl;
  ss << "        const double err = r->rms();" << endl;
  ss << "        r->zero();" << endl;
  ss << "        this->energy_ = energy(energ);" << endl;
  ss << "        this->print_iteration(iter, this->energy_, err);" << endl;
  ss << "        if (err < ref_->thresh()) break;" << endl;
  ss << "      }" << endl;
  ss << "      this->print_iteration(iter == ref_->maxiter());" << endl;
  ss << endl;
  ss << "      std::cout << \" === Computing correlated overlap, <1|1> ===\" << std::endl;" << endl;
  // using norm in various places, eg  y-=Nf<I|Eij|0> and dm1 -= N*rdm1
  ss << "      correlated_norm = correction(correct);" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << "      std::cout << \"      Norm  = \" << std::setprecision(10) << correlated_norm << std::endl;" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << endl;
  ss << "      std::cout << \" === Computing unrelaxed one-body density matrix, dm2, <1|E_pq|1>  ===\" << std::endl;" << endl;
  ss << "      while (!dens2->done())" << endl;
  ss << "        dens2->next_compute();" << endl;
  ss << "#if 0" << endl;
  ss << "      den2->print2(\"smith d1 correlated one-body density matrix dm2 \", 1.0e-5);" << endl;
  ss << "#endif" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << "      std::cout << \" === Computing unrelaxed one-body density matrix, dm1, 2<0|E_pq|1> ===\" << std::endl;" << endl;
  ss << "      while (!dens1->done())" << endl;
  ss << "        dens1->next_compute();" << endl;
  ss << "#if 0" << endl;
  ss << "      den1->print2(\"smith d1 correlated one-body density matrix dm1\", 1.0e-5);" << endl;
  ss << "#endif" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << "      std::cout << \" === Computing unrelaxed two-body density matrix, D1, <0|E_pqrs|1>  ===\" << std::endl;" << endl;
  ss << "      while (!Dens1->done())" << endl;
  ss << "        Dens1->next_compute();" << endl;
  ss << "#if 0" << endl;
  ss << "      Den1->print4(\"smith d2 correlated two-body density matrix D1\", 1.0e-5);" << endl;
  ss << "#endif" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << endl;
  ss << "      std::cout << \" === Computing cI derivative dE/dcI ===\" << std::endl;" << endl;
  ss << "      while (!dec->done())" << endl;
  ss << "        dec->next_compute();" << endl;
  ss << "      deci->print1(\"cI derivative tensor: \", 1.0e-15);" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << "      std::cout << \"      cI derivative * cI    = \" << std::setprecision(10) <<  deci->dot_product(this->rdm0deriv_) << std::endl;" << endl;
  ss << "      std::cout << \"      Expecting 2E          = \" << std::setprecision(10) <<  2.0*this->energy_ << std::endl;" << endl;
  ss << "      std::cout << std::endl;" << endl;
  ss << "" << endl;
  ss << "    };" << endl;
  ss << "" << endl;
  ss << "    double energy(std::shared_ptr<Queue<T>> energ) {" << endl;
  ss << "      double en = 0.0;" << endl;
  ss << "      while (!energ->done()) {" << endl;
  ss << "        std::shared_ptr<Task<T>> c = energ->next_compute();" << endl;
  ss << "        en += c->energy();" << endl;  // prefactors included in main.cc
  ss << "      }" << endl;
  ss << "      return en;" << endl;
  ss << "    }" << endl;
  ss << endl;
  ss << "    double correction(std::shared_ptr<Queue<T>> correct) {" << endl;
  ss << "      double n = 0.0;" << endl;
  ss << "      while (!correct->done()) {" << endl;
  ss << "        std::shared_ptr<Task<T>> c = correct->next_compute();" << endl;
  ss << "        n += c->correction();" << endl;
  ss << "      }" << endl;
  ss << "      return n;" << endl;
  ss << "    }" << endl;
  ss << endl;  // end comparison correction
  ss << "    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }" << endl;
  ss << "    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }" << endl;
  ss << "    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }" << endl;
  ss << endl;
  ss << "    double rdm1_correction() const { return correlated_norm; }" << endl;
  ss << endl;
  ss << "    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }" << endl;
  ss << endl;
  ss << "};" << endl;
  ss << endl;
  ss << "}" << endl;
  ss << "}" << endl;
  ss << "}" << endl;

  ss << "#endif" << endl << endl;

  tt << endl;
  tt << "}" << endl;
  tt << "}" << endl;
  tt << "}" << endl;
  tt << "#endif" << endl << endl;

  return make_pair(ss.str(), tt.str());
}



