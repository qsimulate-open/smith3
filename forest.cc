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


OutStream Forest::generate_code() const {
  OutStream out, tmp;
  string depends, tasks, specials;

  out << generate_headers();
  out << generate_gammas();

  for (auto& i : trees_) {
    tie(tmp, icnt, i0, itensors_) = i->generate_task_list(icnt, i0, gamma_, itensors_);
    out << tmp;
  }

  out << generate_algorithm();

  return out; 
}


OutStream Forest::generate_headers() const {
  OutStream out;
  string indent = "      ";
  // save task zero
  icnt = 0;
  i0 = icnt;

  out.ss << header(forest_name_ + ".h");
  out.tt << header(forest_name_ + "_tasks.h");
  out.cc << header(forest_name_ + "_tasks.cc");

  out.ss << "#ifndef __SRC_SMITH_" << forest_name_ << "_H" << endl;
  out.ss << "#define __SRC_SMITH_" << forest_name_ << "_H" << endl;
  out.ss << "" << endl;
  out.ss << "#include <iostream>" << endl;
  out.ss << "#include <tuple>" << endl;
  out.ss << "#include <iomanip>" << endl;
  out.ss << "#include <src/smith/spinfreebase.h>" << endl;
  out.ss << "#include <src/scf/hf/fock.h>" << endl;
  out.ss << "#include <src/util/f77.h>" << endl;
  out.ss << "#include <src/smith/queue.h>" << endl;
  out.ss << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl;
  out.ss << "#include <src/smith/smith_info.h>" << endl;
  out.ss << "" << endl;
  out.ss << "namespace bagel {" << endl;
  out.ss << "namespace SMITH {" << endl;
  out.ss << "namespace " << forest_name_ << "{" << endl;
  out.ss << "" << endl;
  out.ss << "class " << forest_name_ << " : public SpinFreeMethod {" << endl;
  out.ss << "  protected:" << endl;
  out.ss << "    using SpinFreeMethod::ref_;" << endl;
  out.ss << "" << endl;
  out.ss << "    std::shared_ptr<Tensor> t2;" << endl;
  out.ss << "    std::shared_ptr<Tensor> r;" << endl;
  out.ss << "    double e0_;" << endl;
  out.ss << "    std::shared_ptr<Tensor> sigma_;" << endl;
  out.ss << "    std::shared_ptr<Tensor> den1;" << endl;
  out.ss << "    std::shared_ptr<Tensor> den2;" << endl;
  out.ss << "    std::shared_ptr<Tensor> Den1;" << endl;
  out.ss << "    double correlated_norm;" << endl;
  out.ss << "    std::shared_ptr<Tensor> deci;" << endl;
  out.ss << "" << endl;
  out.ss << "    std::tuple<std::shared_ptr<Queue>, std::shared_ptr<Queue>, std::shared_ptr<Queue>,  std::shared_ptr<Queue>,  std::shared_ptr<Queue>, std::shared_ptr<Queue>, std::shared_ptr<Queue>> make_queue_() {" << endl;
  out.ss << "      auto queue_ = std::make_shared<Queue>();" << endl;
  out.ss << indent << "std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};" << endl;
  out.ss << indent << "std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};" << endl << endl;

  out.tt << "#ifndef __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
  out.tt << "#define __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
  out.tt << "" << endl;
  out.tt << "#include <memory>" << endl;
  out.tt << "#include <algorithm>" << endl;
  out.tt << "#include <src/smith/indexrange.h>" << endl;
  out.tt << "#include <src/smith/tensor.h>" << endl;
  out.tt << "#include <src/smith/task.h>" << endl;
  out.tt << "#include <src/smith/subtask.h>" << endl;
  out.tt << "#include <src/smith/storage.h>" << endl;
  out.tt << "#include <vector>" << endl;
  out.tt << "" << endl;
  out.tt << "namespace bagel {" << endl;
  out.tt << "namespace SMITH {" << endl;
  out.tt << "namespace " << forest_name_ << "{" << endl;
  out.tt << "" << endl;

  out.cc << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl << endl; 
  out.cc << "using namespace std;" << endl;
  out.cc << "using namespace bagel;" << endl;
  out.cc << "using namespace bagel::SMITH;" << endl;
  out.cc << "using namespace bagel::SMITH::" << forest_name_ << ";" << endl << endl;

  // virtual function, generate Task0 which zeros out the residual and starts zero level dependency queue.
  out << trees_.front()->create_target(indent, icnt);
  ++icnt;


  return out; 
}


OutStream Forest::generate_gammas() const {
  OutStream out;
  string indent = "      ";

  // All the gamma tensors (for all trees) should be defined here. Only distinct Gammas are computed.
  for (auto& i : gamma_) {

    i->set_num(icnt);
    assert(i->label().find("Gamma") != string::npos);

    out.ss << i->constructor_str(indent) << endl;

    // switch for blas, if true merged rdm*f1 tensor multiplication will use blas
    bool use_blas = false;
    out << i->generate_gamma(icnt, use_blas, i->der());

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
      out << trees_.front()->generate_task(indent, 0, icnt, tmp, "", 0, true);
    } else {
      out << trees_.front()->generate_task(indent, 0, icnt, tmp);
    }
    ++icnt;
  }

  return out; 
}


OutStream Forest::generate_algorithm() const {
  OutStream out;
  string indent = "      ";

  // generate computational algorithm
  out.ss << "      return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);" << endl;
  out.ss << "    };" << endl;
  out.ss << endl;
  out.ss << "  public:" << endl;
  out.ss << "    " << forest_name_ << "(std::shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {" << endl;
  out.ss << "      this->eig_ = this->f1_->diag();" << endl;
  out.ss << "      t2 = this->v2_->clone();" << endl;
  out.ss << "      e0_ = this->e0();" << endl;
  out.ss << "      sigma_ = this->sigma();" << endl;
  out.ss << "      this->update_amplitude(t2, this->v2_, true);" << endl;
  out.ss << "      t2->scale(2.0);" << endl;
  out.ss << "      r = t2->clone();" << endl;
  out.ss << "      den1 = this->h1_->clone();" << endl;
  out.ss << "      den2 = this->h1_->clone();" << endl;
  out.ss << "      Den1 = this->v2_->clone();" << endl;
  out.ss << "      deci = this->rdm0deriv_->clone();" << endl;
  out.ss << "    };" << endl;
  out.ss << "    ~" << forest_name_ << "() {};" << endl;
  out.ss << "" << endl;
  out.ss << "    void solve() {" << endl;
  out.ss << "      this->print_iteration();" << endl;
  out.ss << "      int iter = 0;" << endl;
  out.ss << "      std::shared_ptr<Queue> queue, energ, correct, dens2, dens1, Dens1, dec;" << endl;
  out.ss << "      for ( ; iter != ref_->maxiter(); ++iter) {" << endl;
  out.ss << "        std::tie(queue, energ, correct, dens2, dens1, Dens1, dec) = make_queue_();" << endl;
  out.ss << "        while (!queue->done())" << endl;
  out.ss << "          queue->next_compute();" << endl;
  out.ss << "        this->update_amplitude(t2, r);" << endl;
  out.ss << "        const double err = r->rms();" << endl;
  out.ss << "        r->zero();" << endl;
  out.ss << "        this->energy_ = energy(energ);" << endl;
  out.ss << "        this->print_iteration(iter, this->energy_, err);" << endl;
  out.ss << "        if (err < ref_->thresh()) break;" << endl;
  out.ss << "      }" << endl;
  out.ss << "      this->print_iteration(iter == ref_->maxiter());" << endl;
  out.ss << endl;
  out.ss << "      std::cout << \" === Computing correlated overlap, <1|1> ===\" << std::endl;" << endl;
  // using norm in various places, eg  y-=Nf<I|Eij|0> and dm1 -= N*rdm1
  out.ss << "      correlated_norm = correction(correct);" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << "      std::cout << \"      Norm  = \" << std::setprecision(10) << correlated_norm << std::endl;" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << endl;
  out.ss << "      std::cout << \" === Computing unrelaxed one-body density matrix, dm2, <1|E_pq|1>  ===\" << std::endl;" << endl;
  out.ss << "      while (!dens2->done())" << endl;
  out.ss << "        dens2->next_compute();" << endl;
  out.ss << "#if 0" << endl;
  out.ss << "      den2->print2(\"smith d1 correlated one-body density matrix dm2 \", 1.0e-5);" << endl;
  out.ss << "#endif" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << "      std::cout << \" === Computing unrelaxed one-body density matrix, dm1, 2<0|E_pq|1> ===\" << std::endl;" << endl;
  out.ss << "      while (!dens1->done())" << endl;
  out.ss << "        dens1->next_compute();" << endl;
  out.ss << "#if 0" << endl;
  out.ss << "      den1->print2(\"smith d1 correlated one-body density matrix dm1\", 1.0e-5);" << endl;
  out.ss << "#endif" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << "      std::cout << \" === Computing unrelaxed two-body density matrix, D1, <0|E_pqrs|1>  ===\" << std::endl;" << endl;
  out.ss << "      while (!Dens1->done())" << endl;
  out.ss << "        Dens1->next_compute();" << endl;
  out.ss << "#if 0" << endl;
  out.ss << "      Den1->print4(\"smith d2 correlated two-body density matrix D1\", 1.0e-5);" << endl;
  out.ss << "#endif" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << endl;
  out.ss << "      std::cout << \" === Computing cI derivative dE/dcI ===\" << std::endl;" << endl;
  out.ss << "      while (!dec->done())" << endl;
  out.ss << "        dec->next_compute();" << endl;
  out.ss << "      deci->print1(\"cI derivative tensor: \", 1.0e-15);" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << "      std::cout << \"      cI derivative * cI    = \" << std::setprecision(10) <<  deci->dot_product(this->rdm0deriv_) << std::endl;" << endl;
  out.ss << "      std::cout << \"      Expecting 2E          = \" << std::setprecision(10) <<  2.0*this->energy_ << std::endl;" << endl;
  out.ss << "      std::cout << std::endl;" << endl;
  out.ss << "" << endl;
  out.ss << "    };" << endl;
  out.ss << "" << endl;
  out.ss << "    double energy(std::shared_ptr<Queue> energ) {" << endl;
  out.ss << "      double en = 0.0;" << endl;
  out.ss << "      while (!energ->done()) {" << endl;
  out.ss << "        std::shared_ptr<Task> c = energ->next_compute();" << endl;
  out.ss << "        en += c->energy();" << endl;  // prefactors included in main.cc
  out.ss << "      }" << endl;
  out.ss << "      return en;" << endl;
  out.ss << "    }" << endl;
  out.ss << endl;
  out.ss << "    double correction(std::shared_ptr<Queue> correct) {" << endl;
  out.ss << "      double n = 0.0;" << endl;
  out.ss << "      while (!correct->done()) {" << endl;
  out.ss << "        std::shared_ptr<Task> c = correct->next_compute();" << endl;
  out.ss << "        n += c->correction();" << endl;
  out.ss << "      }" << endl;
  out.ss << "      return n;" << endl;
  out.ss << "    }" << endl;
  out.ss << endl;  // end comparison correction
  out.ss << "    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }" << endl;
  out.ss << "    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }" << endl;
  out.ss << "    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }" << endl;
  out.ss << endl;
  out.ss << "    double rdm1_correction() const { return correlated_norm; }" << endl;
  out.ss << endl;
  out.ss << "    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }" << endl;
  out.ss << endl;
  out.ss << "};" << endl;
  out.ss << endl;
  out.ss << "}" << endl;
  out.ss << "}" << endl;
  out.ss << "}" << endl;

  out.ss << "#endif" << endl << endl;

  out.tt << endl;
  out.tt << "}" << endl;
  out.tt << "}" << endl;
  out.tt << "}" << endl;
  out.tt << "#endif" << endl << endl;

  return out; 
}



