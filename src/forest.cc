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


#include <tuple>
#include "forest.h"
#include "constants.h"

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
    out.ss << "    std::shared_ptr<Queue> make_" << i->label() << "q();" << endl;
    out.ee << "shared_ptr<Queue> " << forest_name_ << "::" << forest_name_ << "::make_" << i->label() << "q() {" << endl << endl;
    if (i->label() != "deci")
      out.ee << "  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};" << endl;
    else
      out.ee << "  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};" << endl << endl;

    tie(tmp, icnt, i0, itensors_) = i->generate_task_list(icnt, i0, gamma_, itensors_);

    out << tmp;
    out.ee << "  return " << i->label() << "q;" << endl;
    out.ee << "}" << endl << endl;
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
  out.cc << header(forest_name_ + "_gen.cc");
  out.dd << header(forest_name_ + "_tasks.cc");
  out.ee << header(forest_name_ + ".cc");
  out.gg << header(forest_name_ + "_gamma.cc");

  out.ss << "#ifndef __SRC_SMITH_" << forest_name_ << "_H" << endl;
  out.ss << "#define __SRC_SMITH_" << forest_name_ << "_H" << endl;
  out.ss << "" << endl;
  out.ss << "#include <iostream>" << endl;
  out.ss << "#include <tuple>" << endl;
  out.ss << "#include <iomanip>" << endl;
  out.ss << "#include <src/smith/spinfreebase.h>" << endl;
  out.ss << "#include <src/smith/futuretensor.h>" << endl;
  out.ss << "#include <src/scf/hf/fock.h>" << endl;
  out.ss << "#include <src/util/f77.h>" << endl;
  out.ss << "#include <src/smith/queue.h>" << endl;
  out.ss << "#include <src/smith/smith_info.h>" << endl;
  out.ss << "" << endl;
  out.ss << "namespace bagel {" << endl;
  out.ss << "namespace SMITH {" << endl;
  out.ss << "namespace " << forest_name_ << "{" << endl;
  out.ss << "" << endl;
  out.ss << "class " << forest_name_ << " : public SpinFreeMethod {" << endl;
  out.ss << "  protected:" << endl;
  out.ss << "    using SpinFreeMethod::ref_;" << endl;
  out.ss << "    using SpinFreeMethod::e0_;" << endl;
  out.ss << "    using SpinFreeMethod::closed_;" << endl;
  out.ss << "    using SpinFreeMethod::active_;" << endl;
  out.ss << "    using SpinFreeMethod::virt_;" << endl;
  out.ss << "    using SpinFreeMethod::ci_;" << endl;
  out.ss << "    using SpinFreeMethod::rclosed_;" << endl;
  out.ss << "    using SpinFreeMethod::ractive_;" << endl;
  out.ss << "    using SpinFreeMethod::rvirt_;" << endl;
  out.ss << "    using SpinFreeMethod::rci_;" << endl;
  out.ss << "    using SpinFreeMethod::h1_;" << endl;
  out.ss << "    using SpinFreeMethod::f1_;" << endl;
  out.ss << "    using SpinFreeMethod::v2_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm1_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm2_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm3_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm4_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm1deriv_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm2deriv_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm3deriv_;" << endl;
  out.ss << "    using SpinFreeMethod::rdm4deriv_;" << endl;
  out.ss << "" << endl;
  out.ss << "    std::shared_ptr<Tensor> t2;" << endl;
  out.ss << "    std::shared_ptr<Tensor> r;" << endl;
  if (forest_name_ == "MRCI")
    out.ss << "    std::shared_ptr<Tensor> s;" << endl;
  out.ss << "    double e0_;" << endl;
  out.ss << "    std::shared_ptr<Tensor> den1;" << endl;
  out.ss << "    std::shared_ptr<Tensor> den2;" << endl;
  out.ss << "    std::shared_ptr<Tensor> Den1;" << endl;
  out.ss << "    double correlated_norm_;" << endl;
  out.ss << "    std::shared_ptr<Tensor> deci;" << endl;
  out.ss << "" << endl;

  out.ee << "#include <src/smith/" << forest_name_ << ".h>" << endl;
  out.ee << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.ee << "using namespace std;" << endl;
  out.ee << "using namespace bagel;" << endl;
  out.ee << "using namespace bagel::SMITH;" << endl << endl;

  out.tt << "#ifndef __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
  out.tt << "#define __SRC_SMITH_" << forest_name_ << "_TASKS_H" << endl;
  out.tt << "" << endl;
  out.tt << "#include <src/smith/indexrange.h>" << endl;
  out.tt << "#include <src/smith/tensor.h>" << endl;
  out.tt << "#include <src/smith/task.h>" << endl;
  out.tt << "#include <src/smith/subtask.h>" << endl;
  out.tt << "#include <src/smith/storage.h>" << endl;
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

  out.dd << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.dd << "using namespace std;" << endl;
  out.dd << "using namespace bagel;" << endl;
  out.dd << "using namespace bagel::SMITH;" << endl;
  out.dd << "using namespace bagel::SMITH::" << forest_name_ << ";" << endl << endl;

  out.gg << "#include <src/smith/" << forest_name_ << ".h>" << endl;
  out.gg << "#include <src/smith/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.gg << "using namespace std;" << endl;
  out.gg << "using namespace bagel;" << endl;
  out.gg << "using namespace bagel::SMITH;" << endl << endl;

  return out;
}


OutStream Forest::generate_gammas() const {
  OutStream out;
  string indent = "      ";

  // All the gamma tensors (for all trees) should be defined here. Only distinct Gammas are computed.
  out.ss << endl;
  for (auto& i : gamma_) {

    i->set_num(icnt);
    assert(i->label().find("Gamma") != string::npos);

    out.ss << "    std::shared_ptr<FutureTensor> " << i->label() << "_();" << endl;

    out.gg << "shared_ptr<FutureTensor> " << forest_name_ << "::" << forest_name_ << "::" << i->label() << "_() {" << endl;
    out.gg << i->constructor_str() << endl;

    if (!i->der())
      out.gg << "  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};" << endl;
    else
      out.gg << "  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};" << endl;

    // switch for blas, if true merged rdm*f1 tensor multiplication will use blas
    bool use_blas = false;
    out << i->generate_gamma(icnt, use_blas, i->der());

    vector<string> tmp = {i->label()};
    vector<int> rdms = i->active()->required_rdm();
    if (i->der()) { // derivative rdm
      for (auto& j : rdms) {
        stringstream zz;
        zz << "rdm" << j << "deriv_";
        tmp.push_back(zz.str());
      }
    } else {  // normal rdms
      for (auto& j : rdms) {
        stringstream zz;
        zz << "rdm" << j << "_";
        tmp.push_back(zz.str());
      }
    }
    if (i->merged()) {
      // 4RDM derivative is a priori contracted with the fock operator
      if (!i->der() || !(rdms.size() == 1 && rdms[0] == 4)) {
        stringstream mm;
        mm << i->merged()->label() << "_";
        tmp.push_back(mm.str());
      }
    }
    // virtual generate_task
    if (i->der()) {
      out << trees_.front()->generate_task(0, icnt, tmp, "", 0, true);
    } else {
      out << trees_.front()->generate_task(0, icnt, tmp);
    }

    out.gg << "  return make_shared<FutureTensor>(*" << i->label() << ", task" << icnt << ");" << endl;
    out.gg << "}" << endl << endl;
    ++icnt;
  }

  return out;
}


OutStream Forest::generate_algorithm() const {
  OutStream out;
  string indent = "      ";

  // generate computational algorithm
  out.ss << endl;
  out.ss << "  public:" << endl;
  out.ss << "    " << forest_name_ << "(std::shared_ptr<const SMITH_Info> ref);" << endl;

  out.ee << forest_name_ << "::" << forest_name_ << "::" << forest_name_ << "(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {" << endl;
  out.ee << "  this->eig_ = f1_->diag();" << endl;
  out.ee << "  t2 = init_amplitude();" << endl;
  out.ee << "  r = t2->clone();" << endl;
  if (forest_name_ == "MRCI")
    out.ee << "  s = t2->clone();" << endl;

  out.ee << "  den1 = h1_->clone();" << endl;
  out.ee << "  den2 = h1_->clone();" << endl;
  out.ee << "  Den1 = v2_->clone();" << endl;
  out.ee << "  deci = make_shared<Tensor>(vector<IndexRange>{ci_});" << endl;
  out.ee << "}" << endl << endl;

  out.ss << "    ~" << forest_name_ << "() {}" << endl;
  out.ss << "" << endl;
  out.ss << "    void solve();" << endl;
  out.ss << "    void solve_deriv();" << endl;

  out.ee << "void " << forest_name_ << "::" << forest_name_ << "::solve() {" << endl;
  out.ee << "  Timer timer;" << endl;
  out.ee << "  this->print_iteration();" << endl;
  if (forest_name_ == "MRCI") {
    out.ee << "  shared_ptr<Queue> sq = make_sourceq();" << endl;
    out.ee << "  while (!sq->done())" << endl;
    out.ee << "    sq->next_compute();" << endl;
  }
  out.ee << "  int iter = 0;" << endl;
  out.ee << "  for ( ; iter != ref_->maxiter(); ++iter) {" << endl;
  // Compute energies here in CASPT2
  if (forest_name_ == "CASPT2") {
    out.ee << "    shared_ptr<Queue> energyq = make_energyq();" << endl;
    out.ee << "    this->energy_ = accumulate(energyq);" << endl;
  }
  out.ee << "    shared_ptr<Queue> queue = make_residualq();" << endl;
  out.ee << "    while (!queue->done())" << endl;
  out.ee << "      queue->next_compute();" << endl;
  // for CASPT2 there is manual implemention for some diagonal part
  if (forest_name_ == "CASPT2") {
    out.ee << "    diagonal(r, t2);" << endl;
    out.ee << "    this->energy_ += dot_product_transpose(r, t2);" << endl;
  } else if (forest_name_ == "MRCI") {
    out.ee << "    r->ax_plus_y(1.0, s);" << endl;
  }
  out.ee << "    const double err = r->rms();" << endl;
  if (forest_name_ == "MRCI") {
    out.ee << "    r->ax_plus_y(1.0, s);" << endl;
    out.ee << "    this->energy_ += dot_product_transpose(r, t2);" << endl;
  }
  out.ee << "    this->print_iteration(iter, this->energy_, err);" << endl;
  out.ee << endl;
  out.ee << "    this->update_amplitude(t2, r);" << endl;
  out.ee << "    r->zero();" << endl;
  out.ee << "    if (err < ref_->thresh()) break;" << endl;
  out.ee << "  }" << endl;
  out.ee << "  this->print_iteration(iter == ref_->maxiter());" << endl;
  out.ee << "  timer.tick_print(\"" << forest_name_ << " energy evaluation\");" << endl;
  out.ee << "}" << endl;
  out.ee << endl;
  out.ee << "void " << forest_name_ << "::" << forest_name_ << "::solve_deriv() {" << endl;
  // derivative is only supported in CASPT2 so far
  if (forest_name_ == "CASPT2") {
    out.ee << "  Timer timer;" << endl;
    // using norm in various places, eg  y-=Nf<I|Eij|0> and dm1 -= N*rdm1
    out.ee << "  shared_ptr<Queue> corrq = make_corrq();" << endl;
    out.ee << "  correlated_norm_ = accumulate(corrq);" << endl;
    out.ee << "  timer.tick_print(\"T1 norm evaluation\");" << endl;
    out.ee << endl;
    out.ee << "  shared_ptr<Queue> dens2 = make_densityq();" << endl;
    out.ee << "  while (!dens2->done())" << endl;
    out.ee << "    dens2->next_compute();" << endl;
    out.ee << "  shared_ptr<Queue> dens1 = make_density1q();" << endl;
    out.ee << "  while (!dens1->done())" << endl;
    out.ee << "    dens1->next_compute();" << endl;
    out.ee << "  shared_ptr<Queue> Dens1 = make_density2q();" << endl;
    out.ee << "  while (!Dens1->done())" << endl;
    out.ee << "    Dens1->next_compute();" << endl;
    out.ee << "  timer.tick_print(\"Correlated density matrix evaluation\");" << endl;
    out.ee << endl;
    out.ee << "  shared_ptr<Queue> dec = make_deciq();" << endl;
    out.ee << "  while (!dec->done())" << endl;
    out.ee << "    dec->next_compute();" << endl;
    out.ee << "  timer.tick_print(\"CI derivative evaluation\");" << endl;
    out.ee << "  cout << endl;" << endl;
  } else {
    out.ee << "  throw std::logic_error(\"Nuclear gradients not implemented for " << forest_name_ << "\");" << endl;
  }
  out.ee << "}" << endl;

  out.ss << "" << endl;
  out.ss << "    double accumulate(std::shared_ptr<Queue> queue) {" << endl;
  out.ss << "      double sum = 0.0;" << endl;
  out.ss << "      while (!queue->done())" << endl;
  out.ss << "        sum += queue->next_compute()->target();" << endl;  // prefactors included in main.cc
  out.ss << "      return sum;" << endl;
  out.ss << "    }" << endl;
  out.ss << endl;  // end comparison correction
  out.ss << "    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }" << endl;
  out.ss << "    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }" << endl;
  out.ss << "    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }" << endl;
  out.ss << endl;
  out.ss << "    double correlated_norm() const { return correlated_norm_; }" << endl;
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



