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
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <tuple>
#include "forest.h"
#include "constants.h"

using namespace std;
using namespace smith;


void Forest::filter_gamma() {
  shared_ptr<Tree> res;

  bool first = true;
  list<shared_ptr<Tensor>> prev;
  for (auto& i : trees_) {
    list<shared_ptr<Tensor>> g;
    if (first) {
      i->sort_gamma();
      res = i;
      first = false;
    } else {
      i->sort_gamma(prev);
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
    prev = i->gamma();
  }

}


OutStream Forest::generate_code() const {
  OutStream out, tmp;
  string depends, tasks, specials;

  out << generate_headers();
  out << generate_gammas();

  for (auto& i : trees_) {
    out.ss << "    std::shared_ptr<Queue> make_" << i->label() << "q(const bool reset = true, const bool diagonal = true);" << endl;

    out.ee << "shared_ptr<Queue> " << forest_name_ << "::" << forest_name_ << "::make_" << i->label() << "q(const bool reset, const bool diagonal) {" << endl << endl;
    if (i->label().find("deci") == string::npos)
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
  string forest_name_lower = forest_name_;
  transform(forest_name_lower.begin(), forest_name_lower.end(), forest_name_lower.begin(), ::tolower);
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
  if (forest_name_ == "MRCI" || forest_name_ == "RelMRCI")
    out.ss << "#include <src/smith/multitensor.h>" << endl;
  out.ss << "#include <src/smith/smith_info.h>" << endl;
  out.ss << "" << endl;
  out.ss << "namespace bagel {" << endl;
  out.ss << "namespace SMITH {" << endl;
  out.ss << "namespace " << forest_name_ << "{" << endl;
  out.ss << "" << endl;
  out.ss << "class " << forest_name_ << " : public SpinFreeMethod<" << DataType << "> {" << endl;
  out.ss << "  protected:" << endl;
  out.ss << "    std::shared_ptr<Tensor> t2;" << endl;
  out.ss << "    std::shared_ptr<Tensor> r;" << endl;
  out.ss << "    std::shared_ptr<Tensor> s;" << endl;

  if (forest_name_ == "MRCI" || forest_name_ == "RelMRCI") {
    out.ss << "    std::shared_ptr<Tensor> n;" << endl << endl;;

    out.ss << "    int nstates_;" << endl;
    out.ss << "    std::vector<double> energy_;" << endl << endl;

    out.ss << "    std::vector<std::shared_ptr<MultiTensor>> t2all_;" << endl;
    out.ss << "    std::vector<std::shared_ptr<MultiTensor>> rall_;" << endl;
    out.ss << "    std::vector<std::shared_ptr<MultiTensor>> sall_;" << endl;
    out.ss << "    std::vector<std::shared_ptr<MultiTensor>> nall_;" << endl;
  }
  if (forest_name_ == "CASPT2") {
    out.ss << "    std::shared_ptr<Tensor> den1;" << endl;
    out.ss << "    std::shared_ptr<Tensor> den2;" << endl;
    out.ss << "    std::shared_ptr<Tensor> Den1;" << endl;
    out.ss << "    double correlated_norm_;" << endl;
    out.ss << "    std::shared_ptr<Tensor> deci;" << endl << endl;
  }
  out.ss << "    void diagonal(std::shared_ptr<Tensor> r, std::shared_ptr<const Tensor> t) const;" << endl;
  out.ss << "" << endl;

  out.ee << "#include <src/util/math/davidson.h>" << endl;
  out.ee << "#include <src/smith/extrap.h>" << endl;
  out.ee << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << ".h>" << endl;
  out.ee << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.ee << "using namespace std;" << endl;
  out.ee << "using namespace bagel;" << endl;
  out.ee << "using namespace bagel::SMITH;" << endl << endl;

  out.tt << "#ifndef __SRC_SMITH_" << forest_name_ << "_" << forest_name_ << "_TASKS_H" << endl;
  out.tt << "#define __SRC_SMITH_" << forest_name_ << "_" << forest_name_ << "_TASKS_H" << endl << endl;

  out.tt << "#include <src/smith/indexrange.h>" << endl;
  out.tt << "#include <src/smith/tensor.h>" << endl;
  out.tt << "#include <src/smith/task.h>" << endl;
  out.tt << "#include <src/smith/subtask.h>" << endl;
  out.tt << "#include <src/smith/storage.h>" << endl << endl;

  out.tt << "namespace bagel {" << endl;
  out.tt << "namespace SMITH {" << endl;
  out.tt << "namespace " << forest_name_ << "{" << endl << endl;

  out.cc << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.cc << "using namespace std;" << endl;
  out.cc << "using namespace bagel;" << endl;
  out.cc << "using namespace bagel::SMITH;" << endl;
  out.cc << "using namespace bagel::SMITH::" << forest_name_ << ";" << endl << endl;

  out.dd << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.dd << "using namespace std;" << endl;
  out.dd << "using namespace bagel;" << endl;
  out.dd << "using namespace bagel::SMITH;" << endl;
  out.dd << "using namespace bagel::SMITH::" << forest_name_ << ";" << endl << endl;

  out.gg << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << ".h>" << endl;
  out.gg << "#include <src/smith/" << forest_name_lower << "/" << forest_name_ << "_tasks.h>" << endl << endl;
  out.gg << "using namespace std;" << endl;
  out.gg << "using namespace bagel;" << endl;
  out.gg << "using namespace bagel::SMITH;" << endl;
  out.gg << "using namespace bagel::SMITH::" << forest_name_ << ";" << endl << endl;

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
    vector<string> rdms = i->active()->required_rdm();
    if (i->der()) { // derivative rdm
      for (auto& j : rdms) {
        stringstream zz;
        zz << j << "deriv_";
        tmp.push_back(zz.str());
      }
    } else {  // normal rdms
      for (auto& j : rdms) {
        stringstream zz;
        zz << j << "_";
        tmp.push_back(zz.str());
      }
    }
    if (i->merged()) {
      // 4RDM derivative is a priori contracted with the fock operator
      if (!i->der() || !(rdms.size() == 1 && rdms[0] == "rdm4")) {
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
  out.ss << "    " << forest_name_ << "(std::shared_ptr<const SMITH_Info<" << DataType << ">> ref);" << endl;

  out.ee << forest_name_ << "::" << forest_name_ << "::" << forest_name_ << "(shared_ptr<const SMITH_Info<" << DataType << ">> ref) : SpinFreeMethod(ref) {" << endl;
  if (DataType == "double") {
    out.ee << "  eig_ = f1_->diag();" << endl;
  } else {
    out.ee << "  auto eig = f1_->diag();" << endl;
    out.ee << "  eig_.resize(eig.size());" << endl;
    out.ee << "  for (int i = 0; i != eig.size(); ++i)" << endl;
    out.ee << "    eig_[i] = real(eig[i]);" << endl;
  }
  if (forest_name_ == "MRCI" || forest_name_ == "RelMRCI") {
    out.ee << "  nstates_ = ref->ciwfn()->nstates();" << endl << endl;

    out.ee << "  for (int i = 0; i != nstates_; ++i) {" << endl;
    out.ee << "    auto tmp = make_shared<MultiTensor>(nstates_);" << endl;
    out.ee << "    for (auto& j : *tmp)" << endl;
    out.ee << "      j = init_amplitude();" << endl;
    out.ee << "    t2all_.push_back(tmp);" << endl << endl;

    out.ee << "    auto tmp2 = make_shared<MultiTensor>(nstates_);" << endl;
    out.ee << "    for (auto& j : *tmp2)" << endl;
    out.ee << "      j = init_residual();" << endl;
    out.ee << "    sall_.push_back(tmp2);" << endl;
    out.ee << "    nall_.push_back(tmp2->copy());" << endl;
    out.ee << "  }" << endl;
  }
  if (forest_name_ == "CASPT2" || forest_name_ == "RelCASPT2") {
    out.ee << "  t2 = init_amplitude();" << endl;
    out.ee << "  r = init_residual();" << endl;
    out.ee << "  s = init_residual();" << endl;
  }
  out.ee << "}" << endl << endl;

  out.ss << "    ~" << forest_name_ << "() {}" << endl;
  out.ss << "" << endl;
  out.ss << "    void solve();" << endl;
  out.ss << "    void solve_deriv();" << endl;

  out.ee << "void " << forest_name_ << "::" << forest_name_ << "::solve() {" << endl;

  if (forest_name_ == "CASPT2" || forest_name_ == "RelCASPT2")
    out.ee << caspt2_main_driver_();
  else if (forest_name_ == "MRCI" || forest_name_ == "RelMRCI")
    out.ee << msmrci_main_driver_();

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
    out.ee << "  den2 = h1_->clone();" << endl;
    out.ee << "  den2->allocate();" << endl;
    out.ee << "  shared_ptr<Queue> dens2 = make_densityq();" << endl;
    out.ee << "  while (!dens2->done())" << endl;
    out.ee << "    dens2->next_compute();" << endl;
    out.ee << endl;
    out.ee << "  den1 = h1_->clone();" << endl;
    out.ee << "  den1->allocate();" << endl;
    out.ee << "  shared_ptr<Queue> dens1 = make_density1q();" << endl;
    out.ee << "  while (!dens1->done())" << endl;
    out.ee << "    dens1->next_compute();" << endl;
    out.ee << endl;
    out.ee << "  Den1 = init_residual();" << endl;
    out.ee << "  shared_ptr<Queue> Dens1 = make_density2q();" << endl;
    out.ee << "  while (!Dens1->done())" << endl;
    out.ee << "    Dens1->next_compute();" << endl;
    out.ee << "  timer.tick_print(\"Correlated density matrix evaluation\");" << endl;
    out.ee << endl;
    out.ee << "  deci = make_shared<Tensor>(vector<IndexRange>{ci_});" << endl;
    out.ee << "  deci->allocate();" << endl;
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
  out.ss << "      mpi__->allreduce(&sum, 1);" << endl;
  out.ss << "      return sum;" << endl;
  out.ss << "    }" << endl;
  out.ss << endl;  // end comparison correction
  if (forest_name_ == "CASPT2") {
    out.ss << "    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }" << endl;
    out.ss << "    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }" << endl;
    out.ss << "    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }" << endl;
    out.ss << endl;
    out.ss << "    double correlated_norm() const { return correlated_norm_; }" << endl;
    out.ss << endl;
    out.ss << "    std::shared_ptr<const Civec> ci_deriv(std::shared_ptr<const Determinants> det) const { return deci->civec(det); }" << endl;
    out.ss << endl;
  }
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


string Forest::caspt2_main_driver_() {
  stringstream ss;

  ss << "  Timer timer;" << endl;
  ss << "  print_iteration();" << endl;

  ss << "  shared_ptr<Queue> sourceq = make_sourceq();" << endl;
  ss << "  while (!sourceq->done())" << endl;
  ss << "    sourceq->next_compute();" << endl;

  ss << "  Timer mtimer;" << endl;
  ss << "  int iter = 0;" << endl;
  ss << "  for ( ; iter != info_->maxiter(); ++iter) {" << endl;
  ss << "    energy_ = detail::real(dot_product_transpose(s, t2));" << endl;

  ss << "    shared_ptr<Queue> queue = make_residualq();" << endl;
  ss << "    while (!queue->done())" << endl;
  ss << "      queue->next_compute();" << endl;
  ss << "    diagonal(r, t2);" << endl;
  ss << "    r->ax_plus_y(1.0, s);" << endl;

  ss << "    energy_ += detail::real(dot_product_transpose(r, t2));" << endl;

  ss << "    const double err = r->rms();" << endl;
  ss << "    print_iteration(iter, energy_, err, mtimer.tick());" << endl;
  ss << endl;
  ss << "    update_amplitude(t2, r);" << endl;
  ss << "    r->zero();" << endl;
  ss << "    if (err < info_->thresh()) break;" << endl;
  ss << "  }" << endl;
  ss << "  print_iteration(iter == info_->maxiter());" << endl;
  ss << "  timer.tick_print(\"CASPT2 energy evaluation\");" << endl;

  ss << "  cout << \"    * CASPT2 energy : \" << fixed << setw(20) << setprecision(10) << energy_+info_->ciwfn()->energy(0) << endl;" << endl;
  return ss.str();
}

string Forest::msmrci_main_driver_() {
  stringstream ss;
  ss << "  Timer timer;" << endl;
  ss << "  print_iteration();" << endl << endl;

  ss << "  const double core_nuc = core_energy_ + info_->geom()->nuclear_repulsion();" << endl << endl;

  ss << "  // target state" << endl;
  ss << "  for (int istate = 0; istate != nstates_; ++istate) {" << endl;
  ss << "    const double refen = info_->ciwfn()->energy(istate) - core_nuc;" << endl;
  ss << "    // takes care of ref coefficients" << endl;
  ss << "    t2all_[istate]->fac(istate) = 1.0;" << endl;
  ss << "    nall_[istate]->fac(istate)  = 1.0;" << endl;
  ss << "    sall_[istate]->fac(istate)  = refen;" << endl << endl;

  ss << "    for (int jst = 0; jst != nstates_; ++jst) {" << endl;
  ss << "      set_rdm(jst, istate);" << endl;
  ss << "      s = sall_[istate]->at(jst);" << endl;
  ss << "      auto queue = make_sourceq(false, jst == istate);" << endl;
  ss << "      while (!queue->done())" << endl;
  ss << "        queue->next_compute();" << endl;
  ss << "    }" << endl;
  ss << "  }" << endl << endl;

  ss << "  DavidsonDiag_<Amplitude<" << DataType << ">, Residual<" << DataType << ">, " << MatType << "> davidson(nstates_, 10);" << endl << endl;

  ss << "  // first iteration is trivial" << endl;
  ss << "  {" << endl;
  ss << "    vector<shared_ptr<const Amplitude<" << DataType << ">>> a0;" << endl;
  ss << "    vector<shared_ptr<const Residual<" << DataType << ">>> r0;" << endl;
  ss << "    for (int istate = 0; istate != nstates_; ++istate) {" << endl;
  ss << "      a0.push_back(make_shared<Amplitude<" << DataType << ">>(t2all_[istate]->copy(), nall_[istate]->copy(), this));" << endl;
  ss << "      r0.push_back(make_shared<Residual<" << DataType << ">>(sall_[istate]->copy(), this));" << endl;
  ss << "    }" << endl;
  ss << "    energy_ = davidson.compute(a0, r0);" << endl;
  ss << "    for (int istate = 0; istate != nstates_; ++istate)" << endl;
  ss << "      assert(fabs(energy_[istate]+core_nuc - info_->ciwfn()->energy(istate)) < 1.0e-8);" << endl;
  ss << "  }" << endl << endl;

  ss << "  // set the result to t2" << endl;
  ss << "  {" << endl;
  ss << "    vector<shared_ptr<Residual<" << DataType << ">>> res = davidson.residual();" << endl;
  ss << "    for (int i = 0; i != nstates_; ++i) {" << endl;
  ss << "      t2all_[i]->zero();" << endl;
  ss << "      update_amplitude(t2all_[i], res[i]->tensor());" << endl;
  ss << "    }" << endl;
  ss << "  }" << endl << endl;

  ss << "  shared_ptr<MultiTensor> rtmp = nall_[0]->copy();" << endl;
  ss << "  rtmp->zero();" << endl << endl;

  ss << "  Timer mtimer;" << endl;
  ss << "  int iter = 0;" << endl;
  ss << "  vector<bool> conv(nstates_, false);" << endl;
  ss << "  for ( ; iter != info_->maxiter(); ++iter) {" << endl << endl;

  ss << "    // loop over state of interest" << endl;
  ss << "    vector<shared_ptr<const Amplitude<" << DataType << ">>> a0;" << endl;
  ss << "    vector<shared_ptr<const Residual<" << DataType << ">>> r0;" << endl;
  ss << "    for (int istate = 0; istate != nstates_; ++istate) {" << endl;
  ss << "      if (conv[istate]) {" << endl;
  ss << "        a0.push_back(nullptr);" << endl;
  ss << "        r0.push_back(nullptr);" << endl;
  ss << "        continue;" << endl;
  ss << "      }" << endl;
  ss << "      // first calculate left-hand-side vectors of t2 (named n)" << endl;
  ss << "      nall_[istate]->zero();" << endl;
  ss << "      for (int ist = 0; ist != nstates_; ++ist) {" << endl;
  ss << "        for (int jst = 0; jst != nstates_; ++jst) {" << endl;
  ss << "          set_rdm(jst, ist);" << endl;
  ss << "          t2 = t2all_[istate]->at(ist);" << endl;
  ss << "          n  = nall_[istate]->at(jst);" << endl;
  ss << "          auto queue = make_normq(false, jst == ist);" << endl;
  ss << "          while (!queue->done())" << endl;
  ss << "            queue->next_compute();" << endl;
  ss << "        }" << endl;
  ss << "      }" << endl << endl;

  ss << "      // normalize t2 and n" << endl;
  ss << "      const double scal = 1.0 / sqrt(detail::real(dot_product_transpose(nall_[istate], t2all_[istate])));" << endl;
  ss << "      nall_[istate]->scale(scal);" << endl;
  ss << "      t2all_[istate]->scale(scal);" << endl << endl;

  ss << "      a0.push_back(make_shared<Amplitude<" << DataType << ">>(t2all_[istate]->copy(), nall_[istate]->copy(), this));" << endl << endl;

  ss << "      // compute residuals (named r)" << endl;
  ss << "      rtmp->zero();" << endl;
  ss << "      for (int ist = 0; ist != nstates_; ++ist) { // ket sector" << endl;
  ss << "        for (int jst = 0; jst != nstates_; ++jst) { // bra sector" << endl;
  ss << "          set_rdm(jst, ist);" << endl;
  ss << "          t2 = t2all_[istate]->at(ist);" << endl;
  ss << "          r = rtmp->at(jst);" << endl;
  ss << "          auto queue = make_residualq(false, jst == ist);" << endl;
  ss << "          while (!queue->done())" << endl;
  ss << "            queue->next_compute();" << endl;
  ss << "          diagonal(r, t2);" << endl;
  ss << "        }" << endl;
  ss << "      }" << endl << endl;

  ss << "      // <ab/ij| T |0_ist> Eref_ist." << endl;
  ss << "      {" << endl;
  ss << "        shared_ptr<MultiTensor> m = t2all_[istate]->copy();" << endl;
  ss << "        for (int ist = 0; ist != nstates_; ++ist) {" << endl;
  ss << "          // First weighted T2 amplitude" << endl;
  ss << "          m->at(ist)->scale(info_->ciwfn()->energy(ist) - core_nuc);" << endl;
  ss << "          // then add it to residual" << endl;
  ss << "          for (int jst = 0; jst != nstates_; ++jst) {" << endl;
  ss << "            set_rdm(jst, ist);" << endl;
  ss << "            t2 = m->at(ist);" << endl;
  ss << "            n  = rtmp->at(jst);" << endl;
  ss << "            auto queue = make_normq(false, jst == ist);" << endl;
  ss << "            while (!queue->done())" << endl;
  ss << "              queue->next_compute();" << endl;
  ss << "          }" << endl;
  ss << "        }" << endl;
  ss << "      }" << endl << endl;

  ss << "      {" << endl;
  ss << "        shared_ptr<MultiTensor> m = rtmp->copy();" << endl;
  ss << "        for (int ist = 0; ist != nstates_; ++ist)" << endl;
  ss << "          m->fac(ist) = dot_product_transpose(sall_[ist], t2all_[istate]);" << endl;
  ss << "        r0.push_back(make_shared<Residual<" << DataType << ">>(m, this));" << endl;
  ss << "      }" << endl;
  ss << "    }" << endl << endl;

  ss << "    energy_ = davidson.compute(a0, r0);" << endl << endl;

  ss << "    // find new trial vectors" << endl;
  ss << "    vector<shared_ptr<Residual<" << DataType << ">>> res = davidson.residual();" << endl;
  ss << "    for (int i = 0; i != nstates_; ++i) {" << endl;
  ss << "      const double err = res[i]->tensor()->rms();" << endl;
  ss << "      print_iteration(iter, energy_[i]+core_nuc, err, mtimer.tick(), i);" << endl << endl;

  ss << "      t2all_[i]->zero();" << endl;
  ss << "      conv[i] = err < info_->thresh();" << endl;
  ss << "      if (!conv[i])" << endl;
  ss << "        update_amplitude(t2all_[i], res[i]->tensor());" << endl;
  ss << "    }" << endl;
  ss << "    if (nstates_ > 1) cout << endl;" << endl << endl;

  ss << "    if (all_of(conv.begin(), conv.end(), [](bool i){ return i;})) break;" << endl;
  ss << "  }" << endl;
  ss << "  print_iteration(iter == info_->maxiter());" << endl;
  ss << "  timer.tick_print(\"MRCI energy evaluation\");" << endl << endl;

  ss << "  // Davidson corrections..." << endl;
  ss << "  {" << endl;
  ss << "    cout << endl;" << endl;
  ss << "    vector<double> energy_q(nstates_);" << endl;
  ss << "    vector<shared_ptr<Amplitude<" << DataType << ">>> ci = davidson.civec();" << endl;
  ss << "    for (int i = 0; i != nstates_; ++i) {" << endl;
  ss << "      const double c = norm(ci[i]->tensor()->fac(i));" << endl;
  ss << "      const double eref = info_->ciwfn()->energy(i);" << endl;
  ss << "      const double eq = energy_[i]+core_nuc + (energy_[i]+core_nuc-eref)*(1.0-c)/c;" << endl;
  ss << "      print_iteration(0, eq, 0.0, 0.0, i);" << endl;
  ss << "    }" << endl;
  ss << "    cout << endl;" << endl;
  ss << "  }" << endl;
  ss << "  timer.tick_print(\"MRCI+Q energy evaluation\");" << endl;


  return ss.str();
}
