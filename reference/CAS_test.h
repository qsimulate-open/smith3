//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test.h
// Copyright (C) 2012 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_CAS_test_H 
#define __SRC_SMITH_CAS_test_H 

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/CAS_test_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

template <typename T>
class CAS_test : public SpinFreeMethod<T>, SMITH_info {
  protected:
    std::shared_ptr<Tensor<T>> t2;
    std::shared_ptr<Tensor<T>> r;
    double e0_;
    std::shared_ptr<Tensor<T>> ci_coeff_;

    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
      std::shared_ptr<Queue<T>> queue_(new Queue<T>());
      std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};
      std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};

      std::vector<std::shared_ptr<Tensor<T>>> tensor0 = {r};
      std::shared_ptr<Task0<T>> task0(new Task0<T>(tensor0));
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index;
      std::shared_ptr<Tensor<T>> Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor1 = {Gamma0, this->rdm1_, this->f1_};
      std::shared_ptr<Task1<T>> task1(new Task1<T>(tensor1, pindex));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma2(new Tensor<T>(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor2 = {Gamma2, this->rdm1_};
      std::shared_ptr<Task2<T>> task2(new Task2<T>(tensor2, pindex));
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> I0_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor3 = {r, I0};
      std::shared_ptr<Task3<T>> task3(new Task3<T>(tensor3, pindex));
      task3->add_dep(task0);
      queue_->add_task(task3);


      std::vector<std::shared_ptr<Tensor<T>>> tensor4 = {I0, t2, this->v2_};
      std::shared_ptr<Task4<T>> task4(new Task4<T>(tensor4, pindex, this->e0_));
      task3->add_dep(task4);
      task4->add_dep(task0);
      queue_->add_task(task4);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor5 = {I0, t2, I1};
      std::shared_ptr<Task5<T>> task5(new Task5<T>(tensor5, pindex));
      task3->add_dep(task5);
      task5->add_dep(task0);
      queue_->add_task(task5);


      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {I1, Gamma0};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, pindex));
      task5->add_dep(task6);
      task6->add_dep(task0);
      queue_->add_task(task6);

      task6->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {I0, t2, I3};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, pindex));
      task3->add_dep(task7);
      task7->add_dep(task0);
      queue_->add_task(task7);


      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {I3, Gamma0};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, pindex));
      task7->add_dep(task8);
      task8->add_dep(task0);
      queue_->add_task(task8);

      task8->add_dep(task1);

      std::vector<IndexRange> I4_index = {this->closed_, this->virt_, this->virt_, this->closed_};
      std::shared_ptr<Tensor<T>> I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {r, I4};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, pindex));
      task9->add_dep(task0);
      queue_->add_task(task9);


      std::vector<IndexRange> I5_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {I4, this->f1_, I5};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, pindex));
      task9->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {I5, t2};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, pindex));
      task10->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<IndexRange> I9_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {I4, this->f1_, I9};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, pindex));
      task9->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I9, t2};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, pindex));
      task12->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I16_index;
      std::shared_ptr<Tensor<T>> I16(new Tensor<T>(I16_index, false));
      std::vector<IndexRange> I17_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I16, t2, I17};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, pindex));
      energy_->add_task(task14);


      std::vector<std::shared_ptr<Tensor<T>>> tensor15 = {I17, this->v2_};
      std::shared_ptr<Task15<T>> task15(new Task15<T>(tensor15, pindex));
      task14->add_dep(task15);
      energy_->add_task(task15);

      // TODO new
      std::shared_ptr<Queue<T>> dedci_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor016 = {this->deci_}; 
      std::shared_ptr<Task016<T>> task016(new Task016<T>(tensor016));
      dedci_->add_task(task016);

      std::vector<IndexRange> I21_index = {this->closed_, this->virt_, this->closed_, this->virt_, this->ci_};
      std::shared_ptr<Tensor<T>> I21(new Tensor<T>(I21_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor16 = {this->deci_, t2, I21};
      std::shared_ptr<Task16<T>> task16(new Task16<T>(tensor16, cindex));
      task16->add_dep(task016);
      dedci_->add_task(task16);

      std::vector<IndexRange> I21b_index = {this->closed_, this->virt_, this->closed_, this->virt_}; 
      std::shared_ptr<Tensor<T>> I21b(new Tensor<T>(I21b_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor017 = {I21, I21b, this->ci_coeff_};
      std::shared_ptr<Task017<T>> task017(new Task017<T>(tensor017, cindex));
      task16->add_dep(task017);
      task017->add_dep(task016);
      dedci_->add_task(task017);

      std::vector<std::shared_ptr<Tensor<T>>> tensor17 = {I21b, this->v2_};
      std::shared_ptr<Task17<T>> task17(new Task17<T>(tensor17, cindex));
      task017->add_dep(task17);
      task17->add_dep(task016);
      dedci_->add_task(task17);
      // end new

      std::shared_ptr<Queue<T>> correction_(new Queue<T>());
      std::vector<IndexRange> I24_index;
      std::shared_ptr<Tensor<T>> I24(new Tensor<T>(I24_index, false));
      std::vector<IndexRange> I25_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor18 = {I24, t2, I25};
      std::shared_ptr<Task18<T>> task18(new Task18<T>(tensor18, pindex));
      correction_->add_task(task18);


      std::vector<std::shared_ptr<Tensor<T>>> tensor19 = {I25, t2};
      std::shared_ptr<Task19<T>> task19(new Task19<T>(tensor19, pindex));
      task18->add_dep(task19);
      correction_->add_task(task19);

      std::shared_ptr<Queue<T>> density_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor20 = {this->den1_};
      std::shared_ptr<Task20<T>> task20(new Task20<T>(tensor20));
      density_->add_task(task20);

      std::vector<IndexRange> I28_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I28(new Tensor<T>(I28_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor21 = {this->den1_, I28};
      std::shared_ptr<Task21<T>> task21(new Task21<T>(tensor21, pindex));
      task21->add_dep(task20);
      density_->add_task(task21);


      std::vector<IndexRange> I29_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor22 = {I28, t2, I29};
      std::shared_ptr<Task22<T>> task22(new Task22<T>(tensor22, pindex));
      task21->add_dep(task22);
      task22->add_dep(task20);
      density_->add_task(task22);


      std::vector<IndexRange> I30_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I30(new Tensor<T>(I30_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor23 = {I29, t2, I30};
      std::shared_ptr<Task23<T>> task23(new Task23<T>(tensor23, pindex));
      task22->add_dep(task23);
      task23->add_dep(task20);
      density_->add_task(task23);


      std::vector<std::shared_ptr<Tensor<T>>> tensor24 = {I30, Gamma2};
      std::shared_ptr<Task24<T>> task24(new Task24<T>(tensor24, pindex));
      task23->add_dep(task24);
      task24->add_dep(task20);
      density_->add_task(task24);

      task24->add_dep(task2);

      std::vector<IndexRange> I33_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I33(new Tensor<T>(I33_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor25 = {I29, t2, I33};
      std::shared_ptr<Task25<T>> task25(new Task25<T>(tensor25, pindex));
      task22->add_dep(task25);
      task25->add_dep(task20);
      density_->add_task(task25);


      std::vector<std::shared_ptr<Tensor<T>>> tensor26 = {I33, Gamma2};
      std::shared_ptr<Task26<T>> task26(new Task26<T>(tensor26, pindex));
      task25->add_dep(task26);
      task26->add_dep(task20);
      density_->add_task(task26);

      task26->add_dep(task2);

      std::vector<IndexRange> I34_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I34(new Tensor<T>(I34_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor27 = {this->den1_, I34};
      std::shared_ptr<Task27<T>> task27(new Task27<T>(tensor27, pindex));
      task27->add_dep(task20);
      density_->add_task(task27);


      std::vector<IndexRange> I35_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor28 = {I34, t2, I35};
      std::shared_ptr<Task28<T>> task28(new Task28<T>(tensor28, pindex));
      task27->add_dep(task28);
      task28->add_dep(task20);
      density_->add_task(task28);


      std::vector<std::shared_ptr<Tensor<T>>> tensor29 = {I35, t2};
      std::shared_ptr<Task29<T>> task29(new Task29<T>(tensor29, pindex));
      task28->add_dep(task29);
      task29->add_dep(task20);
      density_->add_task(task29);


      std::vector<IndexRange> I38_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I38(new Tensor<T>(I38_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor30 = {this->den1_, I38};
      std::shared_ptr<Task30<T>> task30(new Task30<T>(tensor30, pindex));
      task30->add_dep(task20);
      density_->add_task(task30);


      std::vector<IndexRange> I39_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I39(new Tensor<T>(I39_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor31 = {I38, t2, I39};
      std::shared_ptr<Task31<T>> task31(new Task31<T>(tensor31, pindex));
      task30->add_dep(task31);
      task31->add_dep(task20);
      density_->add_task(task31);


      std::vector<std::shared_ptr<Tensor<T>>> tensor32 = {I39, t2};
      std::shared_ptr<Task32<T>> task32(new Task32<T>(tensor32, pindex));
      task31->add_dep(task32);
      task32->add_dep(task20);
      density_->add_task(task32);


// dens2 start
      // task 0
      std::shared_ptr<Queue<T>> density2_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor33 = {this->den2_};
      std::shared_ptr<Task33<T>> task33(new Task33<T>(tensor33));
      density2_->add_task(task33);

      std::vector<IndexRange> I42_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I42(new Tensor<T>(I42_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor34 = {this->den2_, I42};
      std::shared_ptr<Task34<T>> task34(new Task34<T>(tensor34, pindex));
      task34->add_dep(task33);
      density2_->add_task(task34);


      std::vector<std::shared_ptr<Tensor<T>>> tensor35 = {I42, t2};
      std::shared_ptr<Task35<T>> task35(new Task35<T>(tensor35, pindex));
      task34->add_dep(task35);
      task35->add_dep(task33);
      density2_->add_task(task35);
// dens2 end

      return make_tuple(queue_, energy_, dedci_, density_, correction_, density2_);
    };

  public:
    CAS_test(std::shared_ptr<const Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info() {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->e0();
      ci_coeff_ = this->ci_coeff();
      this->update_amplitude(t2, this->v2_, true);
      t2->scale(2.0);
      r = t2->clone();
      this->den1_ = this->h1_->clone();
      this->den2_ = this->v2_->clone();
      this->deci_ = this->dci_->clone();
    };
    ~CAS_test() {}; 

    void solve() {
      this->print_iteration();
      int iter = 0;
      std::shared_ptr<Queue<T>> queue, energ, dec, dens, correct, dens2;
      for ( ; iter != maxiter_; ++iter) {
        std::tie(queue, energ, dec, dens, correct, dens2) = make_queue_();
        while (!queue->done())
          queue->next_compute();
        this->update_amplitude(t2, r);
        const double err = r->rms();
        r->zero();
        const double en = energy(energ);
        this->print_iteration(iter, en, err);
        if (err < thresh_residual()) break;
      }
      this->print_iteration(iter == maxiter_);
      std::cout << " === Calculating CI derivative dE/dcI ===" << std::endl; 
      while (!dec->done())
        dec->next_compute();
      this->deci_->print1("CI derivative tensor: ", 1.0e-15);
      std::cout << "CI derivative norm: " << std::setprecision(10) <<  this->deci_->norm() << std::endl;
      std::cout << std::endl;

      std::cout << " === Unrelaxed density matrix, dm1, <1|E_pq|1> + 2<0|E_pq|1> ===" << std::endl; 
      while (!dens->done())
        dens->next_compute();
#if 0
      this->den1_->print2("density matrix", 1.0e-5);
#endif
      this->correct_den1_ = correction(correct);
      std::cout << "Unlinked correction term, <1|1>*rdm1 = " << std::setprecision(10) << this->correct_den1_ << "*rdm1" << std::endl;
#if 0
      this->rdm1_->print2("rdm1", 1.0e-5);
#endif
      std::cout << " === Unrelaxed density matrix, dm2, <0|E_pqrs|1>  ===" << std::endl; 
      while (!dens2->done())
        dens2->next_compute();
#if 0
      this->den2_->print4("density matrix", 1.0e-5);
#endif
    };

    double energy(std::shared_ptr<Queue<T>> energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task<T>> c = energ->next_compute();
        en += c->energy();
      }   
      return en; 
    };  


    double correction(std::shared_ptr<Queue<T>> correct) {
      double n = 0.0;
      while (!correct->done()) {
        std::shared_ptr<Task<T>> c = correct->next_compute();
        n += c->correction();
      }   
      return n; 
    };  

};

}
}
}
#endif

