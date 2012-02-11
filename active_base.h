//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

// this will be moved to the newint package.

#ifndef __ACTIVE_BASE
#define __ACTIVE_BASE

class Active_base {
  protected:
    std::unique_ptr<double[]> data_;

  public:
    Active_base(const int n) {};

    double& data() {};

};

#endif
