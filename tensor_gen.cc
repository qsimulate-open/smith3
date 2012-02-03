//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "tensor.h"
#include <sstream>

using namespace std;

string Tensor::generate() const {
  stringstream ss;
  // active should be generated here.
  if (active_) ss << active_->generate();

  return ss.str();
}
