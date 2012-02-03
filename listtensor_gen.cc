//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "listtensor.h"
#include <sstream>

using namespace std;

string ListTensor::generate() const {
  stringstream ss;
  for (auto i = list_.begin(); i != list_.end(); ++i) {
    ss << (*i)->generate();
  }
  return ss.str();
}
