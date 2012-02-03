//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "tree.h"
#include <sstream>

using namespace std;

string Tree::generate() const {
  stringstream ss;
  for (auto i = tensor_.begin(); i != tensor_.end(); ++i) {
    ss << (*i)->generate();
  }
  return ss.str();
}

