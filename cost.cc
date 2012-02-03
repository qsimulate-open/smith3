//
// Author : Toru Shiozaki
// Date   : Feb 2009
//

#include <algorithm>
#include <sstream>
#include "cost.h"

using namespace std;


const string PCost::show() const {
  stringstream out;
  auto j = pcost_.begin();
  for (auto i = indmap_.begin(); i != indmap_.end(); ++i, ++j)
    out << i->first <<  *j; 
  return out.str();
}


void Cost::sort_pcost() {
  sort(cost_.rbegin(), cost_.rend());
}


const string Cost::show() const {
  string out;
  for (auto i = cost_.begin(); i != cost_.end(); ++i) out += i->show() + " ";
  return out;
}


