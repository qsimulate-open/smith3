//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "constants.h"
#include "active.h"
#include <sstream>
#include <algorithm>
#include <iomanip>

using namespace std;

static int count_;

// generator for Active part
string Active::generate(shared_ptr<Tensor> merged) const {
  stringstream ss;

  if (!count_) ss << header << endl << endl;
  count__ = count_;

  // get target indices
  const list<shared_ptr<Index> > in = index(); 
  list<shared_ptr<Index> > m;
  if (merged) m = merged->index();

  const int rank = in.size() - m.size();

  if (!count_) ss << "#include \"active_base.h\"" << endl;
  ss << "\
class Active_" << count_ << " : public Active_base {\n\
  public:\n\
    Active_" << count_ << "() : Active_base(" << rank << ") {\n\
      // Here nact**" << rank << " storage is already created in unique_ptr<double[]> data_.\n\
      // Index ordering is ";

  stringstream target;
  for (auto i = in.begin(); i != in.end(); ++i) {
    if (find(m.begin(), m.end(), *i) == m.end()) {
      if (target.str().size()) target << ", ";
      target << (*i)->str(false);
    }
  }
  ss << target.str() << endl;

  // making loops
  string indent = "      ";
  vector<string> tt;
  for (auto i = in.begin(); i != in.end(); ++i, indent += "  ") {
    string cin = (*i)->str(false);
    ss << indent << "for (int " << cin << " = 0; " << cin << " != nact_; ++" << cin << ") {" << endl; 
    tt.push_back(indent + "}\n");
  }

  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) { 
    ss << indent << (*i)->str(target.str(), merged) << endl;
  }
  for (auto i = tt.rbegin(); i != tt.rend(); ++i) ss << (*i);
  ss << "\
    };\n\
    ~Active_" << count_ << "() {};\n";
  
  ss << "\
};" << endl << endl;

  ++count_;
  return ss.str();
}



string RDM::str(string target, shared_ptr<Tensor> m) const {
  stringstream ss;
  if (delta_.size()) {
    ss << "if (";
    int j = 0;
    for (auto i = delta_.begin(); i != delta_.end(); ++i, ++j) {
      if (j) ss << " && ";
      ss << i->first->str(false) << " == " << i->second->str(false);
    }
    ss << ") ";
  }
  ss << "data(" << target << ") += " << fixed << setw(5) << setprecision(2) << fac_ << " * gamma";
  ss << index_.size()/2 << "(";
  int j = index_.size()-1;
  for (auto i = index_.begin(); i != index_.end(); ++i, --j) {
    ss << (*i)->str(false);
    if (j) ss << ", ";
  }
  if (m) {
    ss << ") * " << m->str();
  } else {
    ss << ");";
  }
  return ss.str();
}

