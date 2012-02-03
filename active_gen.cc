//
// Author : Toru Shiozaki
// Date   : Feb 2012
//

#include "constants.h"
#include "active.h"
#include <sstream>
#include <iomanip>

using namespace std;

static int count;

// generator for Active part
string Active::generate() const {
  stringstream ss;

  if (!count) ss << header << endl << endl;
  count_ = count;

  // get target indices
  const list<shared_ptr<Index> > in = index(); 
  const int rank = in.size();

  if (!count) ss << "#include \"active_base.h\"" << endl;
  ss << "\
class Active_" << count << " : public Active_base {\n\
  public:\n\
    Active_" << count << "() : Active_base(" << rank << ") {\n\
      // Here nact**" << count << " storage is already created in unique_ptr<double[]> data_.\n\
      // Index ordering is ";
  for (auto i = in.begin(); i != in.end(); ++i) ss << (*i)->str(false) << " "; 
  ss << endl;

  // making loops
  string indent = "      ";
  vector<string> tt;
  for (auto i = in.begin(); i != in.end(); ++i, indent += "  ") {
    string cin = (*i)->str(false);
    ss << indent << "for (int " << cin << " = 0; " << cin << " != nact_; ++" << cin << ") {" << endl; 
    tt.push_back(indent + "}\n");
  }

  for (auto i = rdm_.begin(); i != rdm_.end(); ++i) { 
    ss << indent << (*i)->str() << endl;
  }
  ss << indent << "++cnt" << endl; 
  for (auto i = tt.rbegin(); i != tt.rend(); ++i) ss << (*i);
  ss << "\
    };\n\
    ~Active_" << count << "() {};\n";
  
  ss << "\
};" << endl << endl;

  ++count;
  return ss.str();
}



string RDM::str() const {
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
  ss << "data_[cnt] += " << fixed << setw(5) << setprecision(2) << fac_ << " * gamma";
  ss << index_.size()/2 << "(";
  int j = index_.size()-1;
  for (auto i = index_.begin(); i != index_.end(); ++i, --j) {
    ss << (*i)->str(false);
    if (j) ss << ", ";
  }
  ss << ");";
  return ss.str();
}

