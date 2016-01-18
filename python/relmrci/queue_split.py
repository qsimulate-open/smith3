#!/opt/local/bin/python
import string
import os
import re


def header(n) :
    out = "//\n\
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: RelMRCI" + n + ".cc\n\
// Copyright (C) 2014 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// This program is free software: you can redistribute it and/or modify\n\
// it under the terms of the GNU General Public License as published by\n\
// the Free Software Foundation, either version 3 of the License, or\n\
// (at your option) any later version.\n\
//\n\
// This program is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU General Public License\n\
// along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\
//\n\
\n\
#include <bagel_config.h>\n\
#ifdef COMPILE_SMITH\n\
\n\
\n"
    if (n == ""):
        out += "\
#include <src/util/math/davidson.h>\n\
#include <src/smith/extrap.h>\n"
    out += "\
#include <src/smith/relmrci/RelMRCI.h>\n"
    return out

def insert():
    return "#include <src/smith/relmrci/RelMRCI_tasks.h>\n"

def header2():
    return "\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
\n\
"

footer = "#endif\n"

f = open('RelMRCI.cc', 'r')
lines = f.read().split("\n")[35:]

tasks = []
tmp = ""

for line in lines:
    if (len(line) >= 17 and (line[0:17] == "shared_ptr<Queue>" or line[0:17] == "RelMRCI::RelMRCI:")):
        if (tmp != ""):
            tasks.append(tmp)
            tmp = ""
    tmp += line + "\n"
    if (line == "}"):
        tmp += "\n"
tasks.append(tmp)

p = re.compile('make_[a-z0-9]+q')
for task in tasks[0:-1]:
    tag = p.search(task).group()[5:]
    fout = open("RelMRCI_" + tag + ".cc", "w")
    out = header("_" + tag + "q") + insert() + header2() + task + footer
    fout.write(out)
    fout.close()

os.remove("RelMRCI.cc")

fout = open("RelMRCI.cc", "w")
out = header("") + header2() + tasks[len(tasks)-1] + footer
fout.write(out)
fout.close()
