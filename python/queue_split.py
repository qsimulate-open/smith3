#!/opt/local/bin/python
import string
import os
import re


def header(n) :
    return "//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: CASPT2" + n + ".cc\n\
// Copyright (C) 2014 Shiozaki group\n\
//\n\
// Author: Shiozaki group <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 3, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
#include <bagel_config.h>\n\
#ifdef COMPILE_SMITH\n\
\n\
\n\
#include <src/smith/CASPT2.h>\n"

def insert():
    return "#include <src/smith/CASPT2_tasks.h>\n"

def header2():
    return "\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
\n\
"

footer = "#endif\n"

f = open('CASPT2.cc', 'r')
lines = f.read().split("\n")[33:]

tasks = []
tmp = ""

for line in lines:
    if (len(line) >= 17 and (line[0:17] == "shared_ptr<Queue>" or line[0:17] == "CASPT2::CASPT2::C")):
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
    fout = open("CASPT2_" + tag + ".cc", "w")
    out = header("_" + tag + "q") + insert() + header2() + task + footer
    fout.write(out)
    fout.close()

os.remove("CASPT2.cc")

fout = open("CASPT2.cc", "w")
out = header("") + header2() + tasks[len(tasks)-1] + footer
fout.write(out)
fout.close()
