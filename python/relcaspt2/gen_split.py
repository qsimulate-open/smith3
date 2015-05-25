#!/opt/local/bin/python
import string
import os


def header(n) :
    return "//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: RelCASPT2_gen" + str(n) + ".cc\n\
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
#include <src/smith/RelCASPT2_tasks" + str(n) + ".h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
using namespace bagel::SMITH::RelCASPT2;\n\
\n\
"

footer = "#endif\n"

f = open('RelCASPT2_gen.cc', 'r')
lines = f.read().split("\n")[33:]

tasks = []
tmp = ""

for line in lines:
    if (line[0:4] == "Task"):
        if (tmp != ""):
            tasks.append(tmp)
            tmp = ""
    if (line != ""):
        tmp += line + "\n"
        if (line == "}"):
            tmp += "\n"
tasks.append(tmp)

tmp = ""
num = 0
chunk = 50
for i in range(len(tasks)):
    if (num != 0 and num % chunk == 0):
        n = num / chunk
        fout = open("RelCASPT2_gen" + str(n) + ".cc", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
    num = num+1
    tmp = tmp + tasks[i];

n = (num-1) / chunk + 1
fout = open("RelCASPT2_gen" + str(n) + ".cc", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("RelCASPT2_gen.cc")
