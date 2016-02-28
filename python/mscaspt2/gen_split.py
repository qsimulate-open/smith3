#!/opt/local/bin/python
import string
import os


def header(n) :
    return "//\n\
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: MSCASPT2_gen" + str(n) + ".cc\n\
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
#include <src/smith/caspt2/MSCASPT2_tasks" + str(n) + ".h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
using namespace bagel::SMITH::MSCASPT2;\n\
\n\
"

footer = "#endif\n"

f = open('MSCASPT2_gen.cc', 'r')
lines = f.read().split("\n")[32:]

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
        fout = open("MSCASPT2_gen" + str(n) + ".cc", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
    num = num+1
    tmp = tmp + tasks[i];

n = (num-1) / chunk + 1
fout = open("MSCASPT2_gen" + str(n) + ".cc", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("MSCASPT2_gen.cc")
