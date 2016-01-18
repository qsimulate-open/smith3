#!/opt/local/bin/python
import string
import os
import re


def header(n) :
    return "//\n\
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: MRCI_tasks" + str(n) + ".cc\n\
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
#include <src/smith/mrci/MRCI_tasks" + str(n) + ".h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
using namespace bagel::SMITH::MRCI;\n\
\n\
"

footer = "#endif\n"

f = open('MRCI_tasks.cc', 'r')
lines = f.read().split("\n")[32:]

tasks = []
tmp = ""

for line in lines:
    if (len(line) >= 9 and line[0:9] == "void Task"):
        if (tmp != ""):
            tasks.append(tmp)
            tmp = ""
    if (line != ""):
        tmp += line + "\n"
        if (line == "}"):
            tmp += "\n"
tasks.append(tmp)

p = re.compile('[0-9]+')
tmp = ""
num = 0
chunk = 50
n = 1
for task in tasks:
    num = int(p.search(task).group())
    if (num != 0 and num >= n*chunk):
        fout = open("MRCI_tasks" + str(n) + ".cc", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
        n = n+1
    tmp = tmp + task;

n = num / chunk + 1
fout = open("MRCI_tasks" + str(n) + ".cc", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("MRCI_tasks.cc")
