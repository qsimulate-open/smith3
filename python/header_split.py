#!/opt/local/bin/python
import string
import os

def header(n) :
    return "//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: CAS_test_tasks" + str(n) + ".h\n\
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
#ifndef __SRC_SMITH_CAS_test_TASKS" + str(n) + "_H\n\
#define __SRC_SMITH_CAS_test_TASKS" + str(n) + "_H\n\
\n\
#include <src/smith/indexrange.h>\n\
#include <src/smith/tensor.h>\n\
#include <src/smith/task.h>\n\
#include <src/smith/subtask.h>\n\
#include <src/smith/storage.h>\n\
\n\
namespace bagel {\n\
namespace SMITH {\n\
namespace CAS_test{\n\
\n"

footer = "\n}\n}\n}\n\
#endif\n\
\n"


header2 = "//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: CAS_test_tasks.h\n\
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
#ifndef __SRC_SMITH_CAS_test_TASKS_H\n\
#define __SRC_SMITH_CAS_test_TASKS_H\n\
\n"

footer2 = "\n#endif\n\n"

f = open('CAS_test_tasks.h', 'r')
lines = f.read().split("\n")[39:][:-6]

tasks = []
tmp = ""

for line in lines:
    if (len(line) >= 10 and line[0:10] == "class Task"):
        if (tmp != ""):
            tasks.append(tmp)
            tmp = ""
    if (line != ""):
        tmp += line + "\n"
        if (len(line) >= 2 and line == "};"):
            tmp += "\n"
tasks.append(tmp)

tmp = ""
num = 0
chunk = 50
for i in range(len(tasks)):
    if (num != 0 and num % chunk == 0): 
        n = num / chunk
        fout = open("CAS_test_tasks" + str(n) + ".h", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
    num = num+1
    tmp = tmp + tasks[i];

n = (num-1) / chunk + 1
fout = open("CAS_test_tasks" + str(n) + ".h", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("CAS_test_tasks.h")
fout = open("CAS_test_tasks.h", "w")
out = header2
for i in range(n+1):
    if (i > 0):
        out += "#include <src/smith/CAS_test_tasks" + str(i) + ".h>\n"
out += footer2
fout.write(out)
fout.close()
