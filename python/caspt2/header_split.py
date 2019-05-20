#!/usr/bin/python
import string
import os
import re

def header(n) :
    return "//\n\
// Copyright (C) 2019 Quantum Simulation Technologies, Inc. - All Rights Reserved\n\
//\n\
\n\
#include <bagel_config.h>\n\
#ifdef COMPILE_SMITH\n\
\n\
#ifndef __SRC_SMITH_CASPT2_TASKS" + str(n) + "_H\n\
#define __SRC_SMITH_CASPT2_TASKS" + str(n) + "_H\n\
\n\
#include <src/smith/indexrange.h>\n\
#include <src/smith/tensor.h>\n\
#include <src/smith/task.h>\n\
#include <src/smith/subtask.h>\n\
#include <src/smith/storage.h>\n\
\n\
namespace bagel {\n\
namespace SMITH {\n\
namespace CASPT2{\n\
\n"

footer = "\n}\n}\n}\n\
#endif\n\
#endif\n\
\n"


header2 = "//\n\
// Copyright (C) 2019 Quantum Simulation Technologies, Inc. - All Rights Reserved\n\
//\n\
\n\
#include <bagel_config.h>\n\
#ifdef COMPILE_SMITH\n\
\n\
#ifndef __SRC_SMITH_CASPT2_TASKS_H\n\
#define __SRC_SMITH_CASPT2_TASKS_H\n\
\n"

footer2 = "\n#endif\n#endif\n\n"

f = open('CASPT2_tasks.h', 'r')
lines = f.read().split("\n")[18:][:-6]

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

p = re.compile('[0-9]+')
tmp = ""
num = 0
chunk = 50
n = 1
for task in tasks:
    num = int(p.search(task).group())
    if (num != 0 and num >= n*chunk):
        fout = open("CASPT2_tasks" + str(n) + ".h", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
        n = n+1
    tmp = tmp + task;

n = (num-1) / chunk + 1
fout = open("CASPT2_tasks" + str(n) + ".h", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("CASPT2_tasks.h")
fout = open("CASPT2_tasks.h", "w")
out = header2
for i in range(n+1):
    if (i > 0):
        out += "#include <src/smith/caspt2/CASPT2_tasks" + str(i) + ".h>\n"
out += footer2
fout.write(out)
fout.close()
