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
#include <src/smith/caspt2/CASPT2_tasks" + str(n) + ".h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
using namespace bagel::SMITH;\n\
using namespace bagel::SMITH::CASPT2;\n\
\n\
"

footer = "#endif\n"

f = open('CASPT2_tasks.cc', 'r')
lines = f.read().split("\n")[12:]

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
        fout = open("CASPT2_tasks" + str(n) + ".cc", "w")
        out = header(n) + tmp + footer
        fout.write(out)
        fout.close()
        tmp = ""
        n = n+1
    tmp = tmp + task;

n = (num-1) / chunk + 1
fout = open("CASPT2_tasks" + str(n) + ".cc", "w")
out = header(n) + tmp + footer
fout.write(out)
fout.close()

os.remove("CASPT2_tasks.cc")
