#!/usr/bin/python
import string
import os

f = open('RelCASA_gamma.cc', 'r')
lines = f.read().split("\n")
f.close()

os.remove("RelCASA_gamma.cc")

f = open('RelCASA_gamma.cc', 'w')
f.write("\n".join(lines[:25]))
add1 = "\n#include <bagel_config.h>\n#ifdef COMPILE_SMITH\n"
f.write(add1)
f.write("\n".join(lines[25:]))
add2 = "#endif\n"
f.write(add2)

