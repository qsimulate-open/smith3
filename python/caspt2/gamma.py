#!/usr/bin/python
import string
import os

f = open('CASPT2_gamma.cc', 'r')
lines = f.read().split("\n")
f.close()

os.remove("CASPT2_gamma.cc")

f = open('CASPT2_gamma.cc', 'w')
f.write("\n".join(lines[:5]))
add1 = "\n#include <bagel_config.h>\n#ifdef COMPILE_SMITH\n"
f.write(add1)
f.write("\n".join(lines[5:]))
add2 = "#endif\n"
f.write(add2)

