#!/opt/local/bin/python
import string
import os

f = open('RelMSCASPT2_gamma.cc', 'r')
lines = f.read().split("\n")
f.close()

os.remove("RelMSCASPT2_gamma.cc")

f = open('RelMSCASPT2_gamma.cc', 'w')
f.write("\n".join(lines[:25]))
add1 = "\n#include <bagel_config.h>\n#ifdef COMPILE_SMITH\n"
f.write(add1)
f.write("\n".join(lines[25:]))
add2 = "#endif\n"
f.write(add2)

