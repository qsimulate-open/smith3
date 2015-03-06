#!/opt/local/bin/python
import string
import os

for i in range(37):
    if (i > 0):
        print "CASPT2_gen" + str(i) + ".cc CASPT2_tasks" + str(i) + ".cc"
