#!/opt/local/bin/python
import string
import os

for i in range(36):
    if (i > 0):
        print "CAS_test_gen" + str(i) + ".cc CAS_test_tasks" + str(i) + ".cc"
