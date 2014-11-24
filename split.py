#!/opt/local/bin/python
import string

f = open('CAS_test_tasks.cc', 'r')
lines = f.read().split("\n")
chunk = 10
n = (len(lines)-34)/chunk + 1

for i in range(0,n):
    cont = lines[0:33]
    for j in range(0,n):
        if i != j:
            cont.append("#if 0")
        for k in range(33+j*chunk,min(33+(j+1)*chunk,len(lines)-1)):
            cont.append(lines[k])
        if i != j:
            cont.append("#endif")
    fout = open('CAS_test_tasks' + str(i) + '.cc', 'w')
    out = '\n'.join(cont)
    fout.write(out)
    fout.close()
