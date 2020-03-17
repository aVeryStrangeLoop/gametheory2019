#Automator for flup.py
import os

gss = [3,4,5]
comps = [5,10,15,20]

uf = 2
eff = 0.3
num = 10
ins = 0.5



for gs in gss:
    for comp in comps:
        ofile = ("out-n_%d-u_%d-t_%d-c_%d-i_%f-e_%f.out" % (gs,uf,num,comp,ins,eff))
        os.system("python flup.py -n %d -u %d -t %d -c %d -i %f -l %f -o %s" % (gs,uf,num,comp,ins,eff,ofile))
