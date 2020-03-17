# Automator for flup.py

import os
exp_eff = [0.,0.25,0.5,0.75,1.0]
insight = [0.,0.25,0.5,0.75,1.0]

uf = 0
gs = 3
num = 10
comp = 10



for eff in exp_eff:
    for ins in insight:
        ofile = ("out-n_%d-u_%d-t_%d-c_%d-i_%f-e_%f.out" % (gs,uf,num,comp,ins,eff))
        os.system("python flup.py -n %d -u %d -t %d -c %d -i %f -l %f -o %s" % (gs,uf,num,comp,ins,eff,ofile))
