import os
import re
repl = 0
directory = 'soft_params/repl'+str(repl) #Directory over which to iterate 

fin = open("fin.csv", "w+")
fin.write("Time, Replicate, UF, ExpEff, CompCap, NumTasks, Grid, Insight, EI, Entropy, Fitness\n") #One time header

uf = None #Declare all variables to be able to access them later
grid = None
comp = None
exp = None
ins = None
rep = repl
run = None
entropy = None
fitarray = []
EIarray = []
fitness = None
EI = None
num = 10

for filename in os.listdir(directory): #Iterate over files in directory

    cur_file = open(filename, "r+")
    filestr = str(filename) #Filename in str
    content = [line.rstrip('\n') for line in open(filename)] #Read lines of output files and strip \n

    a = re.search('u_(.+?)-t', filestr) #Searches for given pattern
    if a:
        uf = int(a.group(1)) #Saves pattern as variable

    b = re.search('n_(.+?)-u', filestr)
    if b:
        grid = int(b.group(1))

    c = re.search('c_(.+?)-i', filestr)
    if c:
        comp = float(c.group(1))

    d = re.search('e_(.+?).out', filestr)
    if d:
        exp = float(d.group(1))

    e = re.search('i_(.+?)-e', filestr)
    if e:
        ins = float(e.group(1))
    
    for i in range(0,len(content)):
        if content[i].startswith('CURRENT_TIME'):

            sumf = 0
            sumEI = 0

            for j in range(i+3, i+3+grid*grid):
                temp = content[j].split()
                fitarray.append(float(temp[2]))
                EIarray.append(float(temp[3]))
                sumf = sumf + float(temp[2])
                sumEI = sumEI + float(temp[3])
                
            runwords = content[i].split() #Splits line into array based on space
            run = int(runwords[1])

            entwords = content[i+1].split()
            entropy = float(entwords[3])

            fitness = float(sumf/len(fitarray))
            EI = float(sumEI/len(EIarray))
                        
            fitarray.clear()
            EIarray.clear()
            
            fin.write(str(run) + "," + str(rep) + "," + str(uf) + "," + str(exp) + "," + str(comp) + "," + str(num) + "," + str(grid) + "," + str(ins) + "," + str(EI) + "," + str(entropy) + "," + str(fitness) + "\n")

    cur_file.close()        
    
    
    
