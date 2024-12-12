import os
from multiprocessing import Pool
cmds=[]
with open('d.sh','r') as f:
    for line in f:
        a=line.split('\n')[0]
        b=a.split(';')
        if(len(b)>1):
            c1=b[0]
            c2=b[1]
            cmds.append(c1)
            cmds.append(c2)
        else:
            c3=b[0]
            cmds.append(c3)

def f(cmd):
    print(cmd)
    os.system(cmd)



p = Pool(len(cmds))
p.map(f,cmds)

