import numpy as np
import math
import dendropy
from math import log
dic = {"A":0, "T":1, "C":2, "G":3}

def read_name(file):
    f = open(file)
    ls=[]
    for line in f:
        if (line.startswith('>')):
            ls.append((line.replace('\n','')).replace('>',''))
    f.close()
    return ls

def read_seq(file):
    f = open(file)
    ls=[]
    spot = ""
    for line in f:
        if (line.startswith('>')):
            if(len(spot)!=0):
                ls.append(spot)
                spot = ""
        else:
            spot = spot + line.replace('\n','')
    ls.append(spot)
    f.close()
    return ls

def compare_seq(seq1,seq2):
    haha = min(len(seq1),len(seq2))
    result = 0
    for i in range(haha):
        if(seq1[i]!=seq2[i]):
            result += 1.0
    return (-0.75*log(1-3.0/4.0*result/haha))

for i in range(5):
    locals()["seq"+str(i)] = read_seq("./100S1/R"+str(i)+"/rose.aln.true.fasta")
distance_m = np.zeros((5, 100, 100))

for k in range(5):
    print("------")
    for i in range(len(locals()["seq"+str(k)])):
        for j in range(i+1,len(locals()["seq"+str(k)])):
            distance_m[k][i][j] = compare_seq(locals()["seq"+str(k)][i],locals()["seq"+str(k)][j])
            distance_m[k][j][i] = distance_m[k][i][j]
            print(i,j,distance_m[k][i][j])

for i in range(5):
    np.savetxt("./100S1/R"+str(i)+"/distance_JC.txt",distance_m[i])