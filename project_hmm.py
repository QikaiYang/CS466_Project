import numpy as np
import math
import dendropy
from math import log
dic = {"A":0, "T":1, "C":2, "G":3}
hmm_t = np.array([[log(6),log(2),log(2)],
                  [log(4),log(6),0],
                  [log(4),0,log(6)]])
hmm_x = np.array([log(2),log(2),log(2),log(2)])
hmm_y = np.array([log(2),log(2),log(2),log(2)])
hmm_m = np.array([[log(7),log(1),log(1),log(1)],
                  [log(1),log(7),log(1),log(1)],
                  [log(1),log(1),log(7),log(1)],
                  [log(1),log(1),log(1),log(7)]])

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
            spot = spot + (line.replace('\n','')).replace('-','')
    ls.append(spot)
    f.close()
    return ls

def compare_seq(seq1,seq2,hmm_x,hmm_y,hmm_m,hmm_t):
    score_m = np.zeros((len(seq1)+1, len(seq2)+1))
    score_x = np.zeros((len(seq1)+1, len(seq2)+1))
    score_y = np.zeros((len(seq1)+1, len(seq2)+1))
    for i in range(1,len(score_m)):
        score_m[i][0] = score_m[i-1][0]+hmm_x[dic[seq1[i-1]]]+hmm_t[1][1]
        score_x[i][0] = score_x[i-1][0]+hmm_x[dic[seq1[i-1]]]+hmm_t[1][1]
    for j in range(1,len(score_m[0])):
        score_m[0][j] = score_m[0][j-1]+hmm_y[dic[seq2[j-1]]]+hmm_t[2][2]
        score_y[0][j] = score_y[0][j-1]+hmm_y[dic[seq2[j-1]]]+hmm_t[2][2]
    for i in range(1, len(score_m)):
        for j in range(1, len(score_m[0])):
            score_m[i][j] = hmm_m[dic[seq1[i-1]]][dic[seq2[j-1]]] + max(score_m[i-1][j-1]+hmm_t[0][0],score_x[i-1][j-1]+hmm_t[1][0],score_y[i-1][j-1]+hmm_t[1][0])
            score_x[i][j] = hmm_x[dic[seq1[i-1]]] + max(hmm_t[0][1]+score_m[i-1][j], hmm_t[1][1]+score_x[i-1][j])
            score_y[i][j] = hmm_y[dic[seq2[j-2]]] + max(hmm_t[0][1]+score_m[i][j-1], hmm_t[1][1]+score_y[i][j-1])
    return (max(score_m[len(score_m)-1][len(score_m[0])-1], score_x[len(score_x)-1][len(score_x[0])-1], score_y[len(score_y)-1][len(score_y[0])-1]))
for i in range(5):
    locals()["seq"+str(i)] = read_seq("./100S1/R"+str(i)+"/rose.aln.true.fasta")
distance_m = np.zeros((5,100, 100))

for k in range(5):
    print("------")
    for i in range(len(locals()["seq"+str(k)])):
        for j in range(i+1,len(locals()["seq"+str(k)])):
            distance_m[k][i][j] = compare_seq(locals()["seq"+str(k)][i],locals()["seq"+str(k)][j],hmm_x,hmm_y,hmm_m,hmm_t)
            distance_m[k][j][i] = distance_m[k][i][j]
            print(i,j,distance_m[k][i][j])

for i in range(5):
    np.savetxt("./100S1/R"+str(i)+"/distance_hmm.txt",distance_m[i])