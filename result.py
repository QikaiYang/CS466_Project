import numpy as np
import matplotlib.pyplot as plt
import Bio
import dendropy
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceMatrix 
import compare_trees as cmp
from Bio import Phylo
tax = dendropy.TaxonNamespace()
#-----------------------------------------------
start = "./100S1/R"
end = "/rose.true.fasta"
tree = "/rose.tt"
tree2 = "/tree.tt"
tax = dendropy.TaxonNamespace()
#-----------------------------------------------

def read_name(file):
    f = open(file)
    ls=[]
    for line in f:
        if(line.startswith('>')):
            ls.append((line.replace('\n','')).replace('>',''))
    f.close()
    return ls

def get_file(address):
    f = open(address)
    l = np.array([])
    for line in f:
        temp = line.split(' ')
        for j in range(len(temp)):
            l = np.append(l,[float(temp[j])])
    l = l.reshape((100,100))
    return l

def transfer(arr):
    result = []
    for i in range(len(arr)):
        temp = []
        for j in range(i+1):
            temp.append(arr[i][j])
        result.append(temp)
    return result

for i in range(5):
    locals()["real_tree"+str(i)] = dendropy.Tree.get(path = start+str(i)+tree,rooting='force-unrooted',taxon_namespace=tax, schema="newick")
    locals()["name"+str(i)] = read_name(start+str(i)+end)

result = np.zeros((3,5)) # FP, FN, ERROR #hmm, hamming, JC

for i in range(5):
    locals()["hmm"+str(i)] = get_file("./100S1/R"+str(i)+"/distance_hmm.txt")
    ###
    locals()["hmm"+str(i)][i][i] = locals()["hmm"+str(i)][1][0]
    max = np.max(locals()["hmm"+str(i)])
    min = np.min(locals()["hmm"+str(i)])
    for p in range(len(locals()["hmm"+str(i)])):
        for q in range(len(locals()["hmm"+str(i)][0])):
            if(p!=q):
                locals()["hmm"+str(i)][p][q] = (locals()["hmm"+str(i)][p][q]-min)/(max-min)
            else:
                locals()["hmm"+str(i)][p][q] = 0 
    ###    
    locals()["hamming"+str(i)] = get_file("./100S1/R"+str(i)+"/distance_hamming.txt") 
    locals()["JC"+str(i)] = get_file("./100S1/R"+str(i)+"/distance_JC.txt") 
    #--------------------------------------------------Hmm
    dm_hmm = DistanceMatrix(locals()["name"+str(i)], transfer(locals()["hmm"+str(i)]))
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm_hmm)
    Phylo.write(tree, "./100S1/R"+str(i)+"/tree.tt", 'newick')
    estimated_tree = dendropy.Tree.get(path = start+str(i)+tree2,rooting='force-unrooted',taxon_namespace=tax, schema="newick")
    print(cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))
    result[0][i] = (cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))[6]
    #--------------------------------------------------Hamming
    dm_hamming = DistanceMatrix(locals()["name"+str(i)], transfer(locals()["hamming"+str(i)]))
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm_hamming)
    Phylo.write(tree, "./100S1/R"+str(i)+"/tree.tt", 'newick')
    estimated_tree = dendropy.Tree.get(path = start+str(i)+tree2,rooting='force-unrooted',taxon_namespace=tax, schema="newick")
    print(cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))
    result[1][i] = (cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))[6]
    #--------------------------------------------------JC
    dm_jc = DistanceMatrix(locals()["name"+str(i)], transfer(locals()["JC"+str(i)]))
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm_jc)
    Phylo.write(tree, "./100S1/R"+str(i)+"/tree.tt", 'newick')
    estimated_tree = dendropy.Tree.get(path = start+str(i)+tree2,rooting='force-unrooted',taxon_namespace=tax, schema="newick")
    print(cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))
    result[2][i] = (cmp.compare_trees(locals()["real_tree"+str(i)], estimated_tree))[6]
    #--------------------------------------------------
name_list = ['R0','R1','R2','R3','R4']
label = ["HMM","Hamming","Jukes & Cantor"]
x =list(range(len(result[0])))
total_width, n = 0.8, 3
width = total_width / n
plt.figure(1)
plt.bar(x, result[0], width=width, label='HMM', fc = 'y')
for i in range(len(x)):
    x[i] = x[i] + width
plt.bar(x, result[1], width=width, label='Hamming',tick_label = name_list, fc = 'r')
for i in range(len(x)):
    x[i] = x[i] + width
plt.bar(x, result[2], width=width, label='Jukes & Cantor', fc = 'g')
plt.legend()
plt.figure(2)
plt.boxplot(result.T, labels = label, sym = "o")
plt.show()