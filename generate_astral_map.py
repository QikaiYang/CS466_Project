import dendropy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Trees", type=str, help="Input Fasta File Name") #args.Trees is the input fasta files
parser.add_argument("Map", type=str, help="Output Map File Name")    #args.Map is the output Map file
args = parser.parse_args()

dic = {} #used to record how many species and how many multi-copy genes for each species 

def tackle(strr): # tackle each string to get rid of the "_"
    cut_loaction = 0
    for i in range(len(strr)):
        if(strr[i]=="_" or strr[i]==" "):
            cut_loaction = i
            break
    return strr[0:(cut_loaction)]

f = open(args.Trees, 'r')
for line in f.readlines():                #get all the species
    temp = "[&R]" + (line.replace("\n","")).replace(" ", "")
    tree = (dendropy.Tree.get(data=temp, schema="newick", rooting='force-rooted'))
    leafs = tree.leaf_nodes()
    for leaf in leafs:
        if(dic.has_key(tackle(leaf.taxon.label)) != True):
            dic[tackle(leaf.taxon.label)] = []

f = open(args.Trees, 'r')
for line in f.readlines():                #get all the multi-copy genes
    temp = "[&R]" + (line.replace("\n","")).replace(" ","")
    tree = (dendropy.Tree.get(data=temp, schema="newick", rooting='force-rooted'))
    leafs = tree.leaf_nodes()
    for leaf in leafs:
        dic[tackle(leaf.taxon.label)].append(leaf.taxon.label)

#get rid of the duplicated multi-copy genes
for item in dic:
    dic[item] = list(set(dic[item]))

#write to the output map file
result = ""
for key in dic:
    result += key + " : "
    for i in range(len(dic[key])):
        dic[key][i] = dic[key][i].replace(" ","_")
        result += dic[key][i] + ", "
    result = result[:len(result)-2]
    result += "\n"
print(result)

f = open(args.Map, 'w')
f.write(result)
f.close
