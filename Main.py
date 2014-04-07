'''
Created on Apr 6, 2014

@author: Thibault
'''
import os
from AbundanceTree import *

os.chdir("C:\\Users\\Thibault\\Desktop\\LBE\\INRA_M1_PythonScript")

N=10**8
tree=AbundanceTree("silva.nh")
s=len(tree)
tree.power_law(N, True)
tree.fill_tree()

diversity=[]
for i in (0,1,2):
    diversity.append(tree.qdiversity(i))

tree.sample(10**4)
goodturing=tree.dico_good_turing()

minmin=[]
for i in (0,1,2):
    btree=tree.copy()
    minmin.append(btree.minminestimate(goodturing, i))

maxmax=[]
for i in (0,1,2):
    btree=tree.copy()
    maxmax.append(btree.maxmaxestimate(N,goodturing, i))

minmax=[]
for i in (0,1,2):
    btree=tree.copy()
    minmax.append(btree.minmaxestimate(N,goodturing, i))

maxmin=[]
for i in (0,1,2):
    btree=tree.copy()
    maxmin.append(btree.maxminestimate(goodturing, i))

for q in (0,1,2):
    print("for q="+str(q))
    print("\t The underlying diversity is "+str(diversity[q][0]))
    print("\t The MINMIN estimate is "+str(minmin[q][0]))
    print("\t The MAXMIN estimate is "+str(maxmin[q][0]))
    print("\t The MINMAX estimate is "+str(minmax[q][0]))
    print("\t The MAXMAX estimate is "+str(maxmax[q][0]))

input("Press RETURN to close...")