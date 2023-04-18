#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:00:05 2020

@author: qiyaozhu
"""

import os
from ClassesFunctions import *
from dualGraphs import *
import numpy as np
import sys
import numpy as np


Motifs = {}


def helixOrder(RNA):
    dic = {}
    for i in range(1,len(RNA.Helices)):
        dic[RNA.Helices[i].start] = i
        dic[RNA.Bases[RNA.Helices[i].end].indexBP] = i
    order = sorted(dic.keys())
    seq = []
    for o in order:
        seq.append(str(dic[o]))
    return seq


def ctToSequence(arg, ILmax):
    RNA = getCTInfo(arg)
    countHelices(RNA)
    changeHelices(RNA, ILmax)
    order = helixOrder(RNA)
    return RNA, ''.join(order)


def graphFinder(seq, DB, ILmax):
    
    with open("tmpRNAfold.db", 'w') as f:
        f.write(">84 MT246482.1 0 13405 84\n")
        f.write(seq + "\n")
        f.write(DB)
    os.system("dot2ct tmpRNAfold.db tmpRNAfold.ct")    
    
    RNA, VO = ctToSequence("tmpRNAfold.ct", ILmax)
    os.system("rm -rf tmpRNAfold.db tmpRNAfold.ct")
    
    RNA.makeMatrices()
    connectHelices(RNA)
    for i in range(0,len(RNA.adjMatrix)):
        vertexOrder.append(0)
        
    success, graph = calcEigen(RNA)
    correctHNumbers(RNA)
    if success:
        return graph, VO
    elif len(RNA.adjMatrix)==0 or len(RNA.adjMatrix)>9:
        print ("No matching graph exists because vertex number is either 0 or greater than 10.")
        return 'None', 'None'
    else: # no graph ID was assigned as eigen values not in the library S.J. 11/09/2017
        print ("No matching graph exists (even if the vertex number is between 2 and 9).")
        return 'None', 'None'
    

# This function calculates e^-beta*G    
def energyToProb(e):
   
    kB = 1.38*10**(-23)
    T = 310.15
    mol = 6.022*10**23
    kcal = 4184
    
    return np.exp(-np.array(e)*kcal/kB/T/mol)
    

# Calculates probability, given a vector of Gibbs free energy
def probCal(E):
    
    p = energyToProb(E)
    partition = np.sum(p)
    
    return p/partition


# For a sequence of N nucleotides with shape reactivity, calculate probabilities for all possible graphs
def graphProb(seq, file, ILmax):
    
    with open(file, 'r') as f:
        lines = f.readlines()
    
    energies = []
    DBs = [] 
    i = 0
    while i < len(lines):
        if lines[i] == str(file.split('/')[1].split('nt')[0])+'\n':
            energies.append(float(lines[i+1].split('\n')[0]))
            DBs.append(lines[i+2].split('\n')[0])
        i += 1
        
    graphs = []
    for db in DBs:
        graph, VO = graphFinder(seq, db, ILmax)
        graphs.append(graph+' ('+VO+')')
    probs = probCal(energies)
    
    for i in range(len(DBs)):
        print(graphs[i])
        print(DBs[i])
        print(probs[i])
        
    dic = {}
    for i in range(len(graphs)):
        g = graphs[i]
        if g not in dic:
            dic[g] = probs[i]
        else:
            dic[g] += probs[i]
    dic = {k: v for k, v in sorted(dic.items(), key=lambda item: item[1], reverse = True)}
    
    for i in range(len(graphs)):
        g = graphs[i]
        if dic[g]*100 >= 1:
            if g not in Motifs:
                Motifs[g] = [seq, DBs[i]]
    
    return dic
        


# python ProbCalc_Length.py FASTA ILmax
# export DATAPATH=/Users/qz886/Downloads/RNAstructure/data_tables/
if __name__== "__main__":    
    
    virus = sys.argv[1]
    ILmax = int(sys.argv[2])
    N = 30 # N upstream and downstream
    
    # add 7-nt slippery site first, then N upstream and downstream residues one-by-one on both sides
    for n in range(8):
        file = virus + '/' + str(77+n) + "nt.subopt"
        with open(file, 'r') as f:
            lines = f.readlines()
        for l in lines:
            if "Sequence:" in l:
                subseq = l.split()[2].split('\n')[0]
                continue
        dic = graphProb(subseq, file, ILmax)
        
        with open(virus+'.prob', 'a+') as f:
            f.write('>' + str(77+n) + 'nt\n')
            for g in dic:
                p = dic[g]*100
                if p >= 1:
                    s = '%.2f' %p
                    f.write(g+':\t'+s+'%\n')
            f.write('\n')
            
    for n in range(1,N+1):    
        file = virus + '/' + str(84+2*n) + "nt.subopt"
        with open(file, 'r') as f:
            lines = f.readlines()
        for l in lines:
            if "Sequence:" in l:
                subseq = l.split()[2].split('\n')[0]
                continue
        dic = graphProb(subseq, file, ILmax)
        
        with open(virus+'.prob', 'a+') as f:
            f.write('>' + str(84+2*n) + 'nt\n')
            for g in dic:
                p = dic[g]*100
                if p >= 1:
                    s = '%.2f' %p
                    f.write(g+':\t'+s+'%\n')
            f.write('\n')
                    
            
    Motifs = {k: v for k, v in sorted(Motifs.items(), key=lambda item: item[0])} # sort by graph ID    
    
    with open(virus+'_Motifs.txt', 'w') as f:
        for VO in Motifs:
            [subseq, db] = Motifs[VO]
            f.write(VO+'\n')
            
            # rewrite the DB structures in Motifs
            with open("tmpRNAfold.db", 'w') as f1:
                f1.write(">84 MT246482.1 0 13405 84\n")
                f1.write(subseq + "\n")
                f1.write(db)
            os.system("dot2ct tmpRNAfold.db tmpRNAfold.ct")  
            os.system("ct2dot tmpRNAfold.ct -1 tmpRNArefold.db")
            with open("tmpRNArefold.db", 'r') as f1:
                db = f1.readlines()[2].split('\n')[0]
            os.system("rm -rf tmpRNAfold.db tmpRNAfold.ct tmpRNArefold.db")
            
            # align all sequences
            if len(subseq) <= 84:
                emptyN = 114-len(subseq)
            else:
                emptyN = (144-len(subseq))//2                
            for e in range(emptyN):
                subseq = ' '+subseq
                db = ' '+db
            f.write(subseq+'\n'+db+'\n\n')

