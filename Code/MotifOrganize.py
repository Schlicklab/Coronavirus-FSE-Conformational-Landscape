#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:37:15 2020

@author: qiyaozhu
"""
from copy import *
from dualGraphs import *
from ClassesFunctions import *
import sys
import os.path
import numpy as np
import pandas as pd
from plotnine import *
from igraph import *


# Find the minimal motif interacting with the sequence covered in low<=residue (starting from 1)<=up
def minimal_motif(arg, low, up, ILmax):
    
    RNA = getCTInfo(arg)
    countHelices(RNA)
    changeHelices(RNA, ILmax)    
    
    dic = {}
    for i in range(1,len(RNA.Helices)):
        dic[RNA.Helices[i].start] = i
        dic[RNA.Bases[RNA.Helices[i].end].indexBP] = i
    order = sorted(dic.keys())
    HO = [dic[o] for o in order]
        
    IH = [] # list of helices involving the sequence
    
    for i in range(1,len(RNA.Helices)):
        start5 = RNA.Helices[i].start
        end5 = RNA.Helices[i].end
        start3 = RNA.Bases[RNA.Helices[i].end].indexBP
        end3 = RNA.Bases[RNA.Helices[i].start].indexBP
        
        seq = set(range(low, up+1))
        strand5 = set(range(start5, end5+1))
        strand3 = set(range(start3, end3+1))
        
        if bool(strand5 & seq) or bool(strand3 & seq):
            IH.append(i)

    LH = [i for i, x in enumerate(HO) if x in IH][0]
    UH = [i for i, x in enumerate(HO) if x in IH][-1]
    S = []
    S_new = HO[LH : UH+1]   
    
    while S != S_new:
        S = S_new
        LH = [i for i, x in enumerate(HO) if x in S][0]
        UH = [i for i, x in enumerate(HO) if x in S][-1]
        S_new = HO[LH : UH+1]

    motif = np.array(HO[LH : UH+1]) - min(HO[LH : UH+1]) + 1 # the minimal motif with helices renumbered from 1
    MVO = [str(s) for s in motif]
    
    return ''.join(MVO)



def unique_motifs(Motifs, ILmax):  

    with open(Motifs, 'r') as f:
        lines = f.readlines()
    
    unique_Motifs = {}
    
    i = 0
    while i < len(lines):
        if lines[i] != '\n':
            graph = lines[i].split('\n')[0]
            seq = lines[i+1].split()[0]
            db = lines[i+2].split()[0]
            
            if len(seq) <= 84:
                low = len(seq)-77+1
                up = len(seq)
            else:
                low = (len(seq)-84)//2+8
                up = len(seq) - (len(seq)-84)//2
            
            with open("tmpRNAfold.db", 'w') as f:
                f.write(">84 MT246482.1 0 13405 84\n")
                f.write(seq + "\n")
                f.write(db)
            os.system("dot2ct tmpRNAfold.db tmpRNAfold.ct")    
            
            MVO = minimal_motif("tmpRNAfold.ct", low, up, ILmax)
            os.system("rm -rf tmpRNAfold.db tmpRNAfold.ct")
            
            if MVO in unique_Motifs:
                unique_Motifs[MVO].append(graph)
            else:
                unique_Motifs[MVO] = [graph]
            
            i += 3
        i += 1
            
    return unique_Motifs

    
    
def rewrite_Prob(Motifs, Prob, ILmax):
    
    unique_Motifs = unique_motifs(Motifs, ILmax)
    print(unique_Motifs)
    motif_list = sorted(unique_Motifs)
    probMat = np.zeros([38, len(motif_list)])
    
    with open(Prob, 'r') as f:
        lines = f.readlines()
        
    Lengths = []
        
    i = 0
    k = 0
    while i < len(lines):
        if list(lines[i])[0] == '>':
            Lengths.append(lines[i].split('>')[1].split('nt')[0])
            i += 1
            while list(lines[i])[0] != '\n':
                j = 0
                Found = False
                while j < len(motif_list) and not Found:
                    if lines[i].split(':')[0] in unique_Motifs[motif_list[j]]:
                        probMat[k,j] += float(lines[i].split('\t')[1].split('%')[0])
                        Found = True
                    j += 1                          
                i += 1                
            k += 1
        i += 1
    
    columns = [] # only write those motif columns that have at least one probability >=5%
    
    for j in range(len(motif_list)):
        significance = False
        i = 0
        while i<38 and not significance:
            if probMat[i,j] >= 5:
                significance = True
            i += 1
        if significance:
            columns.append(j)
            
    Motif_list = [motif_list[k] for k in columns]
    ProbMat = probMat[:,columns]
    
    # rank the motifs
    weights = np.zeros([1, len(Motif_list)])
    for j in range(len(Motif_list)):
        weights[0,j] = np.sum(np.exp(2*ProbMat[:,j]/100)-1)
    rank = np.flip(np.argsort(weights))
    rank = list(rank[0])
    
    motif_list = [Motif_list[k] for k in rank]
    probMat = ProbMat[:,rank]
    Conformers = []
    
    with open('group.csv', 'w') as f:
        f.write('Length,')
        for j in range(len(motif_list)-1):
            m = motif_list[j]
            conf = unique_Motifs[m][0]
            Conformers.append(conf)
            f.write('Motif %d: ' %(j+1) + conf.split()[0] + ',')
        conf = unique_Motifs[motif_list[-1]][0]
        Conformers.append(conf)
        f.write('Motif %d: ' %len(motif_list) + conf.split()[0] + '\n')
        
        for le in range(38):
            f.write(Lengths[le]+',')
            for j in range(len(motif_list)-1):
                f.write('%.2f,' %probMat[le,j])
            f.write('%.2f\n' %probMat[le,-1])
            
    return Conformers



def plot_Landscape(group):
        
    probMat = pd.read_csv(group)    
    data = pd.melt(probMat, id_vars=['Length'], var_name='Conformations')
    g = ggplot(data, aes(x='Length', y='value', fill='Conformations', order='value')) + geom_bar(stat = "identity", position = position_stack(reverse = True)) #+ scale_fill_manual(values=['green', 'violet', 'red'])
    g = g + xlab("Upstream Residues") + ylab("Probability (%)") + ggtitle("Conformational Landscape") + theme_bw(base_family = "Arial", base_size=15) \
        #+ scale_x_continuous(breaks = (1,2,3,4,5,6,7,8,9), labels = [1,2,3,4,5,6,7,8,9])
    ggsave(plot = g, filename = 'Landscape.pdf', path = "/Users/qz886/Desktop/Coronavirus_Landscape/FS_Testing", width = 15, height = 15, units = "cm")



def plotDual(g):
    
    A = np.array(g.adjMatrix)
    g1 = Graph()
    n = A[0].size
    g1.add_vertices(n)
        
    mylayout=g1.layout_circle()
    bbox=BoundingBox(400,400)

    for i in range(0,n):    
        for j in range(i,n):   
            for k in range(0,int(A[i][j])):
                g1.add_edge(i,j)
    
    mylayout=g1.layout_fruchterman_reingold()
                
    # vlabel = np.argsort(FV)+1
    vlsize = 30
    vldist = 1.5
    
    plot(g1, "%s.eps"%g.graphID,vertex_color="#FF0000", vertex_frame_color="#FF0000", vertex_size=40, edge_color="#FF0000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
 


# export DATAPATH=/Users/qz886/Downloads/RNAstructure/data_tables/
if __name__== "__main__":
    
    Motifs = sys.argv[1]
    Prob = sys.argv[2]
    ILmax = int(sys.argv[3])
    Conformers = rewrite_Prob(Motifs, Prob, ILmax)
    # plot_Landscape('group.csv')
    
    
    IDs = []
    DBs = []
    
    for c in Conformers:
        IDs.append(c.split()[0])
        with open(Motifs, 'r') as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].split('\n')[0] == c:
                DBs.append([c, lines[i+1].split()[0], lines[i+2].split()[0]])
                break
            i += 4            
        
    # for graph in IDs:
    #     n = int(graph.split('_')[0])
    #     Graphs = []
    #     eigen_file = "DualEig/%dEigen"%n
    #     adj_file = "DualAdj/V%dAdjDG"%n
    #     loadEigenvalues(Graphs,n,eigen_file)
    #     loadAdjMatrices(Graphs,n,adj_file)        
    #     for g in Graphs:
    #         if g.graphID == graph:
    #             plotDual(g)
    
    
    with open('Conformers.txt', 'w') as f:
        for structure in DBs:
            [c, seq, db] = structure
            # align all sequences
            if len(seq) <= 84:
                emptyN = 114-len(seq)
            else:
                emptyN = (144-len(seq))//2                
            for e in range(emptyN):
                seq = ' '+seq
                db = ' '+db
            f.write('>'+c+'\n'+seq+'\n'+db+'\n\n')
    
    
    # with open('ArcPlot.R', 'w') as f:
    #     f.write('library(R4RNA)\n\n')
    # for c in Conformers:
    #     with open('ArcPlot.R', 'a+') as f:
    #         f.write('FSE <- readVienna("'+c+'.db")\n')
    #         f.write('FSE$col <- "black"\n')
    #         f.write('setEPS()\n')
    #         f.write('postscript("'+c+'.eps")\n')
    #         f.write('plotHelix(FSE, add = FALSE, shape = "circle", scale = FALSE, line = TRUE, lwd = 3, png = NA, pdf = NA)\n')
    #         f.write('dev.off()\n\n')
        
    
    
    
    
    
    
    
    
