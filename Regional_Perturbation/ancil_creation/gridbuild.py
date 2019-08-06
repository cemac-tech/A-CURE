#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 15:20:49 2019

@author: chmcsy
"""
"""

with open('./landsea_dump.dat','r') as fin:
    linelst=fin.readlines()
    
grid=[]

long=""
lat=""

for i in range (45,72):
    long=long+linelst[i].strip()

for i in range(96,1247,8):
    catline=""
    for j in range (0,8):
        catline = catline+linelst[i+j].strip()+" "
    grid.append(catline)

with open ('./grid.dat','w') as fout:
    for line in grid:
        fout.write("%s\n" % line)
"""

grid=[]

with open('./ResEmiss.csv','r') as fin:
    grid = fin.readlines()

fout=open('./Resreg.dat','w')
    
for line in grid:
    groups=line.split(",")
    fout.write("  "+", ".join(groups[:25])+",\n")
    fout.write("    "+", ".join(groups[25:49])+",\n")
    fout.write("    "+", ".join(groups[49:73])+",\n")
    fout.write("    "+", ".join(groups[73:97])+",\n")
    fout.write("    "+", ".join(groups[97:121])+",\n")
    fout.write("    "+", ".join(groups[121:145])+",\n")
    fout.write("    "+", ".join(groups[145:169])+",\n")
    fout.write("    "+", ".join(groups[169:192]).strip()+",\n")
    
fout.close()
    