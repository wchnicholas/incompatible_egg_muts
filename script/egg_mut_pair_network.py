#!/usr/bin/python
import os
import sys
import operator
import colorsys
import networkx as nx
import numpy as np
from math import exp, log
from itertools import imap
from collections import defaultdict

def CountFreq(single_mut_file):
  single_mut_dict = defaultdict(int)
  infile = open(single_mut_file)
  for line in infile.xreadlines():
    if 'count' in line: continue
    line = line.rstrip().rsplit("\t")
    mut    = line[0]
    count  = line[1]
    single_mut_dict[mut] = float(count)
  infile.close()
  return single_mut_dict

def CountingPairs(double_mut_file):
  double_mut_dict = {}
  infile = open(double_mut_file,'r')
  for line in infile.readlines():
    if 'count' in line: continue
    line = line.rstrip().rsplit("\t")
    pair  = line[0]
    count = float(line[1])
    mut1 = pair.rsplit('-')[0]
    mut2 = pair.rsplit('-')[1]
    aa1  = mut1[0]
    aa2  = mut2[0]
    pos1 = mut1[1::]
    pos2 = mut2[1::]
    if '-' in mut1 or '-' in mut2: continue 
    if pos1 == '0' or pos2 == '0': continue 
    if int(pos1) < 117 or int(pos1) > 265: continue 
    if int(pos2) < 117 or int(pos2) > 265: continue 
    if double_mut_dict.has_key(pair):
      double_mut_dict[pair] += count
    else:
      double_mut_dict[pair] = count
  return double_mut_dict

def FilteringPairs(double_mut_dict,Count_cutoff):
  new_double_mut_dict = {}
  for pair in double_mut_dict.keys():
    count = double_mut_dict[pair]
    if count >= Count_cutoff:
      new_double_mut_dict[pair] = count
  return new_double_mut_dict

def buildgraph(double_mut_dict):
  G=nx.Graph()
  for pair in double_mut_dict: 
    mut1 = pair.rsplit('-')[0]
    mut2 = pair.rsplit('-')[1]
    count = double_mut_dict[pair]
    G.add_edge(mut1,mut2,weight=count)
  return G

def labelnode(deg):
  high = float(222)
  low  = float(15)
  mid  = float(high-low)/2+low
  if deg > high:  return colorsys.rgb_to_hsv(1, 0, 0)
  elif deg > mid: return colorsys.rgb_to_hsv(1,(high-deg)/(high-mid),0)
  elif deg > low: return colorsys.rgb_to_hsv(1,1,(mid-deg)/(mid-low))
  elif deg <= low: return colorsys.rgb_to_hsv(1,1,1)
  else: print 'Something is wrong with the coloring function'; print deg; sys.exit()

def drawgraph(G,outfile,single_mut_dict,edge_rescale,node_rescale):
  outfile = open(outfile,'w')
  outfile.write('strict graph{'+"\n"+"\t"+'rankdir=LR'+
                                "\n"+"\t"+'node [shape=oval,fontname="arial bold",fontsize=18,margin=0]'+
                                "\n"+"\t"+'edge [arrowhead="none"]'+
                                "\n"+"\t"+'overlap = false'+
                                "\n"+"\t"+'splines = true'+"\n")
  for var in G.nodes():
    deg  =  G.degree(var,weight='weight')
    freq = float(single_mut_dict[var])
    size = str(int(log(deg)*node_rescale))
    size = str(0.8)
    col = labelnode(freq)
    col = ','.join(map(str,list(col)))
    outfile.write("\t"+var+' [fillcolor="'+col+'", color=black, fixedsize=shape, \
                              style="filled", width='+size+'];'+"\n")
  for E in G.edges():
    var1 = E[0]
    var2 = E[1]
    weight = float(G[var1][var2]['weight'])*edge_rescale
    outfile.write("\t"+var1+' -- '+var2+' [style=solid, color=black, penwidth='+str('{0:.20f}'.format(abs(weight)))+'];'+"\n")
  outfile.write('}'+"\n")
  outfile.close()

def main():
  single_mut_file  = 'table/mut_single_freq.tsv'
  double_mut_file  = 'table/mut_double_freq.tsv'
  outfile   = 'graph/EggMutPair_network.dot'
  Count_cutoff = 2
  edge_rescale = 0.25
  node_rescale = 1.5
  single_mut_dict = CountFreq(single_mut_file)
  double_mut_dict = CountingPairs(double_mut_file) 
  double_mut_dict = FilteringPairs(double_mut_dict,Count_cutoff)
  G = buildgraph(double_mut_dict)
  print "max count: %i" % max([single_mut_dict[node] for node in G.nodes()])
  drawgraph(G,outfile,single_mut_dict,edge_rescale,node_rescale)
  os.system('fdp -Tpng -Gdpi=300 %s -o %s' % (outfile,outfile.replace('.dot','.png')))

if __name__ == '__main__':
  main()
