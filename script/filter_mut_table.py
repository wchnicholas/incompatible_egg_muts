#!/usr/bin/python
import os
import re
import sys
import glob
from collections import Counter, defaultdict

def ReadingEggMutFreq(mutfreqfile,years):
  EggMutdict = defaultdict(dict)
  infile = open(mutfreqfile,'r')
  for line in infile.xreadlines():
    if 'year' in line: continue
    line = line.rstrip().rsplit("\t")
    year = line[0]
    aa   = line[2]
    eggcount = float(line[3])
    oricount = float(line[4])
    eggfreq = float(line[5])
    orifreq = float(line[6])
    pos  = line[7]
    mut  = aa+pos
    if int(year) in years:
      if orifreq == 0:
        EggMutdict[year][mut] = eggfreq
  infile.close()
  return EggMutdict

def FilterEggMut(EggMutdict, minOccur):
  muts = [mut for year in EggMutdict.keys() for mut in EggMutdict[year].keys()]
  mutOccur = Counter(muts)
  mutofinterest = [mut for mut in mutOccur if mutOccur[mut] >= minOccur]
  return mutofinterest, mutOccur

def WritingMut(Mutdict,mutofinterest,Eggoutfile):
  outfile = open(Eggoutfile,'w')
  outfile.write("\t".join(['year','pos','mut','freq'])+"\n")
  for mut in mutofinterest:
    pos = mut[1::] 
    for year in sorted(Mutdict.keys()):
      freq = Mutdict[year][mut] if Mutdict[year].has_key(mut) else 0
      outfile.write("\t".join(map(str,[year, pos, mut, freq]))+"\n")
  outfile.close()

def WritingEggMut(mutOccur_dict, Egg_mutlist_file):
  outfile = open(Egg_mutlist_file, 'w')
  outfile.write('eggmut'+"\t"+'count'+"\n")
  for mut in sorted(mutOccur_dict.keys(), key=lambda x:int(mutOccur_dict[x]), reverse=True):
    outfile.write(mut+"\t"+str(mutOccur_dict[mut])+"\n")
  outfile.close()

def main():
  mutfreqfile = 'table/egg_vs_ori_bases.tsv'
  Eggoutfile     = 'table/mut_freq_egg_filtered.tsv'
  Egg_mutlist_file  = 'table/mut_count_year_egg.tsv'
  years = range(2003,2019)
  minOccur  = 5
  assert(len(years)>minOccur)
  print "Analyzing %i years: %s" % (len(years), ' '.join(map(str,years)))
  EggMutdict = ReadingEggMutFreq(mutfreqfile,years)
  Eggmutofinterest, mutOccur_dict = FilterEggMut(EggMutdict,minOccur)
  print "Mutations of interest: %s" % ' '.join(sorted(Eggmutofinterest,key=lambda x:int(x[1::])))
  WritingEggMut(mutOccur_dict, Egg_mutlist_file)
  WritingMut(EggMutdict,Eggmutofinterest,Eggoutfile)
  
if __name__ == "__main__":
  main()
