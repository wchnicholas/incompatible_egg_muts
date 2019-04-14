#!/usr/bin/python
import os
import sys
import glob
import string
import operator
from Bio import SeqIO
from itertools import imap
from collections import Counter, defaultdict

def reading_count(file_count, count_dict, sample):
  infile = open(file_count,'r')
  count_dict[sample] = defaultdict(dict)
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    if 'mut' in line: continue
    mut        = line[0]
    count_P0  = line[1]
    count_P1   = line[2]
    count_P2   = line[3]
    count_P3   = line[4]
    count_P4   = line[5]
    count_P5   = line[6]
    count_dict[sample][mut] = {'P0':count_P0, 'P1':count_P1, 'P2':count_P2,
                               'P3':count_P3, 'P4':count_P4, 'P5':count_P5}
  infile.close()
  return count_dict

def count_to_freq(count_dict, outfile):
  freq_dict = {}
  total_count_dict = defaultdict(dict)
  samples   = count_dict.keys()
  muts = []
  [muts.extend(mut.rsplit('-')) for sample in samples for mut in count_dict[sample].keys()]
  muts = sorted(list(set(muts)),key=lambda x:int(x[1:-1]) if x != 'WT' else x)
  passages = ['P0','P1','P2','P3','P4','P5']
  print "writing: %s" % outfile
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['mut', 'passage', 'sample', 'mut_freq'])+"\n")
  for sample in samples:
    for passage in passages:
      total_count = sum([int(count_dict[sample][m1][passage]) for m1 in count_dict[sample].keys()])
      total_count_dict[sample][passage] = total_count
  for mut in muts:
    for sample in samples:
      for passage in passages:
        mut_count = 0
        for mutant in count_dict[sample].keys():
          if mut in mutant:
            mut_count += int(count_dict[sample][mutant][passage])
        mut_freq = float(mut_count)/float(total_count_dict[sample][passage]) if total_count_dict[sample][passage] != 0 else 'NA'
        outfile.write("\t".join(map(str,[mut, passage, sample, mut_freq]))+"\n")
  outfile.close()

def main():
  filenames   = glob.glob('result/*_count.tsv')
  outfile = "result/mut_freq.tsv" 
  count_dict = {}
  for filename in filenames:
    sample = filename.rsplit('result/Bris07_Egg_')[1].rsplit('_count.tsv')[0]
    count_dict = reading_count(filename, count_dict, sample)
  freq_dict  = count_to_freq(count_dict, outfile)

if __name__ == "__main__":
  main()
