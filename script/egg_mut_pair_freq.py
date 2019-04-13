#!/usr/bin/python
import os
import re
import sys
from Bio import SeqIO
from collections import Counter, defaultdict

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def reading_mut_list(infile, count_cutoff):
  infile = open(infile,'r')
  egg_mut_list = []
  for line in infile.xreadlines():
    if 'count' in line: continue
    line = line.rstrip().rsplit("\t")
    mut   = line[0]
    count = line[1]
    if int(count) >= count_cutoff:
      egg_mut_list.append(mut)
  infile.close()
  return egg_mut_list

def extracting_egg_seq(records):
  egg_seq_list = []
  excludepattern = re.compile ("UNKNOWN_1_RHMK|TMK1_MDCK|AMNIOTIC_1_PRHMK_2|M1_RII1,C5|R1_C|R1_S|RII1_C|RII1_S|RIIX_C|RX_C|MDCK_1_RHMK|NC|_MK1")
  unpassagedpattern = re.compile("LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT")
  eggpattern = re.compile("AM[1-9]|E[1-7]|AMNIOTIC|EGG|EX|AM_[1-9]")
  cellpattern = re.compile("S[1-9]|SX|SIAT|MDCK|C[1-9]|CX|C_[1-9]|M[1-9]|MX|X[1-9]|^X_$")
  siatpattern = re.compile("^S[1-9]_$|SIAT2_SIAT1|SIAT3_SIAT1")
  monkeypattern=re.compile("TMK|RMK|RHMK|RII|PMK|R[1-9]|RX")
  siatexcludepattern=re.compile("SIAT|SX|S[1-9]")
  for record in records:
    header = str(record.id)
    if header.count('|')!=4: continue
    ID   = header.rsplit('|')[0]
    PSG  = header.rsplit('|')[1].upper()
    year = header.rsplit('|')[-1][1:5]
    seq  = str(record.seq)
    if ID=='A/Aichi/2/1968_': refseq = seq
    assert(isInt(year))
    if 'X' in seq: continue
    if excludepattern.search(PSG): continue
    elif eggpattern.search(PSG):
      egg_seq_list.append(seq)
  return egg_seq_list

def counting_pair(double_out_file, egg_seq_list, egg_mut_list, posoffset):
  double_count_dict = defaultdict(int)
  single_count_dict = defaultdict(int)
  for seq in egg_seq_list:
    for i in range(len(egg_mut_list)):
      mut_i     = egg_mut_list[i]
      adj_pos_i = int(mut_i[1::])-posoffset-1
      aa_i      = mut_i[0]     
      if seq[adj_pos_i] == aa_i:
        single_count_dict[mut_i] += 1
      for j in range(len(egg_mut_list)):
        if i < j: 
          mut_j     = egg_mut_list[j]
          adj_pos_j = int(mut_j[1::])-posoffset-1
          aa_j      = mut_j[0]     
          pair_ID = mut_i+'-'+mut_j
          if seq[adj_pos_i] == aa_i and seq[adj_pos_j] == aa_j:
            double_count_dict[pair_ID] += 1
  return single_count_dict, double_count_dict

def writing_double_count(double_count_dict, double_out_file):
  print "writing %s" % double_out_file
  outfile = open(double_out_file, 'w')
  outfile.write('double'+"\t"+'count'+"\n")
  for double in sorted(double_count_dict.keys(), key=lambda x:double_count_dict[x], reverse=True):
    outfile.write(double+"\t"+str(double_count_dict[double])+"\n")
  outfile.close()

def writing_single_count(single_count_dict, single_out_file):
  print "writing %s" % single_out_file
  outfile = open(single_out_file, 'w')
  outfile.write('single'+"\t"+'count'+"\n")
  for single in sorted(single_count_dict.keys(), key=lambda x:single_count_dict[x], reverse=True):
    outfile.write(single+"\t"+str(single_count_dict[single])+"\n")
  outfile.close()

def main():
  egg_mut_list_file = 'table/mut_count_year_egg.tsv'
  double_out_file     = 'table/mut_double_freq.tsv'
  single_out_file     = 'table/mut_single_freq.tsv'
  seq_aln_file      = 'Fasta/HumanH3N2_All_2018.aln'
  count_cutoff = 5
  posoffset    = -16
  records         = [record for record in SeqIO.parse(seq_aln_file,"fasta")]
  egg_seq_list    = extracting_egg_seq(records)
  egg_mut_list    = reading_mut_list(egg_mut_list_file, count_cutoff)
  single_count_dict, double_count_dict = counting_pair(double_out_file, egg_seq_list, egg_mut_list, posoffset)
  writing_double_count(double_count_dict, double_out_file)
  writing_single_count(single_count_dict, single_out_file)

if __name__ == "__main__":
  main()
