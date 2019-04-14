#!/usr/bin/python
import os
import sys
import glob
import string
import operator
from Bio import SeqIO
from itertools import imap
from collections import Counter, defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessEgglib(R1file, refseq):
  print "Reading %s" % R1file
  R2file = R1file.replace('_R1','_R2')
  indexlength  = 3
  Primerlength = 23
  R1frameshift = 0
  R2frameshift = 0
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  variant_dict = defaultdict(int)
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    index  = R1seq[0:indexlength]
    R1roi = R1seq[indexlength+Primerlength:-2]
    R2roi = R2seq[Primerlength:-2]
    if 'N' in R1roi or 'N' in R2roi: continue
    R1pep = translation(R1roi)
    R2pep = translation(rc(R2roi))
    if R1pep[57::] == R2pep[0:34]:
      pep = R1pep[0:57]+R2pep
      if '_' in pep: continue
      mut = Bris07mut2ID(pep, refseq)
      variant_dict[mut] += 1
  return variant_dict

def Bris07mut2ID(mut,refseq):
  shift = 117
  haplo = []
  assert(len(mut)==len(refseq))
  for n in range(len(mut)):
    pos = n+shift
    if refseq[n]!=mut[n]:
       haplo.append(refseq[n]+str(pos)+mut[n])
  return '-'.join(haplo)

def extract_count(mut, count_dict):
  return 0 if mut not in count_dict.keys() else count_dict[mut]

def Output(count_dict, passages, outfile):
  print "Compiling results into %s" % outfile
  outfile = open(outfile,'w')
  muts = list(set([mut for passage in count_dict.keys() for mut in count_dict[passage].keys()]))
  outfile.write("\t".join(['mut','P0','P1','P2','P3','P4','P5'])+"\n")
  for mut in muts:
    ID  = 'WT' if mut=='' else mut
    out = [ID]
    for passage in passages:
      count    = 0 if passage not in count_dict.keys() else extract_count(mut, count_dict[passage])
      out.append(count)
    outfile.write("\t".join(map(str,out))+"\n")
  outfile.close()

def main():
  file_refseq = 'Fasta/Bris07_RBD.pep'
  filenames   = glob.glob('fastq/*R1*.fastq')
  refseq      = open(file_refseq,'r').readlines()[1].rstrip()
  experiments = ['G186V_rep1','G186V_rep2','G186V_rep3','L194P_rep1']
  passages    = ['P0','P1','P2','P3','P4','P5']
  for experiment in experiments: 
    outfile  = 'result/Bris07_Egg_'+experiment+'_count.tsv'
    count_dict  = {}
    for passage in passages:
      filename = 'fastq/'+experiment+'_'+passage+'_R1.fastq'
      count_dict[passage] = ProcessEgglib(filename,refseq)
    Output(count_dict, passages, outfile)

if __name__ == "__main__":
  main()
