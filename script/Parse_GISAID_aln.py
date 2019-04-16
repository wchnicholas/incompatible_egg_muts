#!/usr/bin/python
import os
import re
import sys
import glob
from scipy import stats
from Bio import SeqIO
from collections import Counter, defaultdict
import pandas as pd
import pickle
from itertools import combinations

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def ExtractEggPSGNumber(PSG):
  if '/' not in PSG and '+' not in PSG and ',' not in PSG: 
    PSG_number = PSG.replace('_','')
    if ':' in PSG_number: PSG_number = PSG_number.rsplit(':')[1]
    if 'E' in PSG_number and isInt(PSG_number[1::]):
      return PSG_number.replace('E','')
    elif 'AM' in PSG_number and isInt(PSG_number[2::]):
      return PSG_number.replace('AM','')

def ParseSeqs(records):
  refseq   = ''
  Egg_dict = defaultdict(list)
  Ori_dict = defaultdict(list)
  Count_dict = {'Ori':defaultdict(int),'Egg':defaultdict(int)}
  excludepattern = re.compile ("UNKNOWN_1_RHMK|TMK1_MDCK|AMNIOTIC_1_PRHMK_2|M1_RII1,C5|R1_C|R1_S|RII1_C|RII1_S|RIIX_C|RX_C|MDCK_1_RHMK|NC|_MK1")
  unpassagedpattern = re.compile("LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT")
  eggpattern = re.compile("AM[1-9]|E[1-7]|AMNIOTIC|EGG|EX|AM_[1-9]")
  cellpattern = re.compile("S[1-9]|SX|SIAT|MDCK|C[1-9]|CX|C_[1-9]|M[1-9]|MX|X[1-9]|^X_$")
  siatpattern = re.compile("^S[1-9]_$|SIAT2_SIAT1|SIAT3_SIAT1")
  monkeypattern=re.compile("TMK|RMK|RHMK|RII|PMK|R[1-9]|RX")
  siatexcludepattern=re.compile("SIAT|SX|S[1-9]")
  CountSeq = 0
  for record in records:
    CountSeq += 1
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
      Count_dict['Egg'][year]+=1
      Egg_dict[year].append(seq)
    elif unpassagedpattern.search(PSG): 
      Count_dict['Ori'][year]+=1
      Ori_dict[year].append(seq)
  return Egg_dict, Ori_dict, refseq, Count_dict

def writing_count_seq(Count_dict, countseqfile):
  outfile = open(countseqfile,'w')
  outfile.write("\t".join(['year','Ori','Egg'])+"\n")
  years = sorted(list(set(Count_dict['Ori'].keys()+Count_dict['Egg'].keys())))
  for year in years:
    outfile.write("\t".join(map(str,[year, Count_dict['Ori'][year], Count_dict['Egg'][year]]))+"\n")
  outfile.close()

def check_seq_length(seqs):
    seq_len = list(set([len(seq) for seq in seqs]))
    return seq_len

def dict_to_df(base_count_dict):
    dfs = []
    for position, position_bases in base_count_dict.iteritems():
        df = pd.DataFrame({'aa':position_bases.keys(),
                    'aa_count':position_bases.values(),
                    'position':position})
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def get_base_count(seqs):
    base_count_dict = defaultdict(lambda: defaultdict(int))
    for seq in seqs:
        for i, base in enumerate(seq):
            base_count_dict[i][base] += 1
    df = dict_to_df(base_count_dict)
    return df

def parse_year(Egg_dict, Ori_dict, year):
    egg_seqs = Egg_dict[year]
    ori_seqs = Ori_dict[year]

    egg_seq_len = check_seq_length(egg_seqs)
    assert len(egg_seq_len) == 1, 'Egg sequences not match'
    
    ori_seq_len = check_seq_length(ori_seqs)
    assert len(ori_seq_len) == 1, 'Ori sequences not match'

    assert ori_seq_len[0] == egg_seq_len[0], 'Ori sequences not matching egg seuqneces'
    egg_seq = get_base_count(egg_seqs) \
            .assign(data_type='egg')
    ori_seq = get_base_count(ori_seqs) \
            .assign(data_type='ori')

    df = pd.concat([egg_seq, ori_seq]) \
            .assign(year= year)
    return df

def count_to_freq(d):
    d['freq']  = d.aa_count/d.aa_count.sum()
    return d
    

def get_mutation_hotspot(tablename, cutoff=0, converter=None):
    df = pd.read_table(tablename)\
            .groupby(['data_type','year','position'],as_index=False) \
            .apply(count_to_freq) \
            .query("aa!='-'")    

    freq_df = df.pipe(pd.pivot_table,
                    index=['year','position','aa'],
                    columns=['data_type'],
                    values='freq', fill_value=0) \
            .rename(columns = {'egg':'egg_freq',
                               'ori':'ori_freq'}) \
            .reset_index() 

    count_df = df.pipe(pd.pivot_table,
                    index=['year','position','aa'],
                    columns=['data_type'],
                    values='aa_count', fill_value=0) \
            .reset_index() \
            .rename(columns = {'egg':'egg_count',
                               'ori':'ori_count'}) 
    df = count_df.merge(freq_df) \
            .assign(adjusted_position = lambda d: d.position.map(converter))
    df.to_csv('table/egg_vs_ori_bases.tsv',sep='\t',index=False)

    ori_df = df.query('ori_freq!=0')  \
        .pipe(lambda d: d[['position','aa','year']]) \
        .drop_duplicates() \
        .assign(label = 'ori')

    hotspot_positions_df = df\
            .merge(ori_df, how='left', on=['position','aa','year'])\
            .pipe(lambda d: d[pd.isnull(d.label)])  \
            .query('ori_freq==0 & egg_freq > %.3f' %(cutoff)) \
            .pipe(lambda d: d[['position','year','aa','egg_freq']])  \
            .rename(columns = {'egg_freq':'egg'})
    return hotspot_positions_df

def pair_wise_mutation(hotspot_positions_df, Egg_dict, RawPosToH3number_dict):
    double_mutants = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for year, year_df in hotspot_positions_df.groupby('year'):
        for i, (mut_pos1, mut_pos2) in enumerate(combinations(year_df.position.unique(),r=2)):
            key = ','.join(map(str,[mut_pos1, mut_pos2]))
            for seq in Egg_dict[str(year)]:
                value = seq[mut_pos1]+seq[mut_pos2]
                double_mutants[year][key][value] += 1
    
    dfs = []
    for year, year_dict in double_mutants.iteritems():
        for positions, position_dict in year_dict.iteritems():
            df = pd.DataFrame({'position':positions,
                'amino_combination':position_dict.keys(),
                'count':position_dict.values(),
                'year':year}) 
            dfs.append(df)
    
    df = pd.concat(dfs,axis=0) \
            .assign(position1 = lambda d: d.position.str.split(',',expand=True).iloc[:,0])\
            .assign(position2 = lambda d: d.position.str.split(',',expand=True).iloc[:,1])\
            .assign(adjusted_position1 = lambda d: d.position1.astype(int).map(RawPosToH3number_dict)) \
            .assign(adjusted_position2 = lambda d: d.position2.astype(int).map(RawPosToH3number_dict))

    df.to_csv('table/double_mutants.tsv',sep='\t', index=False)
    return 0

def parse_dictionaries(Egg_dict, Ori_dict, years, tablename, RawPosToH3number_dict):
    dfs = [parse_year(Egg_dict, Ori_dict, year) for year in list(years)]
    df = pd.concat(dfs) \
            .assign(adjusted_position = lambda d: d.position.map(RawPosToH3number_dict))
    df.to_csv(tablename, sep='\t', index=False)
    print 'Written %s' %tablename
    return 0

def write_dict(dictionary, dict_name):
    with open(dict_name, 'w') as pkl:
        pickle.dump(dictionary, pkl)
    return 0

def pos_adjust(refseq,posoffset):
  RawPosToH3number_dict = {}
  position = posoffset
  for n in range(len(refseq)):
    if refseq[n]=='-': 
      RawPosToH3number_dict[n] = -99
    else:
      position += 1
      RawPosToH3number_dict[n]=position
  return RawPosToH3number_dict

def make_single_mutation(df, position=1):
    pos_column = 'adjusted_position' + str(position)
    aa = 'aa' + str(position)

    single_df = df\
        .groupby(['year',pos_column, aa], as_index=False)\
        .agg({'count':'sum'})\
        .assign(aa_freq = lambda d: d\
                            .groupby([pos_column,'year'], as_index=False)['count']\
                            .transform(lambda x: x/x.sum())) \
        .rename(columns={'aa_freq':aa+'_freq'}) \
        .drop(['count'],axis=1)
    return single_df

def get_human():
    ori_df = pd.read_table('table/base_count_table.tsv') \
        .query("data_type!='egg'") \
        .pipe(lambda d: d[['aa','adjusted_position','year']]) \
        .drop_duplicates() \
        .assign(label = 'appeared_in_human')\
        .drop('year', axis=1)
    return ori_df

def rename_column(d, position=1):
    return d.rename(columns = {'aa': 'aa' + str(position),
                        'adjusted_position':'adjusted_position'+str(position)})

def filter_human(df):
    return df\
        .merge(rename_column(get_human(),position=1),
               how='left', on = ['adjusted_position1', 'aa1'])\
        .query("label!='appeared_in_human'")\
        .drop('label',axis=1) \
        .merge(rename_column(get_human(),position=2),
               how='left', on = ['adjusted_position2', 'aa2']) \
        .query('label!="appeared_in_human"')\
        .drop('label',axis=1)

def compile_to_mut_freq_table(doublemut_file):
  df = pd.read_table(doublemut_file) \
      .assign(aa1 = lambda d: d.amino_combination.str.slice(0,1))\
      .assign(aa2 = lambda d: d.amino_combination.str.slice(1,2)) \
      .drop_duplicates()
  df.head()
  fdf = df.merge(make_single_mutation(df,position=1),on=['adjusted_position1','year','aa1'],how='inner')\
    .merge(make_single_mutation(df, position=2),on=['adjusted_position2','year','aa2'],how='inner') \
    .assign(combine_freq = lambda d: d \
                .groupby(['year','adjusted_position1','adjusted_position2'], as_index=False)['count']\
                .transform(lambda x: x/x.sum()))  \
    .pipe(filter_human)\
    .assign(dependent_index=lambda d: d.combine_freq/(d.aa1_freq*d.aa2_freq))
  fdf.to_csv('table/mutation_freq.tsv',sep='\t',index=False)

def wrapper(alnfilename):
    posoffset = -16
    alnfile   = alnfilename
    tablename    = 'table/base_count_table.tsv'
    countseqfile = 'table/seq_count.tsv'
    records  = [record for record in SeqIO.parse(alnfile,"fasta")]
    Egg_dict, Ori_dict, refseq, Count_dict = ParseSeqs(records)
    writing_count_seq(Count_dict, countseqfile)
    write_dict(refseq, 'table/ref_seq.pkl')
    write_dict(Egg_dict, 'table/egg_dict.pkl')
    write_dict(Ori_dict, 'table/ori_dict.pkl')
    Egg_dict = pickle.load(open('table/egg_dict.pkl'))
    Ori_dict = pickle.load(open('table/ori_dict.pkl'))
    refseq = pickle.load(open('table/ref_seq.pkl'))
    RawPosToH3number_dict = pos_adjust(refseq,posoffset)
    write_dict(RawPosToH3number_dict, 'table/converter.pkl')
    years = set(Egg_dict.keys()).intersection(set(Ori_dict.keys()))
    parse_dictionaries(Egg_dict, Ori_dict, years, tablename, RawPosToH3number_dict)
    hotspot_positions_df = get_mutation_hotspot(tablename, cutoff=0, converter=RawPosToH3number_dict)

def main():
  wrapper('Fasta/HumanH3N2_All_2018.aln')

if __name__ == "__main__":
  main()
