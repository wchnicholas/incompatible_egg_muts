The sequence analysis is divided into two parts. [Part I](https://github.com/wchnicholas/incompatible_egg_muts#part-i-co-occurrence-analysis-of-egg-mutations) describes the analysis for co-occuring egg-adaptive mutations in sequence database. [Part II](https://github.com/wchnicholas/incompatible_egg_muts#part-ii-deep-sequencing-analysis-of-egg-mutations) describes the analysis for deep sequencing of egg-passaging experiments.
## PART I: CO-OCCURRENCE ANALYSIS OF EGG MUTATIONS
### CO-OCCURRENCE ANALYSIS: INPUT FILE
[./Fasta/HumanH3N2\_All\_2018.aln](./Fasta/HumanH3N2_All_2018.aln): Human H3N2 HA sequences downloaded from [GISAID](https://www.gisaid.org/)

### CO-OCCURRENCE ANALYSIS: ANALYSIS PIPELINE
1. [./script/Parse\_GISAID\_aln.py](./script/Parse_GISAID_aln.py): Identify egg-adaptive mutations and analyze their co-occurrence frequencies
    - Input file: [./Fasta/HumanH3N2\_All\_2018.aln](./Fasta/HumanH3N2_All_2018.aln)
    - Output file:
      - [./table/base\_count\_table.tsv](./table/base_count_table.tsv)
      - [./table/converter.pkl](./table/converter.pkl)
      - [./table/egg\_dict.pkl](./table/egg_dict.pkl)
      - [./table/egg\_vs\_ori\_bases.tsv](./table/egg_vs_ori_bases.tsv)
      - [./table/ori\_dict.pkl](./table/ori_dict.pkl)
      - [./table/ref\_seq.pkl](./table/ref_seq.pkl)
      - [./table/seq\_count.tsv](./table/seq_count.tsv)
2. [./script/filter\_mut\_table.py](./script/filter_mut_table.py): Identify major egg-adaptive mutations
    - Input file: [./table/egg\_vs\_ori\_bases.tsv](./table/egg_vs_ori_bases.tsv)
    - Output file: 
      - [./table/mut\_freq\_egg\_filtered.tsv](./table/mut_freq_egg_filtered.tsv)
      - [./table/mut\_count\_year\_egg.tsv](./table/mut_count_year_egg.tsv)
3. [./script/egg\_mut\_pair\_freq.py](./script/egg_mut_pair_freq.py): Compute the co-occurrence frequencies of egg-adaptive mutations
    - Input file: 
      - [./table/mut\_count\_year\_egg.tsv](./table/mut_count_year_egg.tsv)
      - [./Fasta/HumanH3N2\_All\_2018.aln](./Fasta/HumanH3N2_All_2018.aln)
    - Output file: 
      - [./table/mut\_double\_freq.tsv](./table/mut_double_freq.tsv)
      - [./table/mut\_single\_freq.tsv](./table/mut_single_freq.tsv)

### CO-OCCURRENCE ANALYSIS: PLOTTING
1. [./script/plot\_mut\_freq\_egg.R](./script/plot_mut_freq_egg.R): Plot the occurrence frequencies of egg-adaptive mutations observed in different year
    - Input file: 
      - [./table/mut\_freq\_egg\_filtered.tsv](./table/mut_freq_egg_filtered.tsv)
      - [./table/seq\_count.tsv](./table/seq_count.tsv)
    - Output file: [./graph/MutFreq\_Egg.png](./graph/MutFreq_Egg.png)
2. [./script/egg\_mut\_pair\_network.py](./script/egg_mut_pair_network.py): Plot the network diagram that describes the co-occurrence frequency of egg-adaptive mutations
    - Input file: 
      - [./table/mut\_double\_freq.tsv](./table/mut_double_freq.tsv)
      - [./table/mut\_single\_freq.tsv](./table/mut_single_freq.tsv)
    - Output file: 
      - [./graph/EggMutPair\_network.dot](graph/EggMutPair_network.dot)
      - [./graph/EggMutPair\_network.png](graph/EggMutPair_network.png)

## PART II: DEEP SEQUENCING ANALYSIS OF EGG MUTATIONS
### DEEP SEQUENCING: INPUT FILE


### DEEP SEQUENCING: ANALYSIS PIPELINE
