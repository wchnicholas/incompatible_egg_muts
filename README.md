The sequence analysis is divided into two parts. [Part I](https://github.com/wchnicholas/incompatible_egg_muts#part-i-co-occurrence-analysis-of-egg-mutations) describes the analysis for co-occuring egg-adaptive mutations in sequence database. [Part II](https://github.com/wchnicholas/incompatible_egg_muts#part-ii-deep-sequencing-analysis-of-egg-mutations) describes the analysis for deep sequencing of egg-passaging experiments.
## PART I: CO-OCCURRENCE ANALYSIS OF EGG MUTATIONS
### CO-OCCURRENCE ANALYSIS: INPUT FILE
* [./Fasta/HumanH3N2\_All\_2018.aln](./Fasta/HumanH3N2_All_2018.aln): Human H3N2 HA sequences downloaded from [GISAID](https://www.gisaid.org/)

### CO-OCCURRENCE ANALYSIS: ANALYSIS PIPELINE
1. [./script/Parse\_GISAID\_aln.py](./script/Parse_GISAID_aln.py): Identify egg-adaptive mutations, written by [Douglas Wu](https://wckdouglas.github.io/)
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
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA532726](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA532726), should be placed in fastq/ folder. The filename for read 1 should match those described in [./data/Sample\_name.tsv](./data/Sample_name.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".

### DEEP SEQUENCING: ANALYSIS PIPELINE
1. [./script/Bris07\_Egg\_read\_to\_count.py](./script/Bris07_Egg_read_to_count.py): Convert raw reads to mutation count
    - Input file: Raw sequencing reads in fastq/ folder
    - Output files: result/\*\_count.tsv
2. [./script/Bris07\_Egg\_count\_to\_freq.py](./script/Bris07_Egg_count_to_freq.py): Convert mutation count to mutation frequency
    - Input file: result/\*\_count.tsv
    - Output files: [./result/mut\_freq.tsv](result/mut_freq.tsv)

### DEEP SEQUENCING: PLOTTING
1. [./script/plot\_deep\_seq\_mut\_freq.R](script/plot_deep_seq_mut_freq.R): Plot the mutation frequency in different passages
    - Inupt file: [./result/mut\_freq.tsv](result/mut_freq.tsv)
    - Output file: graph/Freq\_\*.png
