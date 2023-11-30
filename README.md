# DH307: Crafting an RNA-Seq Pipeline for Analysis of Gene Expression Profiles of Parkinsonian Samples
Installation of tools, usage and troubleshooting documentation: [Google Drive](https://docs.google.com/document/d/1p-3MiJXUPp2tVudpx6Qg0e5Mcc1piEOz_jWTYXM1k_c/edit?usp=sharing)

## Usage
### For paired ended reads
To generate gtf alignment file for SRR19538455 run<br>
```python3 pe_rna_seq.py --runall 1 --SRR SRR19538455```<br>
To generate gtf alignment files for all runs from SRR19538455 to SRR19538490<br>
```python3 pe_rna_seq.py --runall 1 --SRR SRR19538455 SRR19538490 1 --i 1```<br>
To generate gtf alignment files for all alternate runs from SRR19538455 to SRR19538490<br>
```python3 pe_rna_seq.py --runall 1 --SRR SRR19548455 SRR19548457 2 --i 1```<br>
To generate gtf alignment files for all runs from a custom list of variable length, say [SRR1, SRR2, SRR3]<br>
```python3 pe_rna_seq.py --runall 1 --SRR SRR1 SRR2 SRR3```<br>
To run any specific step of the pipeline for SRR19538455, say stringtie<br>
```python3 pe_rna_seq.py --stringtie 1 --SRR SRR19538455```<br>

### For single ended reads
To generate gtf alignment files for all runs from a custom list of variable length, say [SRR1, SRR2, SRR3]<br>
```python3 se_rna_seq.py --runall 1 --SRR SRR1 SRR2 SRR3```<br>
To run any specific step of the pipeline for SRR1, say stringtie<br>
```python3 pe_rna_seq.py --stringtie 1 --SRR SRR1```<br>

### Generating genes and transcript counts
After generating all gtf alignment files for GSE205450 in the stringtie folder, first create lst_GSE205450.txt file by running the following command<br>
```python3 create_lst.py GSE205450```<br>
We then run the prepDE.py3 file to generate csv counts.<br>
```python3 prepDE.py3 -i lst_GSE205450.txt```<br>

### DESeq
First generate metadata for the SRR files(of lets say GSE169755 dataset) by first downloading the GSE169755_series_matrix.txt and then run metadata.py<br>
```python3 metadata.py GSE169755_series_matrix.txt```<br>
The above saves metadata in metadata.csv, now open an R workspace and run the following files<br>
```setwd("/path/to/your/working/directory")```<br>
```source("DESeq2.R")```<br>
```source("plot_generation.R")```<br>
The above generates plots and save them in pdf files.<br>
