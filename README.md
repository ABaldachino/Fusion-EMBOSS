# Fusion-EMBOSS

The fusion script is able to obtain the GC percentages and the codon usage, and to search the restriction site in the DNA sequence if a file contening the restriction enzyme is given.

## Code exemple

Obtention of an output file showing the GC percentages and the codon usage in the ORFs found in the sequences from the file(s).
There can be multiple files containing sequence(s). Accepted formats are: (Multi)FASTA (.fasta/.fa) and EMBL (.embl).
If working on operon structures, you may consider to seperate the different ORFs from one another as the program can isolate only one ORF per sequence.
If working on spliced genes, you may consider to give the mature-RNA sequences as the ORF isolation relies on the reading frame.
```{bash}
python3 fusion.py sequence.fa (sequence2.fa, ...)
```
Obtention of the restriction site, option of  
```{bash}
python3 fusion.py sequence.fa -e enz.txt
```
## Motivation
For our Master degree in Bioinformatics, we had to merge two individually written scripts, and we decided to do it on git for training. 

## Installation
```{bash}
git clone git@github.com:ABaldachino/Fusion-EMBOSS.git
```
