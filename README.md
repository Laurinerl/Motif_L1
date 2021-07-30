# Motif_L1
Work done during the internship from  01/02/2021 to 30/07/2021

## Installation

### Dependencied 

- bedtools (v2.30.0)
- homer (v4.11)
- Rstudio (Version 1.4.1103)
- deepTools (v3.5.0)

### Other requirements

- A human reference genome sequence : hg38.fa
- Annotation file : gencode.v37.annotation.gtf

## Usage 

The commands lines made are detailed in the folder "Commands". 

### Python program 

Four python programs have been created:

- creation_PWM.py</li

Allows to generate frequency matrices. It takes as input :
	- a bed file which contains the positions of the motifs 
	- a fasta file which contains the nucleotide sequence of the genome
By default the output matrix will be in homer format but it is possible to use other formats like Jaspar, transfac or raw.

- score_PWM.py

Allows to calculate the score of the motifs of a bed file from a motif matrix. It takes as input :
	- the bed file which presents the positions of the patterns 
	- a fasta file which contains the nucleotide sequence of the genome
	- the matrix 

- filtration.annotation.file.py 

Filtration of annotation files according to the annotation level

>level 1 (verified loci)
>level 2 (manually annotated loci)
>level 3 (automatically annotated loci)

and according to the gene type and transcript type (protein_coding, lncRNA...)

- creation.table.py 

Allows to put annotation files in the form of a table, the informations are separated in columns. This program also adds a column with the length of the annotated region

- position.chr.py

Returns a bed file containing the positions of all N in the file
