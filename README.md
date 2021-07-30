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

They are contained in the file hg38_files_used

## Usage 

The commands lines made are detailed in the folder "Commands". 

### 1. Generate a file presenting the positions and scores of L1 motifs across the human genome 

The file "File filtration.md" details how to find the positions of the motifs that match the L1 motif frequency matrix on the human genome. 

Python programs creation_PWM.py and score_PWM.py are used here to find the frequency matrix of the L1 motif and make some controls

### 2. Filtration of the annotation file

The file "File filtration.md" details how the annotation file was filtered.

	- Intergenic and intronic regions were identified
	- Tables that will be imported into R are created 
	- Annotation files that will be used with deepTools are created 

Python programs filtration.annotation.file.py, creation.table.py and position.chr.py are used here.

### 3. Generation of graphs comparing L1 motif density between regions (gene/intergenic, intron/exon)

The file "R_table.md" details the command lines used to generate these comparative graphs on R

### 4. Generation of graphs comparing the density of L1 patterns at borders in different regions (gene/intergenic, intron/exon)

The file "deeptools.md" details the command lines used to generate these comparative graphs with deeptools

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

Filtration of annotation files according to the annotation level and according to the gene type and transcript type (protein_coding, lncRNA...)

>level 1 (verified loci)
>level 2 (manually annotated loci)
>level 3 (automatically annotated loci)

- creation.table.py 

Allows to put annotation files in the form of a table, the informations are separated in columns. This program also adds a column with the length of the annotated region

- position.chr.py

Returns a bed file containing the positions of all N in the file
