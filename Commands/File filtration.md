# Filtration

## Python program 

The two following python programs are used :

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

### filtration gene

`awk -F "\t" -v OFS="\t" '$3=="gene" {print}' gencode.v37.annotation.gtf > gencode.annotation.filtered.by.gene.gtf`

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/creation.table.py --annotations gencode.annotation.filtered.by.gene.gtf --output table.gencode.v37.annotation.gene.bed`

- all type

level = 1

`awk -F "\t" -v OFS="\t" '$17=="1" {print}' table.gencode.v37.annotation.gene.bed > table.gencode.v37.annotation.gene.level1.bed`

level = 2

`awk -F "\t" -v OFS="\t" '$17=="2" {print}' table.gencode.v37.annotation.gene.bed > table.gencode.v37.annotation.gene.level2.bed`

`cat table.gencode.v37.annotation.gene.level1.bed table.gencode.v37.annotation.gene.level2.bed > table.gencode.v37.annotation.gene.level1.2.bed`

- protein coding

level = 1
gene_type = protein_coding

`awk -F "\t" -v OFS="\t" '$11=="protein_coding" {print}' table.gencode.v37.annotation.gene.bed | awk -F "\t" -v OFS="\t" '$17=="1" {print}' > table.gencode.v37.annotation.gene.level1.coding.bed`

level = 2
gene_type = protein_coding

`awk -F "\t" -v OFS="\t" '$11=="protein_coding" {print}' table.gencode.v37.annotation.gene.bed | awk -F "\t" -v OFS="\t" '$17=="2" {print}' > table.gencode.v37.annotation.gene.level2.coding.bed`

- lncRNA

level = 1
gene_type = lncRNA

`awk -F "\t" -v OFS="\t" '$11=="lncRNA" {print}' table.gencode.v37.annotation.gene.bed | awk -F "\t" -v OFS="\t" '$17=="1" {print}' > table.gencode.v37.annotation.gene.level1.lncRNA.bed`

level = 2
gene_type = lncRNA

`awk -F "\t" -v OFS="\t" '$11=="lncRNA" {print}' table.gencode.v37.annotation.gene.bed | awk -F "\t" -v OFS="\t" '$17=="2" {print}' > table.gencode.v37.annotation.gene.level2.lncRNA.bed`

### filtration transcript + exon

- protein coding

level = 1
TSL = 1,2,3
gene_type = protein_coding
transcript_type = protein_coding

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/filtration.annotation.file.py --annotations gencode.v37.annotation.gtf --output gencode.v37.annotation.transcript.exon.filtered.level1.coding.gtf --level 1 --type protein_coding --transcript yes --exon yes --TSL 1,2,3`

`awk -F "\t" -v OFS="\t" '$3=="exon" {print}' gencode.v37.annotation.transcript.exon.filtered.level1.coding.gtf | wc -l`

> transcript : 3450
> exon : 26556


level = 1,2
TSL = 1,2,3
gene_type = protein_coding
transcript_type = protein_coding

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/filtration.annotation.file.py --annotations gencode.v37.annotation.gtf --output gencode.v37.annotation.transcript.exon.filtered.level1.2.coding.gtf --level 1,2 --type protein_coding --transcript yes --exon yes --TSL 1,2,3`

`awk -F "\t" -v OFS="\t" '$3=="transcript" {print}' gencode.v37.annotation.transcript.exon.filtered.level1.2.coding.gtf | wc -l`

> transcript : 17512

- lncRNA

level = 1
TSL = 1,2,3
gene_type = lncRNA
transcript_type = lncRNA

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/filtration.annotation.file.py --annotations gencode.v37.annotation.gtf --output gencode.v37.annotation.transcript.exon.filtered.level1.lncRNA.gtf --level 1 --type lncRNA --transcript yes --exon yes --TSL 1,2,3`

level = 1,2
TSL = 1,2,3
gene_type = lncRNA
transcript_type = lncRNA

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/filtration.annotation.file.py --annotations gencode.v37.annotation.gtf --output gencode.v37.annotation.transcript.exon.filtered.level1.2.lncRNA.gtf --level 1,2 --type lncRNA --transcript yes --exon yes --TSL 1,2,3`

# Creation table

Creation of tables that will be imported into R

## gene

- all type 

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.gencode.v37.annotation.gene.level1.bed > table.gene.annotation.level1.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.gencode.v37.annotation.gene.level1.2.bed > table.gene.annotation.level1.2.bed`

- protein_coding
- lncRNA

`sort -k1,1 -k2,2n table.gene.annotation.level1.bed > sorted.table.gene.annotation.level1.bed`
`sort -k1,1 -k2,2n table.gene.annotation.level1.2.bed > sorted.table.gene.annotation.level1.2.bed`

- count 

`bedtools map -a sorted.table.gene.annotation.level1.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.gene.annotation.occur.level1.bed`

`bedtools map -a sorted.table.gene.annotation.level1.2.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.gene.annotation.occur.level1.2.bed`


## intergenic 

file containing the sizes of each chromosome : 
`awk -F "\t" -v OFS="\t" '{print($1"\t1\t"$2)}' hg38.chrom.sizes > hg38.size.bed`

bedtools subtract:

`bedtools subtract -a hg38.size.bed -b table.gene.annotation.level1.bed > intergenicN.level1.bed`

sorted the file : 

`sort -k1,1 -k2,2n intergenicN.level1.bed > sorted.intergenicN.level1.bed`

This file contains the intergenic sequences and the N

- Recover the position of N in the genome fasta file with a python program

position.chr.py : returns a bed file containing the positions of all N in the file -> positionN.bed 

This file is sorted : 

`sort -k1,1 -k2,2n positionN.bed > sorted.positionN.bed`

- Subtract 

`bedtools subtract -a sorted.intergenicN.level1.bed -b sorted.positionN.bed > intergenic.level1.bed`

- count 

`bedtools map -a intergenic.level1.2.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.intergenic.annotation.level1.2.bed`

- put in the form of a table

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"$4)}' table.intergenic.annotation.level1.bed > table.intergenic.level1.bed`

## transcript + exon 

- protein coding

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/creation.table.py --annotations Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.coding.gtf --output table.transcript.exon.level1.coding.bed`

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/creation.table.py --annotations Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.coding.gtf --output table.transcript.exon.level1.2.coding.bed`

- lncRNA

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/creation.table.py --annotations Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.lncRNA.gtf --output table.transcript.exon.level1.lncRNA.bed`

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/creation.table.py --annotations Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.lncRNA.gtf --output table.transcript.exon.level1.2.lncRNA.bed`

- put in the form of a table

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.transcript.exon.level1.coding.bed > transcript.exon.level1.coding.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.transcript.exon.level1.2.coding.bed > transcript.exon.level1.2.coding.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.transcript.exon.level1.lncRNA.bed > transcript.exon.level1.lncRNA.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$4"\t"$5"\t"$14"\t.\t"$7"\t"$3"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19)}' table.transcript.exon.level1.2.lncRNA.bed > transcript.exon.level1.2.lncRNA.bed`

`sort -k1,1 -k2,2n transcript.exon.level1.coding.bed > sorted.transcript.exon.level1.coding.bed`
`sort -k1,1 -k2,2n transcript.exon.level1.2.coding.bed > sorted.transcript.exon.level1.2.coding.bed`
`sort -k1,1 -k2,2n transcript.exon.level1.lncRNA.bed > sorted.transcript.exon.level1.lncRNA.bed`
`sort -k1,1 -k2,2n transcript.exon.level1.2.lncRNA.bed > sorted.transcript.exon.level1.2.lncRNA.bed`

- count 

`bedtools map -a sorted.transcript.exon.level1.coding.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.transcript.exon.level1.coding.occur.bed`

`bedtools map -a sorted.transcript.exon.level1.2.coding.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.transcript.exon.level1.2.coding.occur.bed`

`bedtools map -a sorted.transcript.exon.level1.lncRNA.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.transcript.exon.level1.lncRNA.occur.bed`

`bedtools map -a sorted.transcript.exon.level1.2.lncRNA.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.transcript.exon.level1.2.lncRNA.occur.bed`

## intron

To do with level 1, level1.2 in coding and lncRNA

- Recover positions of transcripts

`awk -F "\t" -v OFS="\t" '$7=="transcript" {print($0)}' sorted.transcript.exon.level1.coding.bed > transcript.level1.coding.bed`
`awk -F "\t" -v OFS="\t" '$7=="transcript" {print($0)}' sorted.transcript.exon.level1.2.coding.bed > transcript.level1.2.coding.bed`
`awk -F "\t" -v OFS="\t" '$7=="transcript" {print($0)}' sorted.transcript.exon.level1.lncRNA.bed > transcript.level1.lncRNA.bed`
`awk -F "\t" -v OFS="\t" '$7=="transcript" {print($0)}' sorted.transcript.exon.level1.2.lncRNA.bed > transcript.level1.2.lncRNA.bed`

- Recover positions of exons 

`awk -F "\t" -v OFS="\t" '$7=="exon" {print($0)}' sorted.transcript.exon.level1.coding.bed > exon.level1.coding.bed`
`awk -F "\t" -v OFS="\t" '$7=="exon" {print($0)}' sorted.transcript.exon.level1.2.coding.bed > exon.level1.2.coding.bed`
`awk -F "\t" -v OFS="\t" '$7=="exon" {print($0)}' sorted.transcript.exon.level1.lncRNA.bed > exon.level1.lncRNA.bed`
`awk -F "\t" -v OFS="\t" '$7=="exon" {print($0)}' sorted.transcript.exon.level1.2.lncRNA.bed > exon.level1.2.lncRNA.bed`

- subtract position of exons to positions of transcripts

`bedtools  subtract -a transcript.level1.coding.bed -b exon.level1.coding.bed > intron.level1.coding.bed`
`bedtools  subtract -a transcript.level1.2.coding.bed -b exon.level1.2.coding.bed > intron.level1.2.coding.bed`
`bedtools  subtract -a transcript.level1.lncRNA.bed -b exon.level1.lncRNA.bed > intron.level1.lncRNA.bed`
`bedtools  subtract -a transcript.level1.2.lncRNA.bed -b exon.level1.2.lncRNA.bed > intron.level1.2.lncRNA.bed`

- put in the form of a table

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tintron\t"$3-$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18)}' intron.level1.coding.bed > table.intron.level1.coding.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tintron\t"$3-$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18)}' intron.level1.2.coding.bed > table.intron.level1.2.coding.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tintron\t"$3-$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18)}' intron.level1.lncRNA.bed > table.intron.level1.lncRNA.bed`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tintron\t"$3-$2"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18)}' intron.level1.2.lncRNA.bed > table.intron.level1.2.lncRNA.bed`

`sort -k1,1 -k2,2n table.intron.level1.coding.bed > sorted.table.intron.level1.coding.bed`
`sort -k1,1 -k2,2n table.intron.level1.2.coding.bed > sorted.table.intron.level1.2.coding.bed`
`sort -k1,1 -k2,2n table.intron.level1.lncRNA.bed > sorted.table.intron.level1.lncRNA.bed`
`sort -k1,1 -k2,2n table.intron.level1.2.lncRNA.bed > sorted.table.intron.level1.2.lncRNA.bed`

- count

`bedtools map -a sorted.table.intron.level1.coding.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.intron.level1.coding.occur.bed`

`bedtools map -a sorted.table.intron.level1.2.coding.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.intron.level1.2.coding.occur.bed`

`bedtools map -a sorted.table.intron.level1.lncRNA.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.intron.level1.lncRNA.occur.bed`

`bedtools map -a sorted.table.intron.level1.2.lncRNA.bed -b sorted.lexicographically.positions.motifL1.filtered.0based.bed -c 5 -o count -sorted > table.intron.level1.2.lncRNA.occur.bed`


# Deeptools

Creation of region files that will be used with deeptools 

## transcript

- protein coding

`awk -F "\t" -v OFS="\t" '$3=="transcript" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.coding.gtf > transcript.coding.level1.bed`

`awk -F "\t" -v OFS="\t" '$3=="transcript" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.coding.gtf > transcript.coding.level1.2.bed`

- lncRNA

`awk -F "\t" -v OFS="\t" '$3=="transcript" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.lncRNA.gtf > transcript.lncRNA.level1.bed`

`awk -F "\t" -v OFS="\t" '$3=="transcript" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.lncRNA.gtf > transcript.lncRNA.level1.2.bed`

## exon

- protein coding

`awk -F "\t" -v OFS="\t" '$3=="exon" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.coding.gtf > exon.coding.level1.bed`

`awk -F "\t" -v OFS="\t" '$3=="exon" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.coding.gtf > exon.coding.level1.2.bed`

- lncRNA

`awk -F "\t" -v OFS="\t" '$3=="exon" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.1/gencode.v37.annotation.transcript.exon.filtered.level1.lncRNA.gtf > exon.lncRNA.level1.bed`

`awk -F "\t" -v OFS="\t" '$3=="exon" {print($1"\t"$4"\t"$5"\t"$7)}' Stage_laurine/Desktop/STAGE/Workspace/deeptools/annotation/level.12/gencode.v37.annotation.transcript.exon.filtered.level1.2.lncRNA.gtf > exon.lncRNA.level1.2.bed`

