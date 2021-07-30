# Research L1 motif

## Checking on chr22

1. Genome-wide motif scanning 

`scanMotifGenomeWide.pl /Users/laurinerolland/Desktop/STAGE/Workspace/Motif/hg19.l1neo.soni.loc.helas3.score.homer.motif /Users/laurinerolland/Desktop/STAGE/Workspace/Fasta/2021-04-29/chr22.fa -bed > site_chr22.bed`

2. Scores under the threshold 4.96256 are filtered out

`awk -F "\t" -v OFS="\t" '$5>=4.96256 {print}' site_chr22.bed > chr22.score.filtered.bed`

3. Transition from 1-based to 0-based 

`awk -F"\t" -v OFS="\t" '{print($1,$2-1,$3,$4,$5,$6)}' chr22.score.filtered.1-based.bed > motifL1.score.filtered.0based.bed`

4. Nucleotide sequence extraction

`bedtools getfasta  -s -fi ~/Desktop/STAGE/Workspace/Fasta/2021-04-29/chr22.fa -bed ~/Desktop/STAGE/Workspace/Bed/2021-04-22/motifL1.score.filtered.0based.bed -tab > seq.site.chr22.fa`

5. Display modification 

`awk -F"\t" -v OFS="\t" '{split($1,a,"("); split(a[1],b,"-"); split(a[2],c,")"); split(b[1],d,":"); print(d[1],d[2],b[2],c[1],toupper($2))}' seq.site.chr22.fa > seq.chr22.maj.bed`

6. The nucleotide column of the fasta file is added to the bed file

`bedtools map  -f 1 -r -a motifL1.score.filtered.0based.bed -b seq.chr22.maj.bed -c 5 -o distinct > chr22.site.seq.motifL1.bed`

7. Generation of the weblogo

`cut -f 7 chr22.site.seq.motifL1.bed > motifL1.seq.chr22.fa`

`weblogo -c classic -s large -i -3 < motifL1.seq.chr22.fa > sequence.WebLogo.motifL1.eps`


## Searching for L1 motifs in the genome

`scanMotifGenomeWide.pl hg19.l1neo.soni.loc.helas3.score.homer.motif hg38.fa -bed -p 20 > postions.motifL1.1-based.bed`

Compte professorx:

`scanMotifGenomeWide.pl hg19.l1neo.soni.loc.helas3.score.homer.motif hg38 -bed > postions.motifL1.1-based.bed`

`scp lrolland@134.59.51.200:\postions.motifL1.1-based.bed .`

- Scores under the threshold 4.96256 are filtered out

`awk '$4>=4.96256 {print}' postions.motifL1.1-based.bed >  postions.motifL1.filtered.1.based.bed`

- Transition from 0-based to 1-based

Filtered 

`awk -F"\t" -v OFS="\t" '{print($1,$2-1,$3,$4,$5,$6)}' postions.motifL1.filtered.1.based.bed > positions.motifL1.filtered.0based.bed`

Not filtered 

`awk -F"\t" -v OFS="\t" '{print($1,$2-1,$3,$4,$5,$6)}' postions.motifL1.1-based.bed > positions.motifL1.0based.bed`
