# Deeptools

## Creation of bigWig file

- creation of the file

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$5)}' positions.motifL1.0based.bed > positions.motifL1.bedgraph`

`sort -k1,1 -k2,2n positions.motifL1.bedgraph > sorted.positions.motifL1.bedgraph`

`bedGraphToBigWig sorted.positions.motifL1.bedgraph hg38.chrom.sizes bigWig.positions.motifL1.bw`

`awk -F "\t" -v OFS="\t" '{print($1"\t"$2"\t"$3"\t"$4)}' transcript.filtered.bed > transcript.filteredbis.bed`

`sort -k1,1 -k2,2n transcript.filteredbis.bed > sorted.transcriptfiltered.bedgraph`

## Creation of bigWig AT%

- dowload "gc5Base.bw"

- bigwig to bedgraph

`bigWigToBedGraph gc5Base.bw gc5Base.bedgraph`

`awk -F "\t" -v OFS="\t" '{print("\""$1"\",")}' hg38.chrom.sizes | tr -d "\n" > chr.hg38.list.txt`

- GC% to AT%

Program python 

`python Stage_laurine/Desktop/STAGE/Workspace/Python/Program/GCtoAT.py --bedgraph gc5Base.bedgraph --output at5Base.bedgraph`


`bedGraphToBigWig at5Base.bedgraph hg38.chrom.sizes bigWig.at5Base.bw`

# Transcript protein coding 

## matrix

### scale-region

- level 1/2 

computeMatrix scale-regions --regionBodyLength 2000 -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.trancript.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions transcript.coding.level1.2.bin40.bed 

- AT%

computeMatrix scale-regions --regionBodyLength 2000 -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.trancript.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions AT.transcript.coding.level1.2.bin40.bed 

- both 

computeMatrix scale-regions --regionBodyLength 2000 -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.TA.L1.trancript.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions TA.L1.transcript.coding.level1.2.bin40.bed 


### reference-point

- level 1

computeMatrix reference-point --referencePoint TSS -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.transcript.TSS.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions TSS.transcript.coding.level1.2.bin40.bed 

computeMatrix reference-point --referencePoint TES -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.transcript.TES.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions TES.transcript.coding.level1.2.bin40.bed 

- AT%

computeMatrix reference-point --referencePoint TSS -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.transcript.TSS.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions TSS.AT.transcript.coding.level1.2.bin40.bed 

computeMatrix reference-point --referencePoint TES -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.transcript.TES.coding.level1.2.bin40.gz -p 6 --outFileSortedRegions TES.AT.transcript.coding.level1.2.bin40.bed 

## heatmap

- scale-region

plotHeatmap -m  matrix.trancript.coding.level1.2.bin40.gz -out heatmap.transcript.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"


- reference-point

plotHeatmap -m matrix.transcript.TSS.coding.level1.2.bin40.gz -out heatmap.transcript.TSS.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

plotHeatmap -m matrix.transcript.TES.coding.level1.2.bin40.gz -out heatmap.transcript.TES.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"


## profile 

### level 1

- scale-region

plotProfile -m matrix.trancript.coding.level1.2.bin40.gz -out profile.transcript.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.trancript.coding.level1.2.bin40.gz -out profile.AT.transcript.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"

- both 

plotProfile -m matrix.TA.L1.trancript.coding.level1.2.bin40.gz -out profile.AT.transcript.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "Transcript type : protein coding" --colors blue red --samplesLabel "L1 motif density" "AT%" --ymax 1

- reference-point

plotProfile -m matrix.transcript.TSS.coding.level1.2.bin40.gz -out profile.transcript.TSS.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

plotProfile -m matrix.transcript.TES.coding.level1.2.bin40.gz -out profile.transcript.TES.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.transcript.TSS.coding.level1.2.bin40.gz -out profile.AT.transcript.TSS.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"

plotProfile -m matrix.AT.transcript.TES.coding.level1.2.bin40.gz -out profile.AT.transcript.TES.coding.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"

# Transcript lncRNA 

## matrix

### scale-region

- level 1

computeMatrix scale-regions --regionBodyLength 2000 -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.trancript.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions transcript.lncRNA.level1.2.bin40.bed 

- AT%

computeMatrix scale-regions --regionBodyLength 2000 -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.trancript.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions AT.transcript.lncRNA.level1.2.bin40.bed 

### reference-point

- level 1

computeMatrix reference-point --referencePoint TSS -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.transcript.TSS.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions TSS.transcript.lncRNA.level1.2.bin40.bed 

computeMatrix reference-point --referencePoint TES -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.transcript.TES.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions TES.transcript.lncRNA.level1.2.bin40.bed 

- AT%

computeMatrix reference-point --referencePoint TSS -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.transcript.TSS.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions AT.TSS.transcript.lncRNA.level1.2.bin40.bed 

computeMatrix reference-point --referencePoint TES -b 5000 -a 5000 --missingDataAsZero --binSize 40 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/transcript.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.transcript.TES.lncRNA.level1.2.bin40.gz -p 6 --outFileSortedRegions AT.TES.transcript.lncRNA.level1.2.bin40.bed 

## heatmap

- scale-region

plotHeatmap -m matrix.trancript.lncRNA.level1.2.bin40.gz -out heatmap.transcript.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- reference-point

plotHeatmap -m matrix.transcript.TSS.lncRNA.level1.2.bin40.gz -out heatmap.transcript.TSS.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

plotHeatmap -m matrix.transcript.TES.lncRNA.level1.2.bin40.gz -out heatmap.transcript.TES.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

## profile 

### level 1

- scale-region

plotProfile -m matrix.trancript.lncRNA.level1.2.bin40.gz -out profile.transcript.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.trancript.lncRNA.level1.2.bin40.gz -out profile.AT.transcript.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

- reference-point

plotProfile -m matrix.transcript.TSS.lncRNA.level1.2.bin40.gz -out profile.transcript.TSS.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

plotProfile -m matrix.transcript.TES.lncRNA.level1.2.bin40.gz -out profile.transcript.TES.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.transcript.TSS.lncRNA.level1.2.bin40.gz -out profile.AT.transcript.TSS.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

plotProfile -m matrix.AT.transcript.TES.lncRNA.level1.2.bin40.gz -out profile.transcript.TES.lncRNA.level1.2.bin40.png --plotTitle "Level 1/2" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

# Exon protein coding

## matrix

### scale-region

- level 1

computeMatrix scale-regions --regionBodyLength 100 -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions exon.coding.level1.bin20.bed 

- AT%

computeMatrix scale-regions --regionBodyLength 100 -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.exon.coding.level1.bin20.bed 

### reference-point

- level 1

computeMatrix reference-point --referencePoint TSS -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.start.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions start.exon.coding.level1.2.bin20.bed 

computeMatrix reference-point --referencePoint TES -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.stop.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions stop.exon.coding.level1.bin20.bed 

- AT%

computeMatrix reference-point --referencePoint TSS -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.start.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.start.exon.coding.level1.2.bin20.bed 

computeMatrix reference-point --referencePoint TES -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.coding.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.stop.coding.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.stop.exon.coding.level1.bin20.bed 

## heatmap

- scale-region

plotHeatmap -m matrix.exon.coding.level1.2.bin20.gz -out heatmap.exon.coding.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

- reference-point

plotHeatmap -m matrix.exon.start.coding.level1.2.bin20.gz -out heatmap.exon.start.coding.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

plotHeatmap -m matrix.exon.stop.coding.level1.2.bin20.gz -out heatmap.exon.stop.coding.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

## profile 

### level 1

- scale-region

plotProfile -m matrix.exon.coding.level1.2.bin20.gz -out profile.exon.coding.level1.bin20.png --plotTitle "Exons" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.exon.coding.level1.2.bin20.gz -out profile.AT.exon.coding.level1.bin20.png --plotTitle "Exons" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"

- reference-point

plotProfile -m matrix.exon.start.coding.level1.2.bin20.gz -out profile.exon.start.coding.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

plotProfile -m matrix.exon.stop.coding.level1.2.bin20.gz -out profile.exon.stop.coding.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.exon.start.coding.level1.2.bin20.gz -out profile.AT.exon.start.coding.level1.2.bin20.png --plotTitle "Exons" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"

plotProfile -m matrix.AT.exon.stop.coding.level1.2.bin20.gz -out profile.AT.exon.stop.coding.level1.2.bin20.png --plotTitle "Exons" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : protein coding" --regionsLabel "AT%"


# Exons lncRNA 

## matrix

### scale-region

- level 1

computeMatrix scale-regions --regionBodyLength 100 -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions exon.lncRNA.level1.2.bin20.bed 

- AT%

computeMatrix scale-regions --regionBodyLength 100 -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.exon.lncRNA.level1.2.bin20.bed 

### reference-point

- level 1

computeMatrix reference-point --referencePoint TSS -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.start.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions start.exon.lncRNA.level1.2.bin20.bed 

computeMatrix reference-point --referencePoint TES -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.positions.motifL1.bw -o matrix.exon.stop.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions stop.exon.lncRNA.level1.2.bin20.bed 

- AT%

computeMatrix reference-point --referencePoint TSS -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.start.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.start.exon.lncRNA.level1.2.bin20.bed 

computeMatrix reference-point --referencePoint TES -b 200 -a 200 --missingDataAsZero --binSize 20 -R Desktop/STAGE/Workspace/deeptools/annotation/level.12/exon.lncRNA.level1.2.bed -S Desktop/STAGE/Workspace/deeptools/bigwig/bigWig.at5Base.bw -o matrix.AT.exon.stop.lncRNA.level1.2.bin20.gz -p 6 --outFileSortedRegions AT.stop.exon.lncRNA.level1.2.bin20.bed 

## heatmap

- scale-region

plotHeatmap -m matrix.exon.lncRNA.level1.2.bin20.gz -out heatmap.exon.lncRNA.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- reference-point

plotHeatmap -m matrix.exon.start.lncRNA.level1.2.bin20.gz -out heatmap.exon.start.lncRNA.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

plotHeatmap -m matrix.exon.stop.lncRNA.level1.2.bin20.gz -out heatmap.exon.stop.lncRNA.level1.2.bin40.png --plotTitle "Exons" -T "L1 motif density" --colorMap Blues --whatToShow 'heatmap and colorbar' --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

## profile 

### level 1

- scale-region

plotProfile -m matrix.exon.lncRNA.level1.2.bin20.gz -out profile.exon.lncRNA.level1.2.bin20.png --plotTitle "Exons" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.exon.lncRNA.level1.2.bin20.gz -out profile.AT.exon.lncRNA.level1.2.bin20.png --plotTitle "Exons" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

- reference-point

plotProfile -m matrix.exon.start.lncRNA.level1.2.bin20.gz -out profile.exon.start.lncRNA.level1.2.bin20.png --plotTitle "Exon" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

plotProfile -m matrix.exon.stop.lncRNA.level1.2.bin20.gz -out profile.exon.stop.lncRNA.level1.2.bin20.png --plotTitle "Level 1" -T "L1 motif density" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "L1 motif"

- AT%

plotProfile -m matrix.AT.exon.start.lncRNA.level1.2.bin20.gz -out profile.AT.exon.start.lncRNA.level1.2.bin20.png --plotTitle "Exon" -T "Percentage of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

plotProfile -m matrix.AT.exon.stop.lncRNA.level1.2.bin20.gz -out profile.AT.exon.stop.lncRNA.level1.2.bin20.png --plotTitle "Level 1" -T "Percenatge of AT" --colors blue --samplesLabel "Transcript type : lncRNA" --regionsLabel "AT%"

