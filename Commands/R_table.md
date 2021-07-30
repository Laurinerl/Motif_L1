# Visualisation R

## tables

gene_intergen_1 <- rbind(table_gene_annotation_occur_level1, table_intergenic_level1)

gene_intergen_1 <- gene_intergen_1 %>% mutate(density = occur/length)
>addition of a density column

gene_intergen_1_2 <- rbind(table_gene_annotation_occur_level1_2, table_intergenic_level1_2)

gene_intergen_1_2 <- gene_intergen_1_2 %>% mutate(density = occur/length)

coding_level1 <- rbind(table_transcript_exon_level1_coding_occur, table_intron_level1_coding_occur)

coding_level1 <- coding_level1 %>% mutate(density = occur/length)

lncRNA_level1 <- rbind(table_transcript_exon_level1_lncRNA_occur, table_intron_level1_lncRNA_occur)

lncRNA_level1 <- lncRNA_level1 %>% mutate(density = occur/length)

lncRNA_level1_2 <- rbind(table_transcript_exon_level1_2_lncRNA_occur, table_intron_level1_2_lncRNA_occur)

lncRNA_level1_2 <- lncRNA_level1_2 %>% mutate(density = occur/length)

level1_2 <- rbind(coding_level1_2, lncRNA_level1_2)

level1 <- rbind(coding_level1, lncRNA_level1)

## Description of the data 

-	Number of each feature 

ggplot(level1_2 %>% filter(length > 50, feature =="exon" | feature == "intron")) + geom_bar(aes(x = feature))

ggplot(level1 %>% filter(length > 50, feature =="exon" | feature == "intron")) + geom_bar(aes(x = feature))

-	Number in base of each feature

ggplot(level1_2 %>% filter(length > 50, feature =="exon" | feature == "intron")) + geom_col(aes(x = feature, y= length))
ggplot(level1 %>% filter(length > 50, feature =="exon" | feature == "intron")) + geom_col(aes(x = feature, y= length))

-	Number of each feature depending on the level 

ggplot(level1_2 %>% filter(length > 50)) + geom_bar(aes(x = feature, fill = level))
ggplot(gene_intergen_1_2 %>% filter(length > 50, feature == "gene")) + geom_bar(aes(x = feature, fill = level))

-	Density of features as a function of length 

ggplot(level1 %>% filter(feature == "exon", length > 50)) + geom_point(aes(x = length, y = density), alpha=0.2)
ggplot(level1_2 %>% filter(feature == "exon", length > 50)) + geom_point(aes(x = length, y = density), alpha=0.2)

## Testing the density of genes in L1 motifs as a function of their importance

-	motif enrichment by genes that takes into account the length of the genes

ggplot(gene_intergen_1 %>% filter(length > 50, feature == "gene")) + geom_density(aes(x = density))
ggplot(gene_intergen_1_2 %>% filter(length > 50, feature == "gene")) + geom_density(aes(x = density))

## Observe differences in density across regions

-	box_plot to compare density between intron/exon 
ggplot(level1_2 %>% filter(length > 50, feature =="intron" | feature == "exon")) +  geom_boxplot(aes(x = feature, y = density))
ggplot(level1 %>% filter(length > 50, feature =="intron" | feature == "exon")) +  geom_boxplot(aes(x = feature, y = density))

-	box_plot to compare density between gene/intergenic region 

ggplot(gene_intergen_1 %>% filter(length > 50)) +  geom_boxplot(aes(x = feature, y = density))
ggplot(gene_intergen_1_2 %>% filter(length > 50)) +  geom_boxplot(aes(x = feature, y = density))

-	boxplot to compare density between gene_type in exons

ggplot(level1 %>% filter(length > 50, feature == "exon")) +  geom_boxplot(aes(x = feature, y = density, color = gene_type))
ggplot(level1_2 %>% filter(length > 50, feature == "exon")) +  geom_boxplot(aes(x = feature, y = density, color = gene_type))

-	boxplot to compare density between gene_type in genes

ggplot(gene_intergen_1 %>% filter(feature == "gene",gene_type == "protein_coding"| gene_type =="lncRNA"| gene_type =="processed_pseudogene"| gene_type =="unprocessed_pseudogene")) + geom_boxplot(aes(x = feature, y = density, color = gene_type))


