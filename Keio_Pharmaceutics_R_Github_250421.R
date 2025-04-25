# This file describes the operation in RStudio following the operation in the terminal described in "Keio_Pharmaceutics_Terminal_Github_250421.txt".
# Versions (Including all dependent packages)

# Mac OS: Sequoia 15.3.2
# R: 4.4.2
# RStudio: 2024.09.1+394
# cowplot: 1.1.3
# ggplot2: 3.5.2
# edgeR: 4.4.2
# limma: 3.62.2
# tximport: 1.34.0
# clusterProfiler: 4.14.6
# dplyr: 1.1.4
# data.table: 1.17.0
# org.Rn.eg.db: 3.20.0
# org.Mm.eg.db: 3.20.0
# org.Hs.eg.db: 3.20.0
# AnnotationDbi: 1.68.0
# IRanges: 2.40.1
# S4Vectors: 0.44.0
# Biobase: 2.66.0
# BiocGenerics: 0.52.0
# others(Packages that depend on R)



# Data Explanation (same as "Keio_Pharmaceutics_Terminal_Github_250421.txt")

# mouse_LZ: paired-end mouse labyrinth RNA-seq data
# rat_LZ: paired-end rat labyrinth RNA-seq data
# SRR15514355: single-end human syncytiotrophoblast RNA-seq data (GSE182381)
# SRR15514356: single-end human syncytiotrophoblast RNA-seq data (GSE182381)
# SRR15514387: paired-end human syncytiotrophoblast RNA-seq data (GSE182381)



# We do not include any places of files and change of directories.
# First, the correspondence between Refseq IDs and gene names was entered.

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(data.table)
library(dplyr)
library(clusterProfiler)
RefseqID_to_gene <- function(data.dir, filename, species = c("mouse", "human", "rat")){
  
  OrgDb <- switch(species,
                  "mouse" = "org.Mm.eg.db",
                  "human" = "org.Hs.eg.db",
                  "rat" = "org.Rn.eg.db"
  )
  
  setwd(data.dir)
  object <- fread(filename)
  Refseq <- gsub(object$Name, pattern = "\\.[0-9]+", replacement = "")
  Symbol <- bitr(Refseq, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = OrgDb, drop = F)$SYMBOL
  
  mid <- object %>%
    dplyr::mutate(Refseq_Symbol = paste(Refseq, Symbol, sep = "~")) %>%
    dplyr::select(-Name)
  
  final <- mid %>%
    dplyr::mutate(Refseq = gsub("~[^~]*$", "", mid$Refseq_Symbol)) %>%
    dplyr::mutate(Symbol = gsub("^[^~]*~", "", mid$Refseq_Symbol)) %>%
    dplyr::select(Refseq_Symbol, Refseq, Symbol, everything())
  
  return(final)
}

mouse_LZ_quant <- RefseqID_to_gene(data.dir = "salmon/output", filename = "mouse_LZ_quant.sf", species = "mouse")
rat_LZ_quant <- RefseqID_to_gene(data.dir = "salmon/output", filename = "rat_LZ_quant.sf", species = "rat")
SRR15514355_quant <- RefseqID_to_gene(data.dir = "salmon/output", filename = "SRR15514355_quant.sf", species = "human")
SRR15514356_quant <- RefseqID_to_gene(data.dir = "salmon/output", filename = "SRR15514356_quant.sf", species = "human")
SRR15514387_quant <- RefseqID_to_gene(data.dir = "salmon/output", filename = "SRR15514387_quant.sf", species = "human")



# Next, only Refseq IDs whose prefixes began with NM were retained and the TPMs were re-calculated (did not use TPM this time).

NM_subset <- function(object){
  object <- object[grep("^NM_", object$Refseq),]
  rpk <- object$NumReads/(object$EffectiveLength/1000)
  object[["TPM"]] <- rpk/sum(rpk)*10^6
  return(object)
}

mouse_LZ_quant <- NM_subset(mouse_LZ_quant)
rat_LZ_quant <- NM_subset(rat_LZ_quant)
SRR15514355_quant <- NM_subset(SRR15514355_quant)
SRR15514356_quant <- NM_subset(SRR15514356_quant)
SRR15514387_quant <- NM_subset(SRR15514387_quant)



# Next, Using Tximport, count values for each Refseq ID were summarized by gene name.

mouse_Refseq <- mouse_LZ_quant$Refseq
mouse_Symbol <- mouse_LZ_quant$Symbol
rat_Refseq <- rat_LZ_quant$Refseq
rat_Symbol <- rat_LZ_quant$Symbol
human_Refseq <- SRR15514355_quant$Refseq
human_Symbol <- SRR15514355_quant$Symbol

shapedata_for_tximport <- function(object){
  object <- object[,c(-1,-3)]
  colnames(object)[1] <- "Name"
  return(object)
}

mouse_LZ_quant <- shapedata_for_tximport(mouse_LZ_quant)
rat_LZ_quant <- shapedata_for_tximport(rat_LZ_quant)
SRR15514355_quant <- shapedata_for_tximport(SRR15514355_quant)
SRR15514356_quant <- shapedata_for_tximport(SRR15514356_quant)
SRR15514387_quant <- shapedata_for_tximport(SRR15514387_quant)

library(tximport)
txi_mouse <- tximport(c("mouse_LZ" = "mouse_LZ_quant_NM.tsv"), type = "salmon", txOut = FALSE, countsFromAbundance = "no", tx2gene = data.frame(tx_id = mouse_Refseq, gene_id = mouse_Symbol))
txi_rat <- tximport(c("rat_LZ" = "rat_LZ_quant_NM.tsv"), type = "salmon", txOut = FALSE, countsFromAbundance = "no", tx2gene = data.frame(tx_id = rat_Refseq, gene_id = rat_Symbol))
txi_human <- tximport(c("SRR15514355" = "SRR15514355_quant_NM.tsv", "SRR15514356" = "SRR15514356_quant_NM.tsv", "SRR15514387" = "SRR15514387_quant_NM.tsv"), type = "salmon", txOut = FALSE, countsFromAbundance = "no", tx2gene = data.frame(tx_id = human_Refseq, gene_id = human_Symbol))

mouse_counts <- txi_mouse$counts
rat_counts <- txi_rat$counts
human_counts <- txi_human$counts



# Next, rows were aligned.
# If there was no corresponding gene, the count was marked as 0 for convenience.
# However, in this case, it was not clear if no gene was detected or if there was no genetic information in the first place,
# so we output the cases where NA was used instead of 0 at the end (df_NA below is "Raw_counts.csv").

SortGenes <- function(genenames){
  genes <- genenames %>%
    toupper() %>%
    unique() %>% 
    sort()
  return(genes)
}

allgenes <- SortGenes(c(rownames(mouse_counts), rownames(rat_counts), rownames(human_counts)))
mousegenes <- SortGenes(rownames(mouse_counts))
ratgenes <- SortGenes(rownames(rat_counts))
humangenes <- SortGenes(rownames(human_counts))

df <- data.frame(Symbol = allgenes, mouse_LZ = rep(0, length(allgenes)), rat_LZ = rep(0, length(allgenes)), human_STB1 = rep(0, length(allgenes)), human_STB2 = rep(0, length(allgenes)), human_STB3 = rep(0, length(allgenes)))
df[df$Symbol %in% humangenes, "human_STB1"] <- human_counts[,1]
df[df$Symbol %in% humangenes, "human_STB2"] <- human_counts[,2]
df[df$Symbol %in% humangenes, "human_STB3"] <- human_counts[,3]
df[df$Symbol %in% mousegenes, "mouse_LZ"] <- mouse_counts[,1]
df[df$Symbol %in% ratgenes, "rat_LZ"] <- rat_counts[,1]

df_NA <- data.frame(Symbol = allgenes, mouse_LZ = rep(NA, length(allgenes)), rat_LZ = rep(NA, length(allgenes)), human_STB1 = rep(NA, length(allgenes)), human_STB2 = rep(NA, length(allgenes)), human_STB3 = rep(NA, length(allgenes)))
df_NA[df_NA$Symbol %in% humangenes, "human_STB1"] <- human_counts[,1]
df_NA[df_NA$Symbol %in% humangenes, "human_STB2"] <- human_counts[,2]
df_NA[df_NA$Symbol %in% humangenes, "human_STB3"] <- human_counts[,3]
df_NA[df_NA$Symbol %in% mousegenes, "mouse_LZ"] <- mouse_counts[,1]
df_NA[df_NA$Symbol %in% ratgenes, "rat_LZ"] <- rat_counts[,1]



# Finally, the TMM normalization coefficients were calculated and the CPM was calculated by using them.

library(edgeR)
group <- factor(c("mouse_LZ", "rat_LZ", "human_STB", "human_STB", "human_STB"), levels = c("mouse_LZ", "rat_LZ", "human_STB"))
y <- DGEList(counts = df, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
CPM <- cpm(y)
rownames(CPM) <- y$genes$Symbol
CPM <- as.data.frame(CPM)



# MAplot of SLC and ABC family genes

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
CPM_mean <- CPM
CPM_mean[["human_STB"]] <- apply(CPM[,3:5], 1, mean)
CPM_mean <- CPM_mean[,-3:-5]
CPM_mean[["log2FC_human_vs_Rat"]] <- log2((CPM_mean$rat_LZ + 1)/(CPM_mean$human_STB + 1))
CPM_mean[["log2FC_mouse_vs_Rat"]] <- log2((CPM_mean$rat_LZ + 1)/(CPM_mean$mouse_LZ + 1))
CPM_mean[["Average_human_and_Rat"]] <- (log2(CPM_mean$rat_LZ + 1) + log2(CPM_mean$human_STB + 1))/2
CPM_mean[["Average_mouse_and_Rat"]] <- (log2(CPM_mean$rat_LZ + 1) + log2(CPM_mean$mouse_LZ + 1))/2
CPM_mean <- CPM_mean[grep("^SLC|^ABC", rownames(CPM_mean)),]

ggplot(CPM_mean, aes(Average_human_and_Rat, log2FC_human_vs_Rat)) + 
  geom_point(pch = 1, cex = 2.5) +
  geom_hline(yintercept = c(0)) +
  xlab("Log2CPM Average") +
  ylab("Log2Fold Change") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2.5), limits = c(0,12)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2.5), limits = c(-11,11))

ggplot(CPM_mean, aes(Average_mouse_and_Rat, log2FC_mouse_vs_Rat)) + 
  geom_point(pch = 1, cex = 2.5) +
  geom_hline(yintercept = c(0)) +
  xlab("Log2CPM Average") +
  ylab("Log2Fold Change") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 15)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2.5), limits = c(0,12)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 2.5), limits = c(-11,11))