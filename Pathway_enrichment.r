#!/bin/bash

################ This script performs pathway enrichment analysis on RNA-Seq differentially expressed (DE) genes using clusterProfiler
################ Inputs: DE gene file
################ Outputs: overall figures

## request small number of resources
srun --mem=32GB --time=6:00:00 --pty --cpus-per-task=4 bash

module load rstudio/2022.07

R

## load required libraries
library(clusterProfiler)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggplot2)

## set output directory 
output.dir <- "/home/abalay/scratch/Cell_Lines/Transcriptomics/Pathway_enrichment_analysis"
if (!dir.exists(output.dir)){
    dir.create(output.dir)
}

setwd(output.dir)

## !!!!!!!!!!!! Specify full path to file with differentially expressed (DE) genes
de.file <- "/home/abalay/scratch/Cell_Lines/Transcriptomics/Count_matrices/Stringtie/DEG.List.EBV+KSHVvsEBVhuNSG.txt"
de <- read.table(de.file, sep = "\t", header=T, row.names=1) ## load DE file
dim(de)
# [1] 40085     6


## create an extra gene symbol column
de$symbol <- gsub(".*[|]", "", rownames(de))

################################ Pathway enrichment analysis
## extract genes with log2FC > 0.25 and adjusted p value < 0.05 for both EBV and EBV+KSHV infected huNSG mice separately
ebv.de <- de[de$log2FoldChange < -0.25 & de$padj < 0.05,]
dim(ebv.de)
# [1] 8742    7

ebv.kshv.de <- de[de$log2FoldChange > 0.25 & de$padj < 0.05,]
dim(ebv.kshv.de)
# [1] 5539    7

## run biological pathway analysis
ebv.out <- enrichGO(gene = ebv.de$symbol,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01)

## save the file
write.table(data.frame(ebv.out), file = "BP.analysis.EBVhuNSGmouse.txt", sep = "\t", col.names=NA)

ebv.kshv.out <- enrichGO(gene = ebv.kshv.de$symbol,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01)

## save the file
write.table(data.frame(ebv.kshv.out), file = "BP.analysis.EBVKSHVhuNSGmouse.txt", sep = "\t", col.names=NA)


### establish enriched pathways
### high in ebv+kshv group
### mitochondrial respiratory chain complex I assembly, ATP synthesis coupled electron transport, mitochondrial translation
### ubiquinone biosynthetic process, ribosome biogenesis, tRNA 5'-end processing, tRNA aminoacylation
### mitotic DNA replication, mitotic sister chromatid segregation, centromere complex assembly
### double-strand break repair via break-induced replication
### ERAD pathway, proteasome-mediated ubiquitin-dependent protein catabolic process
### retrograde protein transport, ER to cytosol, 
### cellular oxidant detoxification

### add all DE genes to list
de.list <- list()

### top DE genes in those pathways
ebv.kshv.out <- data.frame(ebv.kshv.out)
### mitochondrial respiratory chain complex I assembly
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "mitochondrial respiratory chain complex I assembly"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["mitochondrial respiratory chain complex I assembly"]] <- out[1:10]
out[1:10]
#  [1] "NDUFAF8" "NDUFAF1" "NDUFB3"  "FOXRED1" "NDUFS8"  "NDUFB6"  "NDUFB11"
#  [8] "NDUFB2"  "NDUFB7"  "NDUFC1"

### ATP synthesis coupled electron transport
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "ATP synthesis coupled electron transport"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["ATP synthesis coupled electron transport"]] <- out[1:10]
out[1:10]
#  [1] "ISCU"    "PINK1"   "NDUFAF1" "NDUFS6"  "NDUFB3"  "NDUFA4"  "CYCS"   
#  [8] "CCNB1"   "NDUFS8"  "UQCRFS1"

### mitochondrial translation
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "mitochondrial translation"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["mitochondrial translation"]] <- out[1:10]
out[1:10]
#  [1] "MPV17L2" "TRUB2"   "COA3"    "MRPL51"  "MRPS7"   "MRPL16"  "MRPS11" 
#  [8] "MRPS6"   "RCC1L"   "TUFM"

### ubiquinone biosynthetic process
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "ubiquinone biosynthetic process"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["ubiquinone biosynthetic process"]] <- out[1:8]
out[1:8]
# [1] "FDXR"   "UBIAD1" "COQ8A"  "PDSS1"  "COQ4"   "COQ7"   "COQ9"   "COQ6"

### tRNA 5'-end processing
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "tRNA 5'-end processing"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["tRNA 5'-end processing"]] <- out[1:9]
out[1:9]
# [1] "RPP21"    "RPP40"    "HSD17B10" "POP7"     "POP5"     "RPP14"    "POP4"    
# [8] "RPP30"    "PRORP"

### tRNA aminoacylation
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "tRNA aminoacylation"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["tRNA aminoacylation"]] <- out[1:10]
out[1:10]
#  [1] "FARS2"  "LARS2"  "GATB"   "VARS1"  "FARSA"  "LRRC47" "MARS2"  "FARSB" 
#  [9] "SARS2"  "DARS2"

### mitotic DNA replication
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "mitotic DNA replication"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["mitotic DNA replication"]] <- out[1:10]
out[1:10]
#  [1] "GINS3" "MCM4"  "CDC45" "RAD51" "LIG1"  "GINS1" "MCM2"  "MCM6"  "RTF2" 
# [10] "MCM3"

### mitotic sister chromatid segregation
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "mitotic sister chromatid segregation"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["mitotic sister chromatid segregation"]] <- out[1:10]
out[1:10]
#  [1] "CDC6"    "ZWINT"   "CDC20"   "MAD1L1"  "TUBG1"   "CDCA5"   "ANAPC15"
#  [8] "CCNB1"   "UBE2C"   "CDT1"

### centromere complex assembly
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "centromere complex assembly"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["centromere complex assembly"]] <- out[1:10]
out[1:10]
#  [1] "CENPS"  "CENPW"  "MIS18A" "CENPO"  "OIP5"   "CENPP"  "CENPX"  "H3-3A" 
#  [9] "H3-3B"  "CENPA"

# ### ERAD pathway
# genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "ERAD pathway"]
# genes <- strsplit(genes, "/")[[1]]

# out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
# out[1:10]
# #  [1] "DERL3"  "WFS1"   "SDF2L1" "CALR"   "AQP11"  "HSPA5"  "RNFT1"  "DERL2" 
# #  [9] "VCP"    "UFD1"

### ubiquitin-dependent ERAD pathway
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "ubiquitin-dependent ERAD pathway"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["ubiquitin-dependent ERAD pathway"]] <- out[1:10]
out[1:10]
#  [1] "DERL3"   "WFS1"    "CALR"    "HSPA5"   "DERL2"   "VCP"     "UFD1"   
#  [8] "SEC61B"  "HERPUD1" "RNF5"


### endoplasmic reticulum to cytosol transport
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "endoplasmic reticulum to cytosol transport"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["endoplasmic reticulum to cytosol transport"]] <- out[1:10]
out[1:10]
#  [1] "DERL3"   "DERL2"   "VCP"     "UFD1"    "SEC61B"  "HM13"    "HERPUD1"
#  [8] "ERLEC1"  "OS9"     "SELENOS"

### cellular oxidant detoxification
genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == "cellular oxidant detoxification"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["cellular oxidant detoxification"]] <- out[1:10]
out[1:10]
#  [1] "NQO1"    "MGST1"   "SRXN1"   "ALOX5AP" "PRDX4"   "GPX1"    "TXNRD1" 
#  [8] "PRDX1"   "UBIAD1"  "DHFR"

### high in ebv group
### antigen receptor-mediated signaling pathway (B cell receptor signaling pathway, T cell receptor signaling pathway)
### cytoplasmic pattern recognition receptor signaling pathway, cytokine production involved in immune response,
### defense response to virus, negative regulation of viral genome replication
### positive regulation of stress-activated MAPK cascade

ebv.out <- data.frame(ebv.out)

### B cell receptor signaling pathway
genes <- ebv.out$geneID[ebv.out$Description == "B cell receptor signaling pathway"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["B cell receptor signaling pathway"]] <- out[1:10]
out[1:10]
#  [1] "LYN"      "CD19"     "MS4A1"    "GCSAM"    "SLC39A10" "TEC"     
#  [7] "FCGR2B"   "BCL2"     "CD300A"   "LCK"

### T cell receptor signaling pathway
genes <- ebv.out$geneID[ebv.out$Description == "T cell receptor signaling pathway"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["T cell receptor signaling pathway"]] <- out[1:10]
out[1:10]
#  [1] "TRAF6"  "MAP3K7" "CSK"    "RBCK1"  "ELF1"   "DUSP22" "BTN3A1" "BTN3A2"
#  [9] "IKBKB"  "PLCG1"

### cytoplasmic pattern recognition receptor signaling pathway
genes <- ebv.out$geneID[ebv.out$Description == "cytoplasmic pattern recognition receptor signaling pathway"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["cytoplasmic pattern recognition receptor signaling pathway"]] <- out[1:10]
out[1:10]
#  [1] "LSM14A"  "TRAF6"   "MAP3K7"  "OTULIN"  "RIOK3"   "NFKBIA"  "ANKRD17"
#  [8] "PUM2"    "ALPK1"   "BIRC2"

# ### cytokine production involved in immune response
# genes <- ebv.out$geneID[ebv.out$Description == "cytokine production involved in immune response"]
# genes <- strsplit(genes, "/")[[1]]

# out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
# out[1:10]
# #  [1] "TRAF6"  "MAP3K7" "HLA-E"  "SLAMF1" "MIF"    "SMAD7"  "TRIM6"  "IRAK3" 
# #  [9] "CD55"   "SEMA7A"

### defense response to virus
genes <- ebv.out$geneID[ebv.out$Description == "defense response to virus"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["defense response to virus"]] <- out[1:10]
out[1:10]
#  [1] "LSM14A"  "DDX17"   "RIOK3"   "NLRC5"   "PPM1B"   "ANKRD17" "PUM2"   
#  [8] "TRIM5"   "MLKL"    "BIRC2"

### negative regulation of viral genome replication
genes <- ebv.out$geneID[ebv.out$Description == "negative regulation of viral genome replication"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["negative regulation of viral genome replication"]] <- out[1:10]
out[1:10]
#  [1] "SETDB1"   "SHFL"     "MPHOSPH8" "RNASEL"   "TRIM6"    "PLSCR1"  
#  [7] "TASOR"    "IFI16"    "OAS1"     "OAS2"

### positive regulation of stress-activated MAPK cascade
genes <- ebv.out$geneID[ebv.out$Description == "positive regulation of stress-activated MAPK cascade"]
genes <- strsplit(genes, "/")[[1]]

out <- de[which(de$symbol %in% genes),] %>% arrange(desc(log2FoldChange)) %>% filter(padj < 0.05) %>% pull(symbol)
de.list[["positive regulation of stress-activated MAPK cascade"]] <- out[1:10]
out[1:10]
#  [1] "TRAF6"   "KLHDC10" "AGER"    "NOD1"    "SPI1"    "DUSP22"  "SLAMF1" 
#  [8] "DUSP19"  "SEMA4C"  "TAOK3"

### put genes into single vector
genes.plot <- c()
for (i in names(de.list)){
    genes.plot <- append(genes.plot, de.list[[i]])
}

### remove duplicate genes
genes.plot <- genes.plot[!duplicated(genes.plot)]
length(genes.plot)
# [1] 157


####################### Heatmap of differentially expressed genes across samples in group 1, 2 and 3
###############
###############
### load normalized gene counts with both discovery + validation samples


### take the genes from enriched pathways and samples from the infected huNSG mouse
genes.plot <- c("CELSR3", "PCDHA6", "PCDHA2", "PCDHB2", "PCDHB3", "PCDHB9", "PCDHGA10", "PCDHGC3", "AMIGO2", "L1CAM",
    "PDPN", "CTSK", "FAP", "COL3A1", "COL5A2", "FSCN1", "PRDX4", "HAS2", "COL5A1", "COL5A2", "SERPINH1", "MMP1", "MMP2", "LOXL4", "POSTN",
    "MMP9", "THBS3", "HOXD8", "MSX1", "BMP6", "SFRP1", "WNT8A", "TSKU", "HOXC9", "ALX1", "CHST11", "NOG", "ACP5",
    "HK2", "IL1B", "CXCL8", "VEGFC", "FGF1", "HGF", "CYBB", "MDK", "ADAM12", "HIF1A", "GRN", "PRKCA", "TBXA2R", "LRG1", "C5AR1",
    "FCGR1A", "SLC11A1", "TAP2", "PSMB8", "TAP1", "RAB32", "CTSD", "LGMN", "PDIA3", "PYCARD", "CALR",
    "ISG15", "IFIH1", "CD14", "OAS3", "OAS2", "PLCG2", "PTPRS", "LILRB1",
    )
temp.counts <- norm.counts[norm.counts$symbols %in% genes.plot, hunsg.cols]
dim(temp.counts)
# [1] 157  14

rownames(temp.counts) <- gsub(".*[|]", "", rownames(temp.counts))

temp.counts$match <- match(rownames(temp.counts), genes.plot)

temp.counts <- temp.counts %>% arrange(match)
temp.counts$match <- NULL

#### min-max normalize count matrix
normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

temp.counts <- as.data.frame(t(apply(temp.counts, 1, normalize)))

### select genes to show
genes.show <- c("NDUFB3", "NDUFS8", "NDUFB7", "NDUFC1", "CYCS", "UQCRFS1", "TRUB2", "COA3", "RCC1L",
    "COQ7", "COQ9", "COQ6", "RPP40", "POP7", "PRORP", "LARS2", "MARS2", 
    "GINS3", "MCM2", "LIG1", "ZWINT", "TUBG1", "ANAPC15", "CENPS", "CENPP",
    "DERL3", "VCP", "UFD1", "HERPUD1", "RNF5", "SELENOS", "NQO1", "SRXN1", "PRDX1",
    "CD19", "MS4A1", "GCSAM", "TRAF6", "MAP3K7", "BTN3A1", "IKBKB", "NFKBIA", "ANKRD17",
    "NLRC5", "SETDB1", "MPHOSPH8", "PLSCR1", "TASOR", "IFI16", "OAS1", "OAS2",
    "KLHDC10", "SLAMF1", "TAOK3")

###  normalized data color palette
colors.minMax.viridis <- colorRamp2(
    seq(from=0, to=1, by=0.1),
    c("black",viridis(n=10)),
    space="RGB")

condition <- c(rep("EBV.huNSG", 6), rep("EBV+KSHV.huNSG",8))

col_levels = list(Group = c("EBV.huNSG" = "red", "EBV+KSHV.huNSG" = "blue")) ### set color levels of group

### Top annotation by group
ha_top <- HeatmapAnnotation(Group=condition, col = col_levels)

marker_anno <- rowAnnotation(
    marker=anno_mark(at=which(genes.plot %in% genes.show),
    labels=genes.show, labels_gp = gpar(fontsize = 7)))

hp <- Heatmap(temp.counts, 
    name = "MinMax", 
    col = colors.minMax.viridis,
    column_names_rot = 45, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    show_row_names = FALSE,
    row_names_side = "right", 
    column_names_side = "bottom", 
    row_names_gp = gpar(fontsize = 6), 
    column_names_gp = gpar(fontsize = 10), 
    top_annotation=ha_top,
    right_annotation=marker_anno, 
    border = TRUE)

png("Heatmap.pathwayenrichment.huNSGmice.051723.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()


####################### Heatmap of genes across pathways
###############
###############
### assign each gene to its pathway
gene.pathway <- list()
for (i in 1:length(genes.plot)){
    temp <- c()
    for (j in 1:length(names(de.list))){
        if (j < 13){
            genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == names(de.list)[j]]
            genes <- strsplit(genes, "/")[[1]]

            if (genes.plot[i] %in% genes){
                temp <- append(temp, "P")
            } else {
                temp <- append(temp, "NP")
            }
        } else {
            genes <- ebv.out$geneID[ebv.out$Description == names(de.list)[j]]
            genes <- strsplit(genes, "/")[[1]]

            if (genes.plot[i] %in% genes){
                temp <- append(temp, "P")
            } else {
                temp <- append(temp, "NP")
            }
        }
    }
    gene.pathway[[i]] <- temp
    names(gene.pathway)[i] <- genes.plot[i]
}

temp.df <- do.call(rbind, gene.pathway); colnames(temp.df) <- names(de.list); rownames(temp.df) <- genes.plot


### select genes to show
genes.show <- c("NDUFB3", "NDUFS8", "NDUFB7", "NDUFC1", "CYCS", "UQCRFS1", "TRUB2", "COA3", "RCC1L",
    "COQ7", "COQ9", "COQ6", "RPP40", "POP7", "PRORP", "LARS2", "MARS2", 
    "GINS3", "MCM2", "LIG1", "ZWINT", "TUBG1", "ANAPC15", "CENPS", "CENPP",
    "DERL3", "VCP", "UFD1", "HERPUD1", "RNF5", "SELENOS", "NQO1", "SRXN1", "PRDX1",
    "CD19", "MS4A1", "GCSAM", "TRAF6", "MAP3K7", "BTN3A1", "IKBKB", "NFKBIA", "ANKRD17",
    "NLRC5", "SETDB1", "MPHOSPH8", "PLSCR1", "TASOR", "IFI16", "OAS1", "OAS2",
    "KLHDC10", "SLAMF1", "TAOK3")

###  normalized data color palette
colors <- c("P" = "tomato1", "NP" = "lavenderblush")

# colors.minMax.viridis <- colorRamp2(
#     seq(from=0, to=1, by=0.1),
#     c("black",viridis(n=10)),
#     space="RGB")


marker_anno <- rowAnnotation(
    marker=anno_mark(at=which(genes.plot %in% genes.show),
    labels=genes.show, labels_gp = gpar(fontsize = 7)))

hp1 <- Heatmap(temp.df, 
    name = "Status", 
    col = colors,
    column_names_rot = 45, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    show_row_names = FALSE,
    row_names_side = "left", 
    column_names_side = "bottom", 
    row_names_gp = gpar(fontsize = 6), 
    column_names_gp = gpar(fontsize = 7),
    rect_gp = gpar(col = "white", lwd = 2),
    right_annotation=marker_anno, 
    border = TRUE,
    use_raster=FALSE)

png("Heatmap.pathwayanno.huNSGmice.051723.png", res=200, unit="in", height=8, width=5)
draw(hp1)
dev.off()


############################## Double-sided barplot
##############
##############
### create dataframe with condition, counts per pathway and pathway info
### count number of genes in each enriched pathway
counts.pathway <- c()
for (j in 1:length(names(de.list))){
    if (j < 13){
        genes <- ebv.kshv.out$geneID[ebv.kshv.out$Description == names(de.list)[j]]
        genes <- strsplit(genes, "/")[[1]]
        counts.pathway <- append(counts.pathway, length(genes))

    } else {
        genes <- ebv.out$geneID[ebv.out$Description == names(de.list)[j]]
        genes <- strsplit(genes, "/")[[1]]
        counts.pathway <- append(counts.pathway, -length(genes))
    }
}

bar.df <- data.frame(Pathway = names(de.list), Condition = c(rep("EBV+KSHV.huNSG",12), rep("EBV.huNSG", 6)), 
    Counts = counts.pathway)

### factorize both pathway and condition columns
bar.df$Pathway <- factor(bar.df$Pathway, levels = bar.df$Pathway)
bar.df$Condition <- factor(bar.df$Condition, levels = c("EBV.huNSG", "EBV+KSHV.huNSG"))

png("Barplot.pathwaystats.huNSGmice.051723.png", res=200, unit="in", height=8, width=11)
ggplot(bar.df,aes(x=Pathway, y=Counts, fill=Condition))+
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("EBV.huNSG" = "red", "EBV+KSHV.huNSG" = "blue")) +
    scale_y_continuous(breaks = seq(-60, 60, 10)) +
    theme_classic() +
    coord_flip()
dev.off()


