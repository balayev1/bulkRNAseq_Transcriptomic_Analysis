#!/bin/bash

################ This script performs differential expression (DE) analysis on RNA-Seq count tables
################ Inputs: Raw count table file + annotation file
################ Outputs: DE table + overall figures

## request small number of resources
srun --mem=32GB --time=6:00:00 --pty --cpus-per-task=4 bash

## Activate conda environment
module load anaconda3
source activate r4_env

R

## load required libraries
library(readxl)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(scatterplot3d)
library(ComplexHeatmap)
library(factoextra)
library(circlize)
library(viridis)

## set output directory 
output.dir <- "/home/abalay/scratch/Cell_Lines/Transcriptomics/Count_matrices/Stringtie"
if (!dir.exists(output.dir)){
    dir.create(output.dir)
}

## !!!!!!!!!!!! Specify full path to raw count table file
counts.file <- "/home/abalay/scratch/Cell_Lines/Transcriptomics/Count_matrices/Stringtie/Cell_Lines_gene_count.csv"
counts <- read.table(counts.file, sep = ",", header=T, row.names=1) ## load counts file
dim(counts)
# [1] 60708    79

## !!!!!!!!!!!! Specify full path to annotation file
anno.file <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/EBV-transformed_lymphoblast_cell_lines.xlsx"
anno <- as.data.frame(read_excel(anno.file, sheet="huNSG mice samples"))[1:79,] ## load annotation file as data frame

## replace any space in the column names with .
colnames(anno) = gsub(" ",".",colnames(anno))

## assign SRA.ID column as the row names
rownames(anno) <- anno[,"SRA-ID"]

## check if all row names in annotation file match column names in counts file
all(rownames(anno) == colnames(counts))
# [1] FALSE

## if FALSE match the names
colnames(counts) <- gsub("_2", "", colnames(counts))
all(rownames(anno) == colnames(counts))
# [1] TRUE

## check for duplicate genes
counts$means <- rowMeans(counts)
counts$symbols <- gsub(".*[|]", "", rownames(counts)) ## create gene symbol column
duplicated.genes <- counts[duplicated(counts$symbols),] ## check for duplicates
table(duplicated.genes$symbols)
#   5_8S_rRNA     5S_rRNA         7SK  AC005476.2  AC008686.1  AC010618.3 
#           5           8           6           1           1           1 
#  AC068587.4  AC073218.1  AC091132.4  AC099520.1  AC104041.1  AC106795.1 
#           1           1           1           1           1           1 
#  AC138627.1  AC145350.2  AJ271736.1     AKAP17A  AL110114.1  AL138756.1 
#           1           1           1           1           1           1 
#  AL139317.3  AL157388.1  AL162253.2  AL355076.2  AL390957.1  AL672277.1 
#           1           1           1           1           1           1 
#  AL683807.1  AL683807.2  AL732314.4  AL732314.6  AL954722.1      AMD1P2 
#           1           1           1           1           1           1 
#  AP001021.1  AP001107.9        ASMT       ASMTL   ASMTL-AS1      ATP11A 
#           1           1           1           1           1           1 
#        CAST      CCDC13     CCDC192        CD99      CD99P1        CDH3 
#           1           1           1           1           1           1 
#   CLCA4-AS1     CNTNAP2        COG7       CRLF2       CRLF3      CSF2RA 
#           1           1           1           1           1           1 
#    DDX11L16       DHRSX   DHRSX-IT1  DNAJC9-AS1      DPH3P2       ELFN2 
#           1           1           1           1           1           1 
#     ELOCP24       ENPP2       ERBB3       EXOC4    FABP5P13     FAM157C 
#           1           1           1           1           1           1 
#     FAM182A   GABARAPL3     GOLGA8M       GOLM2        GPX3        GRK4 
#           1           1           1           1           1           1 
#      GTPBP6    H2AZ1-DT     HERC2P7    HMGN2P46      IL1RAP       IL3RA 
#           1           1           1           1           1           1 
#        IL9R    IQCH-AS1 KBTBD11-OT1       KCNA2      KCNMA1    KRT18P53 
#           1           1           1           1           1           1 
#   LINC00102   LINC00106   LINC00299   LINC00342   LINC00484   LINC00486 
#           1           1           1           1           1           1 
#   LINC00685   LINC00882   LINC01088   LINC01115   LINC01238   LINC02203 
#           1           1           1           1           1           1 
#   LINC02245   LINC02256   LINC02885   MEF2C-AS2 Metazoa_SRP    MGC32805 
#           1           1           1           1         169           1 
#       MGST3     MIR3690     MIR6089       MPRIP        MSRA       NDST1 
#           1           1           1           1           1           1 
#  NUTM2A-AS1         OAF       P2RY8        PBX1        PBX3       PDE3B 
#           1           1           1           1           1           1 
#       PIAS2     PIP4K2A      PLCXD1        PLD1       POLA1     PPP2R3B 
#           1           1           1           1           2           1 
#        PRPH       PWRN1        REST        RGS6   RNA5SP498        ROM1 
#           1           1           1           1           1           1 
#     RPL14P5       SCN2A      SCNN1B       SGMS2        SHOX     SLC19A1 
#           1           1           1           1           1           1 
#     SLC25A6     SNORA62     SNORA63     SNORA70     SNORA71     SNORA72 
#           1           6           6          26           2           7 
#     SNORA73     SNORA74     SNORA75    SNORD115    SNORD116     SNORD18 
#           3           3           6           1           2           1 
#     SNORD27     SNORD30     SNORD33     SNORD39     SNORD42     SNORD63 
#           1           1           1           3           1           2 
#     SNORD81   SOX21-AS1   STX17-AS1       SYNE2 TBC1D26-AS1        TLN2 
#           2           1           1           1           1           1 
#     TMSB15B        TNKS        TP53      TRPC6P          U1          U2 
#           1           1           1           1           2          18 
#          U3          U4          U6          U7          U8       VAMP7 
#          49           6          32           6          21           1 
#        VAPA       Vault    VCAN-AS1      WASH6P      WASIR1       WNT5A 
#           1           3           1           1           1           1 
#       Y_RNA       ZBED1     ZFYVE28     ZSCAN5A 
#         753           1           1           1 

## check the distribution of means of duplicate genes
# summary(counts$$means)
# summary(duplicated.genes$means)

counts$means <- NULL; counts$symbols <- NULL
## remove unannotated MSTRG genes from row names
length(intersect(rownames(counts), grep("|", rownames(counts), value=TRUE))) == dim(counts)[1]
# [1] TRUE
## no unannotated genes

## before creating DESeqDataset: add two columns Infection.agent/Origin + Cell.line/Replicate
anno$Origin.infection=paste(anno$`Infection.agent`, anno$`Origin`, sep=".")
anno$Cell.replicate=paste(anno$`Cell.line`, anno$`Replicate.#`, sep=".")

## put infection by cell origin as main source of variation
gene_dds <- DESeqDataSetFromMatrix(countData = counts, colData = anno, design = ~ Origin.infection)

## remove all genes with total counts < 10 across samples
genes.to.keep <- rowSums(counts(gene_dds)) > 10; gene_dds <- gene_dds[genes.to.keep,]
dim(gene_dds)
# [1] 40085    79

gene_dds$Origin.infection <- factor(gene_dds$Origin.infection, levels=unique(gene_dds$Origin.infection))
gene_dds$Cell.replicate <- factor(gene_dds$Cell.replicate, levels=unique(gene_dds$Cell.replicate))

## collapse technical replicates to one biological replicate per sample
gene.dds.collapsed = collapseReplicates(gene_dds, gene_dds$Cell.replicate)

## check if sum of counts for sample before collapse == counts for sample post collapse
matchFirstLevel <-  gene_dds$Cell.replicate == levels(gene_dds$Cell.replicate)[1]
all(rowSums(counts(gene_dds[,matchFirstLevel])) == counts(gene.dds.collapsed[,1]))
# [1] TRUE

## list of samples
colnames(gene.dds.collapsed)
#  [1] "E4-3.R1"   "E4-3.R2"   "E4-3.R3"   "E4-9.R1"   "E4-9.R2"   "E4-9.R3"  
#  [7] "EK4-11.R1" "EK4-11.R2" "EK4-11.R3" "EK3-11.R1" "EK3-11.R2" "EK3-11.R3"
# [13] "EK3-13.R1" "EK3-13.R2" "LCL1.R1"   "LCL1.R2"   "LCL1.R3"   "LCL2.R1"  
# [19] "LCL2.R2"   "LCL2.R3"   "LCL3.R1"   "LCL3.R2"   "LCL3.R3"   "BCBL1.R1" 
# [25] "AP2.R1"    "AP3.R1"    "AP5.R1"    "HBL6.R1"

## save the whole collapsed object as Robj
save(gene.dds.collapsed, file = file.path(output.dir, "Collapsed.counts.huNSGcelllines.Robj"))

## collapsed count matrix
gene.dds.collapsed.count <- data.frame(counts(gene.dds.collapsed))

## save count matrix as tab-delimited text file
write.table(gene.dds.collapsed.count, file.path(output.dir, "Collapsed.raw.counts.cell.lines.txt"), sep="\t", col.names=NA)

## create a folder to put the figures
figures.folder <- file.path(output.dir, "figures")
if (!dir.exists(figures.folder)){
    dir.create(figures.folder)
}

## change working directory to figures folder
setwd(figures.folder)

## density plot of raw counts across samples
Sample <- list(); Counts <- list() 
for (i in 1:length(colnames(gene.dds.collapsed.count))){
    Sample[[i]] <- data.frame(sample = rep(colnames(gene.dds.collapsed.count)[i], dim(gene.dds.collapsed.count)[1]))
    Counts[[i]] <- data.frame(count = gene.dds.collapsed.count[,i])
}
sample <- as.data.frame(do.call(rbind, Sample))
ncounts <- as.data.frame(do.call(rbind, Counts))

gene.dds.collapsed.countm <- cbind(sample = sample, counts = ncounts)

colnames(gene.dds.collapsed.countm) <- c("Sample", "Counts")

png("Density.plot.gene.count.051223.png", res=200, unit="in", height=8, width=11)
ggplot(gene.dds.collapsed.countm, aes(x=Counts, color=Sample, add="median")) + geom_density() +
    scale_x_log10() + xlab("Counts")
dev.off()

## fit NGB model to gene counts from different factors of infection+origin (i.e. EBV:huNSG, EBV+KSHV:huNSG and et.c.)
gene.dds.collapsed=DESeq(gene.dds.collapsed, test="Wald")

resultsNames(gene.dds.collapsed)
# [1] "Intercept"                                   
# [2] "Origin.infection_EBV.KSHV.huNSG_vs_EBV.huNSG"
# [3] "Origin.infection_EBV.CL_vs_EBV.huNSG"        
# [4] "Origin.infection_KSHV.CL_vs_EBV.huNSG"       
# [5] "Origin.infection_EBV.KSHV.CL_vs_EBV.huNSG"

ref.group = "EBV.huNSG"
res = results(gene.dds.collapsed, contrast=c("Origin.infection", "EBV.KSHV.huNSG", ref.group), 
    independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)

# Save DEG list for EBV+KSHV.huNSG vs EBV.huNSG as tab-delimited text file
write.table(data.frame(res), file=file.path(output.dir, "DEG.List.EBV+KSHVvsEBVhuNSG.txt"), sep = "\t", col.names=NA)

################################## Analysis ########################################
#####################
##### Normalize the raw count data 
### Using VST transformation
gene.dds.collapsed.norm <- vst(gene.dds.collapsed, blind=TRUE)
write.table(assay(gene.dds.collapsed.norm), file.path(output.dir, "Collapsed.vstnorm.counts.cell.lines.txt"), sep="\t", col.names=NA)

#### Compute PCA
res.pca <- prcomp(t(assay(gene.dds.collapsed.norm)), scale = FALSE)

### Scree plot for the variance explained by eigenvectors
png("Screeplot.PCA.cell.lines.051223.png", res=200, unit="in", height=8, width=11)
fviz_eig(res.pca)
dev.off()

### 3D PCA PLOT
variance.explained <- (res.pca$sdev)^2/(sum((res.pca$sdev)^2))*100

png("3DPCA.plot.cell.lines.051223.png", res=200, unit="in", height=8, width=11)
pch.shapes <- c(4, 15, 16, 17, 18); pch.shapes <- pch.shapes[gene.dds.collapsed$Origin.infection]
colors <- c("red", "blue", "green", "orange", "purple"); colors <- colors[gene.dds.collapsed$Origin.infection]
s3d <- scatterplot3d(res.pca$x[,1:3], main = "3D PCA Plot for Cell Lines", 
    xlab = paste("PC1 (", round(variance.explained[1],1), " %)", sep=""),
    ylab = paste("PC2 (", round(variance.explained[2],1), " %)", sep=""),
    zlab = paste("PC3 (", round(variance.explained[3],1), " %)", sep=""),
    pch = pch.shapes,
    color=colors,
    grid=TRUE,
    box=TRUE)
legend("bottom", pch = c(4, 15, 16, 17, 18), col = c("red", "blue", "green", "orange", "purple"),
    legend = unique(gene.dds.collapsed$Origin.infection), xpd = TRUE, horiz = TRUE, inset = -0.15)
dev.off()


### Heatmap: expression of gene-editing enzymes in all samples
norm.gene.count <- data.frame(assay(gene.dds.collapsed.norm))

#### min-max normalize count matrix
normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

norm.gene.count$symbols <- gsub(".*[|]", "", rownames(norm.gene.count))

###### APOBEC3A/B/C/D/F/G/H, AICDA, POLH, ADAR1, ADARB1, ADARB2
#### Note: APOBEC3E NOT IN EXPRESSION MATRIX
editing.genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",  "APOBEC3F", "APOBEC3G", "APOBEC3H",
    "AICDA", "POLH", "ADAR", "ADARB1", "ADARB2")

temp.count <- norm.gene.count[norm.gene.count$symbols %in% editing.genes, ]
rownames(temp.count) <- temp.count$symbols; temp.count$symbols <- NULL

temp.count <- as.data.frame(apply(temp.count, 1, normalize))

### factorize gene order
col.order <- c("ADAR", "AICDA", "APOBEC3A", "APOBEC3B", "APOBEC3C", 
    "APOBEC3G", "APOBEC3D", "APOBEC3F", "APOBEC3H", "POLH", "ADARB1", "ADARB2")

temp.count <- temp.count[, col.order]

#### Heatmap
###  normalized data color palette
colors.minMax.viridis <- colorRamp2(
    seq(from=0, to=1, by=0.1),
    c("black",viridis(n=10)),
    space="RGB")

# col_levels = list(Group = c("EBV-huNSG" = "green", "EBV+KSHV-huNSG" = "blue", "EBV-LCL" = "red", 
#     "EBV+KSHV-Cell.Line" = "purple", "KSHV-Cell.Line" = "yellow")) ### set color levels of group

col_levels = list(Group = c("EBV.huNSG" = "red", "EBV+KSHV.huNSG" = "blue", "EBV.CL" = "green", 
    "EBV+KSHV.CL" = "orange", "KSHV.CL" = "purple")) ### set color levels of group

### Left annotation by group
ha_right <- rowAnnotation(Group=gene.dds.collapsed$Origin.infection, col = col_levels)

hp <- Heatmap(temp.count, 
    name = "z-score", 
    col = colors.minMax.viridis, 
    row_title = "Samples",
    row_title_rot = 90, 
    column_names_rot = 45, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    row_names_side = "left", 
    column_names_side = "top", 
    row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), 
    row_split = gene.dds.collapsed$Origin.infection, 
    column_split=rep(c("A", "B", "C"), 
    c(2,4,6)), 
    right_annotation=ha_right, border = TRUE)

png("Heatmap.editing_enzyme.expression_allsamples.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()

### Heatmap: expression of gene-editing enzymes in samples from huNSG mouse
editing.genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",  "APOBEC3F", "APOBEC3G", "APOBEC3H",
    "AICDA", "POLH", "ADAR", "ADARB1", "ADARB2")

temp.count <- norm.gene.count[norm.gene.count$symbols %in% editing.genes, ]
rownames(temp.count) <- temp.count$symbols; temp.count$symbols <- NULL

hunsg.cols <- c(grep("E4", colnames(temp.count)), grep("EK4", colnames(temp.count)), grep("EK3", colnames(temp.count)))
temp.count <- temp.count[, hunsg.cols]

temp.count <- as.data.frame(apply(temp.count, 1, normalize))

### factorize gene order
col.order <- c("ADAR", "AICDA", "APOBEC3A", "APOBEC3B", "APOBEC3C", 
    "APOBEC3G", "APOBEC3D", "APOBEC3F", "APOBEC3H", "POLH", "ADARB1", "ADARB2")

temp.count <- temp.count[, col.order]

#### Heatmap
###  normalized data color palette
colors.minMax.viridis <- colorRamp2(
    seq(from=0, to=1, by=0.1),
    c("black",viridis(n=10)),
    space="RGB")

# col_levels = list(Group = c("EBV-huNSG" = "green", "EBV+KSHV-huNSG" = "blue", "EBV-LCL" = "red", 
#     "EBV+KSHV-Cell.Line" = "purple", "KSHV-Cell.Line" = "yellow")) ### set color levels of group

col_levels = list(Group = c("EBV.huNSG" = "red", "EBV+KSHV.huNSG" = "blue")) ### set color levels of group

### Left annotation by group
ha_right <- rowAnnotation(Group=gene.dds.collapsed$Origin.infection[1:14], col = col_levels)

hp <- Heatmap(temp.count, 
    name = "z-score", 
    col = colors.minMax.viridis, 
    row_title = "Samples",
    row_title_rot = 90, 
    column_names_rot = 45, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    row_names_side = "left", 
    column_names_side = "top", 
    row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), 
    row_split = gene.dds.collapsed$Origin.infection[1:14], 
    column_split=rep(c("A", "B", "C"), 
    c(2,4,6)), 
    right_annotation=ha_right, border = TRUE)

png("Heatmap.editing_enzyme.expression_hunsgsamples.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()



