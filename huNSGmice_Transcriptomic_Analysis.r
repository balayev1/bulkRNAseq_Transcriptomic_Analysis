library("ggpubr")
library("DESeq2")
library("factoextra")
library("readxl")
library("scatterplot3d")
library("ComplexHeatmap")
library("circlize")
library("dplyr")
library("GOplot")


## Load gene matrix of infected huNSG mice samples
parent.directory="/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/Count_matrices/Stringtie/"
gene_count=read.csv(paste(parent.directory,'Cell_Lines_gene_count.csv',sep=""), sep=",", header=T, row.names=1); colnames(gene_count)=gsub("_cutadapt_sorted", "", colnames(gene_count))
nrow(gene_count)
# [1] 60713

## Load annotation table 
reference.dir="/home/cluster/abalay/scratch/Pipelines_and_Files/Files_for_scripts"
anno.table=data.frame(read_excel(paste(reference.dir, "/EBV-transformed_lymphoblast_cell_lines.xlsx",sep=""), sheet="huNSG mice samples"))
anno.table$Origin.infection=paste(anno.table$`Infection.agent`, anno.table$`Origin`, sep=".") # create artificial column with infection agent + sample cellular origin
anno.table$Cell.replicate=paste(anno.table$`Cell.line`, anno.table$`Replicate..`, sep=".") # create artificial column with cell line + replicate no
colnames(anno.table) = gsub(" ",".",colnames(anno.table)); anno.table$Origin.infection=factor(anno.table$Origin.infection); anno.table$Cell.replicate=factor(anno.table$Cell.replicate)

# Subset rows with valid samples and remove the samples of no interest
anno.table = anno.table[1:79,]; rownames(anno.table) <- anno.table[,"SRA.ID"]
gene_count <- gene_count[, colnames(gene_count) %in% rownames(anno.table)]

# Make sure each row in annotation table correspond to each column in gene matrix
all(rownames(anno.table) == colnames(gene_count))
# [1] TRUE

## Remove genes with zero counts across samples
gene_count=gene_count[rowSums(gene_count[,])!=0,] # remove genes with zero counts across all samples
nrow(gene_count)
# [1] 46162

## Check for duplicate genes in gene matrix
gene_count$gene.symbols = gsub(".*[|]", "", rownames(gene_count))
duplicate.gene.count = gene_count[duplicated(gene_count$gene.symbols),]
table(duplicate.gene.count$gene.symbols)
# 5_8S_rRNA     5S_rRNA         7SK  AC005476.2  AC008686.1  AC010618.3 
#           2           5           6           1           1           1 
#  AC068587.4  AC091132.4  AC099520.1  AC104041.1  AC106795.1  AC145350.2 
#           1           1           1           1           1           1 
#  AJ271736.1     AKAP17A  AL110114.1  AL138756.1  AL139317.3  AL157388.1 
#           1           1           1           1           1           1 
#  AL162253.2  AL355076.2  AL390957.1  AL672277.1  AL683807.1  AL683807.2 
#           1           1           1           1           1           1 
#  AL732314.4  AL732314.6      AMD1P2  AP001107.9        ASMT       ASMTL 
#           1           1           1           1           1           1 
#   ASMTL-AS1      ATP11A        CAST      CCDC13        CD99      CD99P1 
#           1           1           1           1           1           1 
#        CDH3   CLCA4-AS1        COG7       CRLF2       CRLF3      CSF2RA 
#           1           1           1           1           1           1 
#    DDX11L16       DHRSX   DHRSX-IT1  DNAJC9-AS1       ELFN2     ELOCP24 
#           1           1           1           1           1           1 
#       ENPP2       EXOC4     FAM157C     FAM182A   GABARAPL3       GOLM2 
#           1           1           1           1           1           1 
#        GPX3        GRK4      GTPBP6    H2AZ1-DT     HERC2P7    HMGN2P46 
#           1           1           1           1           1           1 
#      IL1RAP       IL3RA        IL9R    IQCH-AS1       KCNA2   LINC00102 
#           1           1           1           1           1           1 
#   LINC00106   LINC00299   LINC00342   LINC00484   LINC00685   LINC01088 
#           1           1           1           1           1           1 
#   LINC01238   LINC02203   LINC02245   LINC02256   LINC02885 Metazoa_SRP 
#           1           1           1           1           1         124 
#       MGST3     MIR6089       MPRIP        MSRA       NDST1  NUTM2A-AS1 
#           1           1           1           1           1           1 
#         OAF       P2RY8        PBX1        PBX3       PIAS2     PIP4K2A 
#           1           1           1           1           1           1 
#      PLCXD1        PLD1       POLA1     PPP2R3B        REST        RGS6 
#           1           1           2           1           1           1 
#        ROM1       SCN2A      SCNN1B       SGMS2     SLC19A1     SLC25A6 
#           1           1           1           1           1           1 
#     SNORA62     SNORA63     SNORA70     SNORA72     SNORA73     SNORA74 
#           4           2          14           4           3           2 
#     SNORA75     SNORD18     SNORD27     SNORD63   SOX21-AS1       SYNE2 
#           4           1           1           1           1           1 
#        TLN2     TMSB15B        TNKS        TP53      TRPC6P          U1 
#           1           1           1           1           1           2 
#          U2          U3          U4          U6          U7          U8 
#          14          15           5          13           2           7 
#       VAMP7        VAPA       Vault      WASH6P      WASIR1       Y_RNA 
#           1           1           2           1           1         475 
#       ZBED1     ZFYVE28     ZSCAN5A 
#           1           1           1  

## Delete rRNA genes + genes with multiple loci + small nucleolar RNA genes
genes.to.remove = c("5_8S_rRNA", "5S_rRNA", "U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8", "Y_RNA", "Metazoa_SRP") # delete rRNA genes and genes with multple loci
gene_count = gene_count[!gene_count$gene.symbols %in% genes.to.remove,]; gene_count = gene_count[!grepl("SNORA", gene_count$gene.symbols) & !grepl("SNORD", gene_count$gene.symbols),]
duplicate.gene.count = gene_count[duplicated(gene_count$gene.symbols),]

## Select gene copy with highest average expression across every sample
indices.to.be.removed = c()
for (gene in unique(duplicate.gene.count$gene.symbols)){
  gene.indices = which(gene_count$gene.symbols==gene)
  expression=c()
  for (index in gene.indices){
    expression=append(expression, mean(as.numeric(gene_count[index, 1:(ncol(gene_count)-1)])))
  }
  if (length(which(expression==max(expression)))>1){
    indices.to.be.removed = append(indices.to.be.removed, gene.indices[which(expression==max(expression))[2:length(expression)]])
  }
  if (length(which(expression==max(expression)))==1){
    indices.to.be.removed = append(indices.to.be.removed, gene.indices[-which(expression==max(expression))])
  }
}
gene_count = gene_count[-indices.to.be.removed,]
duplicate.gene.count = gene_count[duplicated(gene_count$gene.symbols),]
table(duplicate.gene.count$gene.symbols)
# table of extent 0 >

# Set rows to gene symbols
rownames(gene_count)=gene_count$gene.symbols; gene_count$gene.symbols=NULL
nrow(gene_count)
# [1] 44993

# Put infection by cell origin as main source of variation
gene_dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = anno.table, design = ~ Origin.infection)

# Remove all genes with total counts < 10 across samples
genes.to.keep = rowSums(counts(gene_dds)) > 10; gene_dds <- gene_dds[genes.to.keep,]
dim(gene_dds)
# [1] 37082    79

gene_dds$Origin.infection = factor(gene_dds$Origin.infection, levels=unique(gene_dds$Origin.infection))

# Collapse technical replicates to one biological replicate per sample
gene.dds.collapsed = collapseReplicates(gene_dds, anno.table$Cell.replicate)

# Check if sum of counts for sample before collapse == counts for sample post collapse
matchFirstLevel <-  anno.table$Cell.replicate == levels(anno.table$Cell.replicate)[1]
all(rowSums(counts(gene_dds[,matchFirstLevel])) == counts(gene.dds.collapsed[,1]))
# [1] TRUE

# List of samples
colnames(gene.dds.collapsed)
# [1] "AP2.R1"    "AP3.R1"    "AP5.R1"    "BCBL1.R1"  "E4-3.R1"   "E4-3.R2"  
#  [7] "E4-3.R3"   "E4-9.R1"   "E4-9.R2"   "E4-9.R3"   "EK3-11.R1" "EK3-11.R2"
# [13] "EK3-11.R3" "EK3-13.R1" "EK3-13.R2" "EK4-11.R1" "EK4-11.R2" "EK4-11.R3"
# [19] "HBL6.R1"   "LCL1.R1"   "LCL1.R2"   "LCL1.R3"   "LCL2.R1"   "LCL2.R2"  
# [25] "LCL2.R3"   "LCL3.R1"   "LCL3.R2"   "LCL3.R3" 

gene.dds.collapsed.count <- data.frame(counts(gene.dds.collapsed))

# Save count matrix as txt file
write.table(gene.dds.collapsed.count, file.path(parent.directory, "Raw.Counts.Cell.Lines.txt"), sep="\t", col.names=NA)

# Density plot of raw counts across samples
png("Density.plot.gene.count.06132022.png", res=200, unit="in", height=8, width=11)
ggdensity(gene.dds.collapsed.count, x=colnames(gene.dds.collapsed.count), color=colnames(gene.dds.collapsed.count), add="median", merge=TRUE) +
    scale_x_log10() + xlab("Counts")
dev.off()

# Fit NGB model to gene counts from different factors of infection+origin (i.e. EBV:huNSG, EBV+KSHV:huNSG and et.c.)
gene.dds.collapsed=DESeq(gene.dds.collapsed, test="Wald")

resultsNames(gene.dds.collapsed)
# [1] "Intercept"                                   
# [2] "Origin.infection_EBV.KSHV.huNSG_vs_EBV.huNSG"
# [3] "Origin.infection_EBV.CL_vs_EBV.huNSG"        
# [4] "Origin.infection_KSHV.CL_vs_EBV.huNSG"       
# [5] "Origin.infection_EBV.KSHV.CL_vs_EBV.huNSG" 

# All groups are compared against reference=EBV.huNSG
# Pairwise comparison of groups using Wald significance test
out.wald.list=list(); ref.group="EBV.huNSG"
for (var in unique(gene.dds.collapsed$Origin.infection)){
    var=gsub("[+]",".",var)
    if (var != ref.group){
        res = results(gene.dds.collapsed, contrast=c("Origin.infection", var, ref.group),
independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
        res = lfcShrink(gene.dds.collapsed, coef=paste("Origin.infection_", var, "_vs_", ref.group, sep=""), type="apeglm", res=res)
        out.wald.list[[paste(var,"vs",ref.group,sep="")]]=res
    }
}

# DEG list from EBV+KSHV.huNSG vs EBV.huNSG test
genes.wald=data.frame(out.wald.list["EBV.KSHV.huNSGvsEBV.huNSG"][1]); sig.genes.wald=genes.wald[genes.wald$EBV.KSHV.huNSGvsEBV.huNSG.padj < 0.05 & genes.wald$EBV.KSHV.huNSGvsEBV.huNSG.log2FoldChange < -1 | genes.wald$EBV.KSHV.huNSGvsEBV.huNSG.log2FoldChange > 1, ]
nrow(sig.genes.wald)
# [1] 751

colnames(genes.wald)=gsub("EBV.KSHV.huNSGvsEBV.huNSG.", "", colnames(genes.wald))

# Save DEG list for EBV+KSHV.huNSG vs EBV.huNSG as text file
write.table(genes.wald, file=file.path(parent.directory, "DEG.List.EBV+KSHVvsEBVhuNSG.txt"), sep = "\t", col.names=NA)



################################## Analysis ########################################
#####################
##### Normalize the raw count data 
### Using VST transformation
gene.dds.collapsed.norm <- vst(gene.dds.collapsed, blind=TRUE)
write.table(assay(gene.dds.collapsed.norm), file.path(parent.directory, "Normalized.Counts.Cell.Lines.txt"), sep="\t", col.names=NA)

#### Compute PCA
res.pca <- prcomp(t(assay(gene.dds.collapsed.norm)), scale = FALSE)

### Scree plot for the variance explained by eigenvectors
png("Scree.plot.PCA.Cell.Lines.06142022.png", res=200, unit="in", height=8, width=11)
fviz_eig(res.pca)
dev.off()

### 3D PCA PLOT
variance.explained <- (res.pca$sdev)^2/(sum((res.pca$sdev)^2))*100

png("3DPCA.plot.Cell.Lines.06142022.png", res=200, unit="in", height=8, width=11)
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


### Heatmap: expression of gene-editing enzymes in samples
norm.gene.count <- assay(gene.dds.collapsed.norm)

###### APOBEC3A/B/C/D/E/F/G/H, AICDA, POLH, ADAR1, ADARB1, ADARB2
#### Note: APOBEC3E NOT IN EXPRESSION MATRIX
editing.genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D",  "APOBEC3F", "APOBEC3G", "APOBEC3H",
    "AICDA", "POLH", "ADAR", "ADARB1", "ADARB2")
temp.count <- norm.gene.count[rownames(norm.gene.count) %in% editing.genes, ]
temp.count = temp.count[editing.genes,]

#### Heatmap
colors <- colorRamp2(c(min(temp.count),max(temp.count)), c('blue', 'red')) ### set color scheme

# col_levels = list(Group = c("EBV-huNSG" = "green", "EBV+KSHV-huNSG" = "blue", "EBV-LCL" = "red", 
#     "EBV+KSHV-Cell.Line" = "purple", "KSHV-Cell.Line" = "yellow")) ### set color levels of group

col_levels = list(Group = c("EBV.huNSG" = "green", "EBV+KSHV.huNSG" = "blue", "EBV.CL" = "red", 
    "EBV+KSHV.CL" = "purple", "KSHV.CL" = "yellow")) ### set color levels of group

### Left annotation by group
ha_right <- rowAnnotation(Group=gene.dds.collapsed$Origin.infection, col = col_levels)

hp <- Heatmap(as.matrix(t(temp.count)), name = "Expression", col = colors, row_title = "Samples",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), row_split = gene.dds.collapsed$Origin.infection, 
    column_split=rep(c("A", "B", "C"), c(8,1,3)), right_annotation=ha_right, border = TRUE)

png("Heatmap.editing_enzyme.expression.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()

#### Load DEG list EBV+KSHV vs EBV
deg1 <- read.table(file.path(parent.directory, "DEG.List.EBV+KSHVvsEBVhuNSG.txt"), sep="\t", header=TRUE)

deg1 <- deg1 %>% mutate(Status = case_when(log2FoldChange >= 1 ~ 'UP', log2FoldChange <= -1 ~ 'DOWN', 
    log2FoldChange > -1 & log2FoldChange < 1 ~ 'NEUTRAL'))

deg1 <- deg1[deg1$X %in% editing.genes, ]

### Scatterplot
p <- ggplot(deg1, aes(x=log2FoldChange, y=-log10(padj), color=Status, na.rm=TRUE)) + geom_point() + 
        geom_text(aes(label = X), hjust=0.5, vjust=1.3) +
        theme_classic() +
        labs(x="log2FC", y="-log10(padj)") + 
        theme(axis.text.x = element_text(angle=90, hjust=1, size=10),  
        plot.title = element_text(size = 13, face = "bold"), axis.title=element_text(size=15),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14)) +
        scale_colour_manual(values = c("red", "grey", "blue")) +
        ggtitle("EBV+KSHV vs EBV huNSG mice DE editing genes") +
        geom_vline(xintercept=-1, colour="green", linetype = "longdash") +
        geom_vline(xintercept=1, colour="green", linetype = "longdash")

png("Scatterplot.DE.editing_enzyme.FC.png", res=200, unit="in", height=8, width=11)
p
dev.off()


###### expression of most commonly known genes in lymphomas
############# 
gene.list <- c("MYCBP", "BCL2", "BCL2A1", "BCL6", "IRF4", "FOXP1", "CXCR4", "CD27", "CD72", "PTPN6", "IFNGR1", "CAMK1", "CD22",
    "CD83", "NFKBIA", "LMO2", "BCL2A1", "MIR155HG", "CCR6", "BANK1", "RASGRP2", "PRDM1", "XBP1", "FKBP11", "BLNK",
        "IRAG2", "MYBL1", "RGS13", "TUBB2A", "UCHL1", "XIST", "STAG3", "TOX2", "CR2", "BATF", "P2RX5", "CCL22", "SH3BP5",
        "JADE3", "DOCK10", "NLRP7", "NIBAN3", "VPREB3", "FCMR", "GPR183", "FOXO1", "CD79A", "CD79B", "MS4A1", "CD19")

temp.count <- norm.gene.count[rownames(norm.gene.count) %in% gene.list, ]
temp.count = temp.count[gene.list,]

#### Heatmap
colors <- colorRamp2(c(min(temp.count),max(temp.count)), c('blue', 'red')) ### set color scheme

# col_levels = list(Group = c("EBV-huNSG" = "green", "EBV+KSHV-huNSG" = "blue", "EBV-LCL" = "red", 
#     "EBV+KSHV-Cell.Line" = "purple", "KSHV-Cell.Line" = "yellow")) ### set color levels of group

col_levels = list(Group = c("EBV.huNSG" = "green", "EBV+KSHV.huNSG" = "blue", "EBV.CL" = "red", 
    "EBV+KSHV.CL" = "purple", "KSHV.CL" = "yellow")) ### set color levels of group

### Left annotation by group
ha_right <- HeatmapAnnotation(Group=gene.dds.collapsed$Origin.infection, col = col_levels)

hp <- Heatmap(as.matrix(temp.count), name = "Expression", col = colors, column_title = "Samples",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), column_split = gene.dds.collapsed$Origin.infection, 
    top_annotation=ha_right, border = TRUE)

png("Heatmap.lymphoma_genes.expression.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()


temp.count <- temp.count[, which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

#### Generate rowwise z-scores
mu <- apply(temp.count, 1, mean); sdv <- apply(temp.count, 1, sd)

temp.count <- (temp.count-mu)/sdv

anno <- gene.dds.collapsed$Origin.infection[which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

#### Heatmap
colors <- colorRamp2(c(min(temp.count),max(temp.count)), c('blue', 'red')) ### set color scheme

# col_levels = list(Group = c("EBV-huNSG" = "green", "EBV+KSHV-huNSG" = "blue", "EBV-LCL" = "red", 
#     "EBV+KSHV-Cell.Line" = "purple", "KSHV-Cell.Line" = "yellow")) ### set color levels of group

col_levels = list(Group = c("EBV.huNSG" = "green", "EBV+KSHV.huNSG" = "blue")) ### set color levels of group

### Left annotation by group
ha_right <- HeatmapAnnotation(Group=anno, col = col_levels)

hp <- Heatmap(as.matrix(temp.count), name = "z-score", col = colors, column_title = "Samples",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), column_split = anno, 
    top_annotation=ha_right, border = TRUE)

png("Heatmap.lymphoma_genes.expression.huNSGonly.png", res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()


############## GENE SET ENRICHMENT ANALYSIS######################
gene.count.matrix <- cbind(NAME = rownames(norm.gene.count), norm.gene.count); rownames(gene.count.matrix) <- NULL

## Subset samples derived from EBV+ and EBV+KSHV+ huNSG mice
samples.keep <- c("NAME", "E4-3.R1", "E4-3.R2", "E4-3.R3", "E4-9.R1", "E4-9.R2", "E4-9.R3", "EK3-11.R1", "EK3-11.R2", "EK3-11.R3", "EK3-13.R1", "EK3-13.R2",
     "EK4-11.R1", "EK4-11.R2", "EK4-11.R3")
gene.count.matrix <- gene.count.matrix[, colnames(gene.count.matrix) %in% samples.keep]
all(colnames(gene.count.matrix)==samples.keep)
# [1] TRUE

# Save the count file
dir.count <- "/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/GSEA/EBV+KSHVvsEBVhuNSG/ExpressionMatrix"
if (!dir.exists(dir.count)){
    dir.create(dir.count, recursive = TRUE)
}

write.table(gene.count.matrix, file = file.path(dir.count, "Normalized_GSEA_EBV+KSHVvsEBV.gene_count_matrix.txt"), sep = "\t", row.names=FALSE)




#### Make phenotype labels file
## First line: number of samples, number of categories and 1
first_line <- paste(length(samples.keep)-1, "2", "1", sep=" ")

## Second line: '#', category 1, category 2
second_line <- paste("#", "EBV", "EBV+KSHV", sep=" ")

## Third line:
ebv_samples <- paste(rep("EBV", 6), sep=" "); ebv_kshv_samples <- paste(rep("EBV+KSHV", 8), sep = " ")
third_line <- paste(c(ebv_samples, ebv_kshv_samples), collapse=" ")

## Put into single file
lines <- c(first_line, second_line, third_line)
# Save the file
dir.pheno <- "/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/GSEA/EBV+KSHVvsEBVhuNSG/PhenotypeLabelsFile"
if (!dir.exists(dir.pheno)){
    dir.create(dir.pheno)
}

writeLines(lines, file.path(dir.pheno, "EBV+KSHVvsEBV.cls"))


######## pathway difference
############### ebv + kshv high
############## high mitochondrial activity, mitotic cell cycle/DNA replication, negative apoptosis,
############## immune response (antigen presentation), ER folding erad

##### Main Genes to plot: NDUFB6, NDUFC2, NDUFA4, NDUFV3, SDHC, SDHAF2, CYCS, ATP5F1B, ATP5PB, 
##### VCP, UFD1, NPLOC4, FAF2, UBE2J2, RNF5, HSPA5, SDF2L1, ANAPC15, ANAPC5, CDC20, CDK2, CDC25A, CDC25C, BAX, DNAJA3, CDKN1A, CDKN2A, BAD,
##### AIFM2, BCL2L12, RNF34, TRIAP1, TRAIP, PRELID1, FYB1, TRGV5, HLA-B, HLA-DRB5, HLA-DRB1, HLA-C, CD7, APOBEC3A, APOBEC3B, IGHLC3, IGLC7
##### mitochondrial electron transport, NADH to ubiquinone: 1-9, ERAD pathway: 10-17, 
##### anaphase-promoting complex-dependent catabolic process: 18-20, G2/M transition of mitotic cell cycle: 21-23, 
##### negative regulators of apoptosis: 24-28, positive regulators of apoptosis: 29-34, immune response: 35-41, IMMUNOGLOBULIN: 2



##### Supplementary Genes to plot: ASF1A, ASF1B, CHAF1A, CHAF1B, GINS4, GINS3, RAD51, SSBP1, RECQL4, MCM4, MCM7,
##### MPV17L2, MRPL51, MRPL43, MRPS31, MRPS11, TUFM, LARS2,  ECHS1, SLC27A4, CPT2, CROT, ACOX1, HRAS, 
##### DNA replication-dependent chromatin assembly: 1-4, DNA unwinding involved in DNA replication: 5-11, 
##### organelle assembly: 12-17, mitochondrial translation: 18-19, fatty acid catabolism: 19-23, OTHER: 24


##### IMMUNOGLOBULINS: IGLV5-37, IGHV4-61, IGLV1-40, IGKV1D-37, IGHV6-1, IGHV1-46, IGHV4-34, IGLV3-25, IGKV1-37,
##### IGLC7, IGLV1-50, IGLV1-44, IGLV3-16, IGLC3


############## ebv high
############## type 1 IFN antiviral response, negative regulation of antigen receptor-mediated signaling pathway, 
############## positive regulation of cytokine production involved in immune response, MAPK cascade,
############## positive regulation of stress-activated MAPK cascade

###### Main Genes to plot: OAS1, IFIT1, IRF9, IFITM3, SP100, TBKBP1, EIF2AK2, PVRIG, FCGR2B, CBLB, CD300A, CD22, DUSP22,
###### TEC, SLAMF1, XCL1, MIF, TRIM6, TNFRSF14, MAPK12, MAP3K8, MAP3K2, MAPKAPK3, MEF2A, GADD45B, TAOK1, TAOK3, RELL2
###### type 1 IFN antiviral response: 1-7, negative regulation of antigen receptor-mediated signaling pathway: 8-13,
###### positive regulation of cytokine production involved in immune response: 14-19, MAPK cascade: 20-24
###### positive regulation of stress-activated MAPK cascade: 25-28



###### Supplementary Genes to plot: SHFL, RESF1, DHX58, DDX58, AICDA, UNC93B1, TLR10, IRAK3, CD40, RIPK2, NOD1, PTEN, SMAD3, SMAD4,
###### KRAS
###### DEFENSE RESPONSE TO VIRUS:1:4, pattern recognition receptor signaling pathway:5:10, OTHER: 11:14

###### IMMUNOGLOBULINS: IGKV1D-33, IGKV2D-28, IGLJ1, IGLV6-57, IGHV3-69-1

### Load GSEA results
dir.gsea <- "/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/GSEA"
ebv.kshv.gsea <- read.table(file.path(dir.gsea, "EBV.KSHVhuNSG.Gene.Set.Enrichment.txt"), header=TRUE, sep="\t", quote="")
ebv.gsea <- read.table(file.path(dir.gsea, "EBVhuNSG.Gene.Set.Enrichment.txt"), header=TRUE, sep="\t", quote="")

column.names <- c("GO.biological.process", "Homo sapiens - REFLIST (20589)", "detected", "expected", "over/under", "fold.Enrichment",
    "raw.P-value", "FDR")
colnames(ebv.kshv.gsea) <- column.names; colnames(ebv.gsea) <- column.names

ebv.kshv.gsea <- ebv.kshv.gsea[-which(ebv.kshv.gsea["fold.Enrichment"] == " < 0.01"),]
ebv.gsea <- ebv.gsea[-which(ebv.gsea["fold.Enrichment"] == " < 0.01"),]

ebv.gsea["fold.Enrichment"] <- as.numeric(unlist(ebv.gsea["fold.Enrichment"]))*-1

go.ebv.kshv.select <- c("mitochondrial respiratory chain complex assembly (GO:0033108)", 
    "proton motive force-driven ATP synthesis (GO:0015986)", "mitochondrial translation (GO:0032543)",
    "ERAD pathway (GO:0036503)", "endoplasmic reticulum unfolded protein response (GO:0030968)", 
    "mitotic cell cycle (GO:0000278)", "DNA replication (GO:0006260)",
    "immune response (GO:0006955)", "organelle assembly (GO:0070925)", "fatty acid catabolic process (GO:0009062)", 
    "apoptotic process (GO:0006915)")
go.ebv.select <- c("innate immune response (GO:0045087)", "response to type I interferon (GO:0034340)", "negative regulation of antigen receptor-mediated signaling pathway (GO:0050858)",
    "pattern recognition receptor signaling pathway (GO:0002221)", "positive regulation of cytokine production involved in immune response (GO:0002720)",
        "positive regulation of stress-activated MAPK cascade (GO:0032874)", "MAPK cascade (GO:0000165)", 
        "negative regulation of viral genome replication (GO:0045071)")


ebv.kshv.gsea <- ebv.kshv.gsea[ebv.kshv.gsea[,"GO.biological.process"] %in% go.ebv.kshv.select, ]

ebv.gsea <- ebv.gsea[ebv.gsea[,"GO.biological.process"] %in% go.ebv.select, ]

gsea.df <- rbind(ebv.kshv.gsea, ebv.gsea)

gsea.df["Group"] <- c(rep("Up", nrow(ebv.kshv.gsea)), rep("Down", nrow(ebv.gsea)))

gsea.df <- gsea.df %>% mutate(log10FDR =
    case_when(Group == "Up" ~ -log10(FDR),
                Group == "Down" ~ log10(FDR)))

gsea.df <- gsea.df %>% arrange(factor(GO.biological.process, levels = c(go.ebv.kshv.select, go.ebv.select)))

#### Create Figures folder
dir.gsea.figures <- "/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/GSEA/Figures"
if (!dir.exists(dir.gsea.figures)){
    dir.create(dir.gsea.figures, recursive = TRUE)
}

##### Figure with enrichment pathway
p <- ggplot(gsea.df, aes(x=GO.biological.process,, y=log10FDR, fill=Group)) + 
    geom_bar(stat='identity') +
    geom_text(aes(label=fold.Enrichment), hjust = ifelse(gsea.df$log10FDR>0, 1.1, -0.1), color = "white") +
    xlab("Biological process") + ylab("log10(adj p-value)") +
    scale_fill_manual(values=c("Up"="red", "Down"="blue")) +
    scale_x_discrete(limits=rev(gsea.df$GO.biological.process)) +
    theme_classic() +
    ggtitle("EBV+KSHV vs EBV huNSG mice GO analysis") 

png(file.path(dir.gsea.figures,"Enrichment.pathways.EBV+KSHVvsEBV.huNSG.png"), res=200, unit="in", height=8, width=11)
p + coord_flip()
dev.off()

##### Main heatmap: 43+28 genes
main.genes <- c('NDUFB6', 'NDUFC2', 'NDUFA4', 'NDUFV3', 'SDHC', 'SDHAF2', 'CYCS', 'ATP5F1B', 'ATP5PB', 
    'VCP', 'UFD1', 'NPLOC4', 'FAF2', 'UBE2J2', 'RNF5', 'HSPA5', 'SDF2L1', 'ANAPC15', 'ANAPC5', 'CDC20', 'CDK2', 'CDC25A', 
    'CDC25C', 'BAX', 'DNAJA3', 'CDKN1A', 'CDKN2A', 'BAD', 'AIFM2', 'BCL2L12', 'RNF34', 'TRIAP1', 'TRAIP', 'PRELID1', 'FYB1', 
    'TRGV5', 'HLA-B', 'HLA-DRB5', 'HLA-DRB1', 'HLA-C', 'CD7', 'APOBEC3A', 'APOBEC3B', 'IGLC3', 'IGLC7', 'OAS1', 'IFIT1', 'IRF9', 'IFITM3', 'SP100', 
    'TBKBP1', 'EIF2AK2', 'PVRIG', 'FCGR2B', 'CBLB', 'CD300A', 'CD22', 'DUSP22', 'TEC', 'SLAMF1', 'XCL1', 'MIF', 'TRIM6', 
    'TNFRSF14', 'MAPK12', 'MAP3K8', 'MAP3K2', 'MAPKAPK3', 'MEF2A', 'GADD45B', 'TAOK1', 'TAOK3', 'RELL2')

temp.count <- data.frame(norm.gene.count[rownames(norm.gene.count) %in% main.genes, ])
temp.count <- temp.count %>% arrange(factor(rownames(temp.count), levels = main.genes))

temp.count <- temp.count[, which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

#### Generate rowwise z-scores
mu <- apply(temp.count, 1, mean); sdv <- apply(temp.count, 1, sd)

temp.count <- (temp.count-mu)/sdv

##### Heatmap
colors <- colorRamp2(c(min(temp.count),max(temp.count)), c('blue', 'red')) ### set color scheme

col_levels = list(Group = c("EBV.huNSG" = "green", "EBV+KSHV.huNSG" = "blue")) ### set color levels of group

anno <- gene.dds.collapsed$Origin.infection[which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

### Left annotation by group
ha_right <- HeatmapAnnotation(Group=anno, col = col_levels)

hp <- Heatmap(as.matrix(temp.count), name = "z-score", col = colors, column_title = "Samples",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), column_split = anno, 
    top_annotation=ha_right, border = TRUE)

png(file.path(dir.gsea.figures, "Heatmap.DE.genes.GO.analysis.expression.png"), res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()
##### Supplementary heatmap: 20+14
sup.genes <- c('ASF1A', 'ASF1B', 'CHAF1A', 'CHAF1B', 'GINS4', 'GINS3', 'RAD51', 'SSBP1', 'RECQL4', 'MCM4', 'MCM7',
    'MPV17L2', 'MRPL51', 'MRPL43', 'MRPS31', 'MRPS11', 'TUFM', 'LARS2',  'ECHS1', 'SLC27A4', 'CPT2', 'CROT', 'ACOX1', 'HRAS',
    'SHFL', 'RESF1', 'DHX58', 'DDX58', 'AICDA', 'UNC93B1', 'TLR10', 'IRAK3', 'CD40', 'RIPK2', 'NOD1', 'PTEN', 'SMAD3', 'SMAD4', 'KRAS',
    'MYC')

temp.count <- data.frame(norm.gene.count[rownames(norm.gene.count) %in% sup.genes, ])
temp.count <- temp.count %>% arrange(factor(rownames(temp.count), levels = sup.genes))

temp.count <- temp.count[, which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

#### Generate rowwise z-scores
mu <- apply(temp.count, 1, mean); sdv <- apply(temp.count, 1, sd)

temp.count <- (temp.count-mu)/sdv

##### Heatmap
colors <- colorRamp2(c(min(temp.count),max(temp.count)), c('blue', 'red')) ### set color scheme

col_levels = list(Group = c("EBV.huNSG" = "green", "EBV+KSHV.huNSG" = "blue")) ### set color levels of group

anno <- gene.dds.collapsed$Origin.infection[which(gene.dds.collapsed$Origin.infection == "EBV.huNSG" |
    gene.dds.collapsed$Origin.infection == "EBV+KSHV.huNSG")]

### Left annotation by group
ha_right <- HeatmapAnnotation(Group=anno, col = col_levels)

hp <- Heatmap(as.matrix(temp.count), name = "z-score", col = colors, column_title = "Samples",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 8), 
    column_names_gp = gpar(fontsize = 10), column_split = anno, 
    top_annotation=ha_right, border = TRUE)

png(file.path(dir.gsea.figures, "Heatmap.DE.genes.GO.analysis.supplementary.expression.png"), res=200, unit="in", height=8, width=11)
draw(hp)
dev.off()

##### GOPlot
deg1 <- read.table(file.path(parent.directory, "DEG.List.EBV+KSHVvsEBVhuNSG.txt"), sep="\t", header=TRUE)
colnames(deg1) <- c("ID", "AveExpr", "logFC", "lfcSE", "P.Value", "adj.P.Val")

go.output <- read.table(file.path(dir.gsea, "go.output.ebv+kshvvsebv.hunsg.txt"), sep="\t", header=TRUE)

circ <- circle_dat(go.output, deg1)
genes <- deg1[deg1$ID %in% c(main.genes, sup.genes), c("ID", "logFC")]
chord <- chord_dat(data = circ, genes = genes)

png(file.path(dir.gsea.figures, "GOChord.huNSGmice.png"), res=200, unit="in", height=20, width=27)
GOChord(chord, title = "Circos plot (EBV+KSHV vs EBV GO analysis)", space = 0.02, gene.order = 'logFC', gene.space = 0.23, gene.size = 7, process.label = 15, lfc.min = -20,
    lfc.max = 20, lfc.col = c("red", "white", "blue"))
dev.off()




















