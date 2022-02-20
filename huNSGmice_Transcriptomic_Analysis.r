################

library("ggpubr")
library("DESeq2")
library("readxl")

## Load gene matrix of infected huNSG mice samples
parent.directory="/mnt/smb_share/AGS_AB/Cell_Lines/Transcriptomic_Profiles/Count.Matrices/"
gene_count=read.csv(paste(parent.directory,'huNSGTissue_gene_count.csv',sep=""), sep=",", header=T, row.names=1); colnames(gene_count)=gsub("[_].*", "", colnames(gene_count))
nrow(gene_count)
#  [1] 93964

## Load annotation table 
reference.dir="/mnt/smb_share/AGS_AB/Pipelines_and_Files/Files_for_scripts"
anno.table=data.frame(read_excel(paste(reference.dir, "/EBV-transformed_lymphoblast_cell_lines.xlsx",sep=""), sheet="EBV KSHV in vivo samples"))
colnames(anno.table)=anno.table[1,]; anno.table=anno.table[-1,]; anno.table=anno.table[1:(dim(anno.table)[1]-8),]; rownames(anno.table)=anno.table[,1]; anno.table=anno.table[,-1]
anno.table$Origin.infection=paste(anno.table$`Infection agent`, anno.table$`Origin`, sep=".") # create artificial column with infection agent + sample cellular origin
anno.table$Cell.replicate=paste(anno.table$`Cell line`, anno.table$`Replicate #`, sep=".") # create artificial column with cell line + replicate no
colnames(anno.table) = gsub(" ",".",colnames(anno.table)); anno.table$Origin.infection=factor(anno.table$Origin.infection); anno.table$Cell.replicate=factor(anno.table$Cell.replicate)

# Make sure each row in annotation table correspond to each column in gene matrix
all(rownames(anno.table) == colnames(gene_count))
# [1] TRUE

## Remove genes with zero counts across samples
gene_count=gene_count[rowSums(gene_count[,])!=0,] # remove genes with zero counts across all samples
nrow(gene_count)
# [1] 76569

## Check for duplicate genes in gene matrix
gene_count$gene.symbols = gsub(".*[|]", "", rownames(gene_count))
duplicate.gene.count = gene_count[duplicated(gene_count$gene.symbols),]
table(duplicate.gene.count$gene.symbols)
#  5_8S_rRNA     5S_rRNA         7SK  AC008686.1  AC010618.3  AC068587.4 
#           3           6           6           1           1           1 
#  AC091132.4  AC099520.1  AC104041.1  AC106795.1  AC145350.2  AJ271736.1 
#           1           1           1           1           1           1 
#     AKAP17A  AL110114.1  AL139317.3  AL157388.1  AL163952.1  AL390957.1 
#           1           1           1           1           1           1 
#  AL672277.1  AL683807.1  AL683807.2  AL732314.4  AL732314.6      AMD1P2 
#           1           1           1           1           1           1 
#  AP001107.9        ASMT       ASMTL   ASMTL-AS1      ATP11A        CAST 
#           1           1           1           1           1           1 
#      CCDC13        CD99      CD99P1        CDH3   CLCA4-AS1     CNTNAP2 
#           1           1           1           1           1           1 
#       CRLF2       CRLF3      CSF2RA    DDX11L16       DHRSX   DHRSX-IT1 
#           1           1           1           1           1           1 
#       ENPP2       EXOC4     F11-AS1     FAM182A   GABARAPL3       GOLM2 
#           1           1           1           1           1           1 
#        GPX3        GRK4      GTPBP6    H2AZ1-DT     HERC2P7      IL1RAP 
#           1           1           1           1           1           1 
#       IL3RA        IL9R KBTBD11-OT1   LINC00102   LINC00106   LINC00299 
#           1           1           1           1           1           1 
#   LINC00484   LINC00486   LINC00685   LINC01088   LINC01238   LINC02203 
#           1           1           1           1           1           1 
#   LINC02245   LINC02256   MEF2C-AS2 Metazoa_SRP       MGST3     MIR6089 
#           1           1           1         119           1           1 
#       MPRIP        MSRA       NDST1         OAF       P2RY8        PBX1 
#           1           1           1           1           1           1 
#       PDE3B       PIAS2     PIP4K2A      PLCXD1        PLD1       POLA1 
#           1           1           1           1           1           2 
#     PPP2R3B        REST        RGS6        ROM1     RPL14P5       SCN2A 
#           1           1           1           1           1           1 
#      SCNN1B       SGMS2     SLC19A1     SLC25A6     SNORA62     SNORA63 
#           1           1           1           1           3           1 
#     SNORA70     SNORA72     SNORA73     SNORA74     SNORA75     SNORD27 
#          12           2           2           1           3           1 
#     SNORD63   SOX21-AS1   STX17-AS1       SYNE2     TMSB15B        TNKS 
#           1           1           1           1           1           1 
#        TP53      TRPC6P          U1          U2          U3          U4 
#           1           1           2          15          11           4 
#          U6          U7          U8       VAMP7       Vault      WASH6P 
#          11           3           5           1           2           1 
#      WASIR1       Y_RNA       ZBED1     ZSCAN5A 
#           1         369           1           1 

## Delete rRNA genes + genes with multple loci + small nucleolar RNA genes
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
    print(expression)
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
# [1] 75538

# Put infection by cell origin as main source of variation
gene_dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = anno.table, design = ~ Origin.infection)

# Remove all genes with total counts < 10 across samples
genes.to.keep = rowSums(counts(gene_dds)) > 10; gene_dds <- gene_dds[genes.to.keep,]
dim(gene_dds)
# [1] 66218    79

gene_dds$Origin.infection = factor(gene_dds$Origin.infection, levels=unique(gene_dds$Origin.infection))

# Collapse technical replicates to one biological replicate per sample
gene.dds.collapsed = collapseReplicates(gene_dds, anno.table$Cell.replicate)

# Check if sum of counts for sample before collapse == counts for sample post collapse
matchFirstLevel <-  anno.table$Cell.replicate == levels(anno.table$Cell.replicate)[1]
all(rowSums(counts(gene_dds[,matchFirstLevel])) == counts(gene.dds.collapsed[,1]))
# [1] TRUE

# List of samples
colnames(gene.dds.collapsed)
#  [1] "AP2.R1"    "AP3.R1"    "AP5.R1"    "BCBL1.R1"  "E4-3.R1"   "E4-3.R2"  
#  [7] "E4-3.R3"   "E4-9.R1"   "E4-9.R2"   "E4-9.R3"   "EK3-11.R1" "EK3-11.R2"
# [13] "EK3-11.R3" "EK3-13.R1" "EK3-13.R2" "EK4-11.R1" "EK4-11.R2" "EK4-11.R3"
# [19] "HBL6.R1"   "LCL1.R1"   "LCL1.R2"   "LCL1.R3"   "LCL2.R1"   "LCL2.R2"  
# [25] "LCL2.R3"   "LCL3.R1"   "LCL3.R2"   "LCL3.R3" 

gene.dds.collapsed.count <- data.frame(counts(gene.dds.collapsed))

# Density plot of raw counts across samples
png("Density.plot.gene.count.02192022.png", res=200, unit="in", height=8, width=11)
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
# [1] 1658

colnames(genes.wald)=gsub("EBV.KSHV.huNSGvsEBV.huNSG.", "", colnames(genes.wald))

# Save DEG list for EBV+KSHV.huNSG vs EBV.huNSG as text file
write.table(genes.wald, file="DEG.List.EBV+KSHVvsEBVhuNSG.txt", sep = "\t", col.names=NA)
system(paste("sudo cp -r DEG.List.EBV+KSHVvsEBVhuNSG.txt ", parent.directory, sep=""))



