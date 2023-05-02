export WD=/home/ubuntu
export NTHREADS=32
export rseqc=$WD/RSeQC-3.0.1/scripts
export Illumina_truseq=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
export multiqc=$WD/MultiQC/multiqc
export server_dir=/mnt/smb_share/RCLab-share-UZH/Agshin\ Data\ analysis/AA_3DPDO_RNAseq/2021\ RAW\ RNASeq
export featureCounts=$WD/subread-2.0.1-Linux-x86_64/bin/featureCounts
export seqtk=$WD/seqtk/seqtk

cd "$server_dir"

# QC
#for file in $(find . -name "*.fq.gz")
#do
#   fastqc $file
#done

#cp -r $(find . -name "*_fastqc.html") .

# Adapter trimming

#for file in $(find ./LncRNA -name "*_1.fq.gz")
#do
#  id=$(echo $file | awk -F '_1.fq.gz' '{print $1}')
#  cutadapt -o "$id"_cutadapt_1.fastq.gz -p "$id"_cutadapt_2.fastq.gz -q 10 -m 35 -j $NTHREADS "$id"_1.fq.gz "$id"_2.fq.gz > "$id"_cutadapt_report.txt
#done

#for file in $(find ./miRNA -name "*.fq.gz")
#do
#  id=$(echo $file | awk -F '.fq.gz' '{print $1}')
#  cutadapt -o "$id"_cutadapt_1.fastq.gz -q 10 -j $NTHREADS "$id".fq.gz  > "$id"_cutadapt_report.txt
#done

# Alignment

#for file in $(find ./LncRNA -name "*cutadapt_1.fastq.gz")
#do
#  id=$(echo $file | awk -F '_cutadapt_1.fastq.gz' '{print $1}')
#  gsnap --gunzip -d human_genome.fa.gz -D $WD/gmap-2020-06-01  -t $NTHREADS  -n 10 -N 1 -s $WD/human_splicesite_file.iit -g $WD/human_annotation.iit  \
#  -A sam "$id"_cutadapt_1.fastq.gz "$id"_cutadapt_2.fastq.gz \
#  | samtools view -b | samtools sort -@ 32 > "$id"_cutadapt_sorted.bam
#  samtools index "$id"_cutadapt_sorted.bam
#done

#for file in $(find ./miRNA/Cells/Batch2/C_915 -name "*cutadapt_1.fastq.gz")
#do
#  id=$(echo $file | awk -F '_cutadapt_1.fastq.gz' '{print $1}')
#  echo $id
#  gsnap --gunzip -d mirna_db.fa.gz -D $WD/gmap-2020-06-01 -t $NTHREADS  -n 10 \
#  -A sam "$id"_cutadapt_1.fastq.gz \
#  | samtools view -b > "$id"_mirna_cutadapt_sorted.bam
#  samtools view -f 4 "$id"_mirna_cutadapt_sorted.bam | samtools sort -@ $NTHREADS | samtools fastq > "$id"_mirna_cutadapt_sorted.fastq.gz 
#  samtools view -b -F 4 "$id"_mirna_cutadapt_sorted.bam > "$id"_mirna_aln_cutadapt_sorted.bam 
#  gsnap --gunzip -d human_genome.fa.gz -D $WD/gmap-2020-06-01  -t $NTHREADS  -n 10 -N 1 -s $WD/human_splicesite_file.iit -g $WD/human_annotation.iit  \
#  -A sam "$id"_mirna_cutadapt_sorted.fastq.gz \
#  | samtools view -b | samtools sort -@ 32 > "$id"_mirna_hcutadapt_sorted.bam
#  samtools view -b -F 4 "$id"_mirna_hcutadapt_sorted.bam > "$id"_mirna_aln_hcutadapt_sorted.bam
#  samtools merge -@ $NTHREADS "$id"_mirna_aln_sorted.bam "$id"_mirna_aln_cutadapt_sorted.bam "$id"_mirna_aln_hcutadapt_sorted.bam
#done

# Mapping QC

#touch ./infer_experiment.txt
#touch ./rrna_contamination.txt
#for file in $(find ./miRNA -name "*_cutadapt_sorted.bam")
#do
#  curr=$(basename "$file" _cutadapt_sorted.bam)
#  echo "Strandedness for $file" >> infer_experiment.txt
#  python3 $rseqc/infer_experiment.py -r $WD/Homo_sapiens.GRCh38.79.bed  -i $file >> infer_experiment.txt # Identify strandedness of reads
#  python3 $rseqc/read_duplication.py -i $file  -o "$curr"
#  samtools view -@ $NTHREADS -hb -U "$curr".rrna_free.bam -L $WD/hg38_rRNA.bed  $file > "$curr".in.bam
#  Total_reads=$(samtools view -@ $NTHREADS -hb -c $file)
#  rRNA_hits=$(samtools view -@ $NTHREADS -hb -c "$curr".in.bam)
#  printf ""$curr"_cutadapt_sorted\t$Total_reads\t$rRNA_hits\n" >> rrna_contamination.txt
#  rm "$curr".in.bam 
#done

#mkdir -p "$server_dir"/LncRNA/Counts/lncRNA
#cd "$server_dir"/LncRNA/Counts/lncRNA

## For LncRNA
#$featureCounts -p -s 2 -T $NTHREADS -t 'exon' -g 'gene_id'  -a $WD/gencode.v34.long_noncoding_RNAs.gtf.gz -o lncRNA_counts.txt  $server_dir/LncRNA/*_cutadapt_sorted.bam

#mkdir -p "$server_dir"/LncRNA/Counts/mRNA
#cd "$server_dir"/LncRNA/Counts/mRNA

## For mRNA
### by gene_id
#$featureCounts -p -s 2 -T $NTHREADS -t 'exon' -g 'gene_id'  -a $WD/gencode.v37.annotation.gtf -o All_RNA_counts_by_geneid.txt  "$server_dir"/LncRNA/*_cutadapt_sorted.bam

### by gene_name
#$featureCounts -p -s 2 -T $NTHREADS -t 'exon' -g 'gene_name'  -a $WD/gencode.v37.annotation.gtf -o All_RNA_counts_by_genename.txt  "$server_dir"/LncRNA/*_cutadapt_sorted.bam
#cd "$server_dir"/miRNA

#for file in $(find ./miRNA/EVs/Batch2/E_915 -name "*_mirna_aln_sorted.bam")#do
#  id=$(echo $file | awk -F '_mirna_aln_sorted.bam' '{print $1}')
#  samtools view -h "$id"_mirna_aln_sorted.bam | sed 's/\tXA\:Z\:[^\t]*//' | samtools view -Sb -h > "$id"_mirna_aln_sorted_reduced.bam
#$featureCounts -s 0 -T $NTHREADS -t 'miRNA' -g 'Name' -a $WD/hsa.gff3 -o mirna_counts.txt  *_mirna_aln_sorted.bam
#done

cd "$server_dir"
mkdir -p "$server_dir"/miRNA/miRNA_fasta_files
for file in $(find ./miRNA -maxdepth 1 -name "*_mirna_aln_sorted.bam")
do
   curr=$(basename "$file" _mirna_aln_sorted.bam)
   samtools bam2fq "$file" | $seqtk seq -A > ./miRNA/miRNA_fasta_files/"$curr"_mirna_aln_sorted.fa
done


