#!/bin/bash

################ This script generates count matrix tables on gene and transcript level using Stringtie
################ It consists of 4 steps: 1) Generate sample GTF file using reference genome GTF, 2) Merge GTFs, 3) Generate sample GTF using merged GTF
################ 4) Make count matrix tables on gene and transcript level
################ Input: BAM files
################ Output: Count matrix tables

### Set required variables: 
export NTHREADS=16
#### Directory with BAM files
export dir_bam=/home/cluster/abalay/scratch/Cell_Lines
#### Set dir for count matrix
export output_dir=/home/cluster/abalay/scratch/Cell_Lines/Transcriptomics/Count_matrices/Stringtie

### Please do not change these variables unless they do not exist
export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files
export stringtie=$script_dir/Packages/stringtie-2.1.4/stringtie 
export ref_gen=$script_dir/Files_for_scripts/GRCh38.p13.genome.fa
export gtf_file=$script_dir/Files_for_scripts/gencode.v37.annotation.gtf
export gff_file=$script_dir/Files_for_scripts/gencode.v37.annotation.gff3

############### Step 1: Generate sample GTF file using reference genome GTF
############### Input: BAM file
############### Output: GTF file

## Create directory for count matrices
if [ ! -d $output_dir ];then
   mkdir -p $output_dir
fi

## Create preliminary text file with list to path to GTF file
touch $output_dir/gtf_list.txt

## Generate sample GTF file 
echo "Generating GFT file"

for file in $(find $dir_bam -name "*.bam")
do
    curr=$(basename "$file" .bam)
    $stringtie $file -G $gtf_file -eB -v -p $NTHREADS -o $output_dir/"$curr".gtf 
    echo $output_dir/"$curr".gtf >> $output_dir/gtf_list.txt
done


################ Step 2: Merge GTFs
################ Input: GTF file
################ Output: Merged GTF file with non-redundant transcripts

## Merge GTFs
echo "Generate merging GTF"
for file in $(find $output_dir -name "*.gtf")
do
   curr=$(basename "$file" .gtf)
   $stringtie --merge -G $gtf_file -p $NTHREADS -o $output_dir/merged.gtf  $output_dir/gtf_list.txt
done

################# Step 3: Generate sample GTF file using merged GTF file
################# Input: BAM file
################# Output: GTF file
## Touch file with sample IDs and path to every sample GTF file
touch $output_dir/transcript_human_sample_list.txt

## Generate GTF file using merged GTF file and BAM file from Step 1 in order to produce read coverage tables per transcript
for file in $(find $dir_bam -name "*.bam")
do
    curr=$(basename "$file" .bam)
    printf "$curr\t$output_dir/$curr.gtf\n" >> $output_dir/transcript_human_sample_list.txt # generate a file with sample IDs and path to every sample GTF file produced by st$
    rm $output_dir/"$curr".gtf
    $stringtie $file -G $output_dir/merged.gtf -eB -p $NTHREADS  -o $output_dir/"$curr".gtf
done

################## Step 4: Generate count matrix tables on gene and transcript levels
################## Input: File with sample ID and path to GTF file
################## Output: Gene and transcript count matrices
export stringtie=$script_dir/Packages/stringtie-2.1.4
$stringtie/prepDE.py -i $output_dir/transcript_human_sample_list.txt  -g $output_dir/Cell_Lines_gene_count.csv  -t $output_dir/Cell_Lines_transcript_count.csv

## Remove unused files
rm $output_dir/transcript_human_sample_list.txt $output_dir/gtf_list.txt $output_dir/*.gtf $output_dir/merged.gtf





