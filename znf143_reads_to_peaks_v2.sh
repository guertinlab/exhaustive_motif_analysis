#!/bin/bash
#SBATCH --job-name=znf143_reads_to_peaks
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load bowtie2
module load py-macs3
module load ucsc_genome
module load sratoolkit
module load bedtools
module load cutadapt
module load samtools
module load meme/5.5.9

#Set up
home=/scratch/afrutos/GSE266489_control3
ncore=32
cd $home

#Get the latest genome
mkdir $home/hs1
cd $home/hs1
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz -P $home/hs1

#Build indices
gunzip hs1.fa.gz
bowtie2-build --threads $ncore hs1.fa hs1

#Get Markov model for background
fasta-get-markov -m 3 hs1.fa > hs1_bkgrnd.txt

#Get blacklist
wget https://hgdownload.cse.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
gunzip hg38ToHs1.over.chain.gz
awk '{OFS="\t";} {print $1, $2, $3}' hg38-blacklist.v2.bed > hg38-blacklist.v2.test.bed
liftOver hg38-blacklist.v2.test.bed hg38ToHs1.over.chain hs1-blacklist.v2.bed hs1-blacklist.v2.unmapped.bed

#Make directories
cd $home
mkdir -p Motif_analysis_automation
cd Motif_analysis_automation
mkdir -p SRA FASTQ BAM BED Peaks Peaks/temp_macs MEME TOMTOM

#Get the SRAs
cd SRA

for i in 68 82 79 78
do
  echo $i
  prefetch SRR288794${i}
done

cd ..

#Make FASTQs
cd FASTQ

for i in 68 82 79 78
do
  echo $i
  fasterq-dump ../SRA/SRR288794${i}/SRR288794${i}.sra
done

#Rename files
mv SRR28879468_1.fastq HEK_CloneZD29_HA_control_rep1_PE1.fastq
mv SRR28879468_2.fastq HEK_CloneZD29_HA_control_rep1_PE2.fastq
mv SRR28879482_1.fastq HEK_CloneZD29_HA_control_rep2_PE1.fastq
mv SRR28879482_2.fastq HEK_CloneZD29_HA_control_rep2_PE2.fastq

mv SRR28879479_1.fastq HEK_CloneZD29_mock_control_rep1_PE1.fastq
mv SRR28879479_2.fastq HEK_CloneZD29_mock_control_rep1_PE2.fastq
mv SRR28879478_1.fastq HEK_CloneZD29_mock_control_rep2_PE1.fastq
mv SRR28879478_2.fastq HEK_CloneZD29_mock_control_rep2_PE2.fastq

cd ..

#Process reads
cd $home/Motif_analysis_automation/FASTQ
for file in *_PE1.fastq
do
    #Extract name (everything before _PE1.fastq)
    name="${file%%_PE1.fastq}"
    echo "Processing sample: $name"
    
    #Cut off adapters (I have to use a specific mamba environment for cutadapt)
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
             --cores=$ncore -m 10 -O 1 \
             -o ${name}_PE1_no_adapt.fastq \
             -p ${name}_PE2_no_adapt.fastq \
             ${name}_PE1.fastq ${name}_PE2.fastq \
             2>&1 | tee ${name}_cutadapt.log

    #Align and sort by name
    bowtie2 -p $ncore \
            --maxins 800 \
            -x $home/hs1/hs1 \
            -1 ${name}_PE1_no_adapt.fastq \
            -2 ${name}_PE2_no_adapt.fastq \
            | samtools sort -@ $ncore -n -o ../BAM/${name}.bw.bam

    #Fix mates, coordinate sort, and mark duplicates
    samtools fixmate -m ../BAM/${name}.bw.bam - \
        | samtools sort -@ $ncore - \
        | samtools markdup -s -r - ../BAM/${name}.hs1.bam
done

cd ..

#Call peaks
cd Peaks
#(I have to use a specific mamba environment for macs3)
macs3 callpeak --call-summits -t ../BAM/*HA_control*hs1.bam -c ../BAM/*mock*hs1.bam -n ZNF143_ChIP -g hs -q 0.01 --keep-dup all -f BAMPE --nomodel --tempdir temp_macs

#Call mock peaks
macs3 callpeak --call-summits -t ../BAM/*mock*hs1.bam -n ZNF143_mock -g hs -q 0.01 --keep-dup all -f BAMPE --nomodel --tempdir temp_macs

#Initialize variables
blacklist=$home/hs1/hs1-blacklist.v2.bed
sizes=$home/hs1/hs1.chrom.sizes.txt
genome=$home/hs1/hs1.fa
dir=$home/Motif_analysis_automation/Peaks
cd $dir

#Remove junk and get a window centered on the summit
name=ZNF143_ChIP
grep -v "random" ${name}_summits.bed | grep -v "chrUn" | grep -v "chrEBV" | grep -v "chrM" | grep -v "alt" | intersectBed -v -a - -b $blacklist > ${name}_summits_final.bed
slopBed -b 50 -i ${name}_summits_final.bed -g $sizes  | sort -k1,1 -k2,2n > ${name}_summit_100window.bed
slopBed -b 250 -i ${name}_summits_final.bed -g $sizes  | sort -k1,1 -k2,2n > ${name}_summit_500window.bed
fastaFromBed -fi $genome -bed ${name}_summit_100window.bed -fo ${name}_summit_100window.fasta

#Sort out the top 1000 peaks
sort -nrk5,5 ${name}_summit_100window.bed | head -n 1000 > ${name}_top1000_summit_100window.bed
fastaFromBed -fi $genome -bed ${name}_top1000_summit_100window.bed -fo ${name}_top1000_summit_100window.fasta

#Repeat for mock peaks (only 976 peaks by the way, so subsetting is pointless)
name=ZNF143_mock
grep -v "random" ${name}_summits.bed | grep -v "chrUn" | grep -v "chrEBV" | grep -v "chrM" | grep -v "alt" | intersectBed -v -a - -b $blacklist > ${name}_summits_final.bed
slopBed -b 50 -i ${name}_summits_final.bed -g $sizes  | sort -k1,1 -k2,2n > ${name}_summit_100window.bed
fastaFromBed -fi $genome -bed ${name}_summit_100window.bed -fo ${name}_summit_100window.fasta
sort -nrk5,5 ${name}_summit_100window.bed | head -n 1000 > ${name}_top1000_summit_100window.bed
fastaFromBed -fi $genome -bed ${name}_top1000_summit_100window.bed -fo ${name}_top1000_summit_100window.fasta
