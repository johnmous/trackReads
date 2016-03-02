#!/bin/bash
# Title: A script to track reads mapped on transcriptome to their location on the genome
# Author: Ioannis Moustakas

# Read from command line
# the gene ID I want to trace
GeneID=$1

# Use the transcriptome file fasta headers info to build a list of tarnscriptIDs
TranscriptIDs=$(grep $GeneID ../Results/HumanGenome_GRCh38_release83/Homo_sapiens.GRCh38.cdna.all.NoAlternates_FixIDs.fa |cut -d ' ' -f1|cut -d '>' -f 2 )

grep -m 1 $GeneID ../Results/HumanGenome_GRCh38_release83/Homo_sapiens.GRCh38.cdna.all.NoAlternates_FixIDs.fa

TranscriptIDs=( $TranscriptIDs )

echo "Number of Transcripts:" ${#TranscriptIDs[@]}

allReadIDs=$(\
for TrID in $TranscriptIDs;
do
  grep $TrID ../Scratch/mapped_Ioannis/S21_tmap_transc.sam | grep -v "@"|cut -f 1 | sort | uniq 
done\
)

countReadIDs=( $allReadIDs )
echo "Number of reads: " ${#countReadIDs[@]}

echo -e "$allReadIDs" | ./readID2GenomicLocation.py ../Scratch/tr_exonmapping_Ioannis/S21_STAR_bowtie.bam  \
../Results/HumanGenome_GRCh38_release83/Homo_sapiens.GRCh38.83.chr.gtf_GeneOnlySimplified.txt