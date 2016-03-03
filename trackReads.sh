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
for TrID in ${TranscriptIDs[@]};
do
  grep $TrID ../Scratch/mapped_Ioannis/S21_tmap_transc.sam | grep -v "@" | cut -f 1 | sort | uniq
  # calculate mismatch statistics on the trascriptome file. 
  # Redirect the STDOUT to STDERR so it is not captured in the variable
  echo -e "\nStats on Reads of Transcript: $TrID" 1>&2 
  grep $TrID ../Scratch/mapped_Ioannis/S21_tmap_transc.sam | grep -v "@" |\
  perl -e 'my $len; 
  $mismatch; 
  my $i=0; 
  while (<STDIN>) {
    $i++; 
    my @field=split("\t", $_); 
    $len+=length(@field[9]); 
    $_ =~ m/NM:i:([0-9]+)/;
    $mismatch+=$1;
    }
  if ($i==0) {print(" Number of Reads: ", $i, "\n")
    }
  else {
    printf(" Number of Reads: %d \n Average MisMatch/ReadLength: %.3f \n Average MisMatch Number: %.1f \n", $i,$mismatch/$len, $mismatch/$i)
    }' 1>&2 
    
done\
)

countReadIDs=( $allReadIDs )
echo -e "\nNumber of reads in all Transcripts: " ${#countReadIDs[@]}

echo -e "$allReadIDs" | ./readID2GenomicLocation.py ../Scratch/tr_exonmapping_Ioannis/S21_STAR_bowtie.bam  \
../Results/HumanGenome_GRCh38_release83/Homo_sapiens.GRCh38.83.chr.gtf_GeneOnlySimplified.txt