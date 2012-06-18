#!/bin/bash
# purpose: Generates BAM file from input fastq from flowsim. Uses flowsim-specific quality control procedures. 
# note: not tested with non-flowsim fastq
# author: Caleb Forbes Davis V
# date: 10/06/11
#
# external requires:
# -bwa
# -samtools
# -java
#
# internal requires:
# -flowsim_qc
# -FastqToTbl
#
# assumes:
# -samtools index for reference is "ref.fa.fai"
# -path to MarkDuplicates = ~/src/picard-tools-1.40/MarkDuplicates.jar
#
# generates:
#            pass_qc.fastq	- adapters trimmed, base qualities scaled to match empirical distributions, duplicates (by ori and start coord) removed
#       pass_qc.fastq.cidx	- cdbfasta index of pass_qc.fastq
#             read_len.txt	- readid, read_length
#                    0.sam	- raw alignment
#                    1.sam	- removed 0-length alignments, if any
#                    2.bam 	- bam for 1.sam
#                    3.bam 	- sorted 2.bam
#                    4.bam	- 3.bam filtered by picard to remove duplicates
#                4.bam.bai	- samtools index of 4.bam

# die if abused
if [[ ( -e $1 ) ]]
    then
	f=$(readlink -f $1)
    else
        echo "flowsim2bam ref.fa"
        exit 0
fi

###### make BAM file from input fastq
#expects fastq input, presumably result of flowsim simulation of reads from a single sample
flowsim_qc | tee pass_qc.fastq | FastqToTbl | awk '{print $1, length($2);}' > read_len.txt;
#make index of passed_qc reads for later
cdbfasta -Q pass_qc.fastq;

#align to ref
bwa bwasw $f pass_qc.fastq > 0.sam;
#remove 0-length alignments, if any
cat <(grep "^@" 0.sam) <(grep -v "^@" 0.sam | awk '{if (NF >= 11) print $0;}') > 1.sam;
#create bam file
samtools view -bt $f.fai -o 2.bam 1.sam;
#sort bam file
samtools sort 2.bam 3;
#remove duplicate reads
java -jar ~/src/picard-tools-1.40/MarkDuplicates.jar INPUT=3.bam OUTPUT=4.bam METRICS_FILE=4_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true;
#index the bam file for samtools
samtools index 4.bam;
