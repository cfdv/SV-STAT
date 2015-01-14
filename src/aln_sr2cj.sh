#!/bin/bash
# purpose: Reads without full-length alignments to the reference genome are aligned to
#          candidate junctions. 
# author: Caleb Forbes Davis V
# date: 10/13/11
# external requires:
# -bwa
# -samtools
#
# internal requires:
#
# assumes:
#
# generates:
#        _FASTA/derAB_##.*	- bwa, samtools, and fasta files for candidate junctions search
#		     1.bed	- annotation for all hits to candidate junctions, regardless of map quality
#          raw_cj_supp.txt	- amounts of junction support contributed by each read
#            cj_report.txt	- summary of candidate junction support metrics ( junctionid, sum_bp, sum_qual, sum_bp*qual, log10(sum_bp*qual) )
#

# die if abused
if [[ ( $# -gt 0 ) && ( -e $1 ) ]]
    then
	srfq=$(readlink -f $1)
    else
        echo "aln_sr2cj stacked_reads.fastq"
        exit 1
fi

#index the candidate junction probes
$BWA_BIN_DIR/bwa index -a bwtsw derAB.fasta;
#index the candidates for samtools
samtools faidx derAB.fasta;
#align to candidate array
$BWA_BIN_DIR/bwa bwasw derAB.fasta $srfq > 0.sam;
#create the bam
samtools view -bt derAB.fasta.fai -o 0.bam 0.sam;
#sorting
samtools sort 0.bam 1;
#indexing
samtools index 1.bam;
#generate text file with full report, and show top hits in stdout. tee -i is an attempt to avoid truncating output
$BEDTOOLS_BIN_DIR/bedtools bamtobed -i 1.bam | tee -i 1.bed | awk '{cjID=$1; start=$2; end=$3; mq=$5; lb=501-start; rb=end-500; if (lb<rb) minb=lb; else minb=rb; if (minb < 0) minb = 0; print $1,$2,$3,$4,$5, minb, mq, minb*mq;}' | tee -i raw_cj_supp.txt | awk '{cjID=$1; b=$6; q=$7; bq=$8; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq;}END{for (cjID in bs){if (bqs[cjID]!=0) print cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10);};}' | sort +3 -4 -rn > cj_report.txt;
#report top 10 hits to user
head cj_report.txt;
