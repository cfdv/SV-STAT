#!/bin/bash
# purpose: Identifies amount of support for each candidate junction given stacked reads
# author: Caleb Forbes Davis V
# date: 10/06/11
# external requires:
# -bwa
# -samtools
#
# internal requires:
# -hypoDBgen.pl
#
# assumes:
# - you're looking for t(4;11), t(12;21), t(9;22), or t(1;19) canonical pediatric ALL markers
#
# generates:
#                  derAB.*	- bwa, samtools, and fasta files for candidate junctions search
#		     1.bed	- annotation for all hits to candidate junctions, regardless of map quality
#          raw_cj_supp.txt	- amounts of junction support contributed by each read
#            cj_report.txt	- summary of candidate junction support metrics ( junctionid, sum_bp, sum_qual, sum_bp*qual, log10(sum_bp*qual) )
#

# die if abused
if [[ ( -e $1 ) ]]
    then
	srtx=$(readlink -f $1)
	srfq=$(readlink -f $2)
    else
        echo "sr2sv stackedreads.txt stackedreads.fastq"
        exit 1
fi

### align stacking reads to candidate junctions
#generate candidate junctions as probes for the array
for d in 4_11 12_21 9_22 1_19; do 
	c1=${d%_*}; 
	c2=${d#*_}; 
	cd $d; 
	#report the type of translocation 
        echo $d;
	hypoDBgen hg18 chr$c1 chr$c2 1100 $srtx derAB.fasta;
	#index the candidate junction probes
	bwa index -a bwtsw derAB.fasta;
	#index the candidates for samtools
	samtools faidx derAB.fasta;
	#align to candidate array
	bwa bwasw derAB.fasta $srfq > 0.sam;
	#create the bam
	samtools view -bt derAB.fasta.fai -o 0.bam 0.sam
	#sorting
	samtools sort 0.bam 1;
	#indexing
	samtools index 1.bam;
	#generate text file with full report, and show top hits in stdout
	bamToBed -i 1.bam | tee 1.bed | awk '{cjID=$1; start=$2; end=$3; mq=$5; lb=501-start; rb=end-500; if (lb<rb) minb=lb; else minb=rb; if (minb < 0) minb = 0; print $1,$2,$3,$4,$5, minb, mq, minb*mq;}' | tee raw_cj_supp.txt | awk '{cjID=$1; b=$6; q=$7; bq=$8; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq;}END{for (cjID in bs){if (bqs[cjID]!=0) print cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10);};}' | sort +3 -4 -rn > cj_report.txt;
	cd ..;
	#report top 10 hits to user
	head cj_report.txt;
done;
