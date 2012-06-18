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
# versions:
# -v1: executes hypoDBgen for each pair of regions connected by putative SV
# -v2: moves to /_FASTA before aligning stacked reads to SV probes, reducing cohesion

### align stacking reads to candidate junctions
#generate candidate junctions as probes for the array
for d in 4_11 12_21 9_22 1_19; do 
	c1=${d%_*}; 
	c2=${d#*_}; 
	cd $d; 
	#report the type of translocation 
        echo $d;

	# die if abused
	if [[ ( -e $1 )  && ( -e $2 ) ]]
	    then
		srtx=$(readlink -f $1)
		srfq=$(readlink -f $2)
	    else
		echo "sr2sv stackedreads.txt stackedreads.fastq"
		exit 0
	fi

	#dedicate folder to candidate junction fasta files
	mkdir _FASTA;
	cd _FASTA;
	#generate candidate junctions
	hypoDBgen hg18 chr$c1 chr$c2 1100 ../stackedreads.txt derAB.fasta;
	#perform alignments of stacked reads fastq against candidate junctions
	#multi_aln_sr2cj.sh #processors stackedreads.fastq
	multi_aln_sr2cj 2 ../stackedreads.txt.fastq;
	cd ../..;
done;
