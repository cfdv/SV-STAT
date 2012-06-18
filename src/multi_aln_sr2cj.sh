#!/bin/bash
# purpose: Reads without full-length alignments to the reference genome are aligned to
#          candidate junctions. Supports parallel processing
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
# versions:
# v1: controls number of alignment threads used per run
# v2: removes dependence on /_FASTA directory structure, now handled upstream

# die if abused
if [[ ( $# -gt 0 ) && ( $1 -gt 0 ) && ( -e $2 ) ]]
    then
        let numprocs=$1
	srfq=$(readlink -f $2)
    else
        echo "multi_aln_sr2cj max_procs stacked_reads.fastq"
        exit 0
fi


### align stacking reads to candidate junctions
#repeat for each candidate junction fasta file
#keep track of the process ids of each alignment
declare -a pids=();
for d in *; do
	#how many running processes do we have?
	sys_procs=${#pids[@]};
	while [ $sys_procs -eq $numprocs ]; do
		#wait a minute...
		sleep 60;
		#are all the processes still running?
		for i in `seq 0 $(( $sys_procs-1 ))`; do
			#ps gives us 2 lines of output if the process is still running
			#otherwise it returns one line
			stopped=`ps -p "${pids[$i]}" | wc -l`;
			if [ $stopped -eq 1 ];
				then
					unset pids[$i];
			fi
		done;
		sys_procs=${#pids[@]};
	done;

	cd $d;
	#run the job
	{ aln_sr2cj $srfq & };
	#add job's processID to list of running processes
	pids=("${pids[@]}" $!);
	cd ..;
done;
