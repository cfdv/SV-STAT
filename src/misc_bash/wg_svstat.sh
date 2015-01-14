#!/bin/bash
# whole-genome SVSTAT. for now, just export the stacked reads, please
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
# -samtools index is "ref.fa.fai"
# -target coordinates =  hg18 { chr1:162870000-163070000
# 				chr11:117815000-117915000
# 				chr12:11871000-11971000
# 				chr19:1544000-1594000
# 				chr21:35135000-35385000
# 				chr22:21790000-21990000
# 				chr4:88095000-88295000
# 				chr9:132560000-132760000 }
# -path to MarkDuplicates = ~/src/picard-tools-1.40/MarkDuplicates.jar

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
#  onehit_stackedreads.txt	- all partially aligned reads
#         rstack_coord.txt	- stacks with tails pointing upstream relative to reference
#         fstack_coord.txt	- stacks with tails pointing downstream
#  {4_11,9_22,12_21,1_19}/	- translocation-specific workspaces (t/)
#       t/stackedreads.txt	- read alignment annotations joined to their stacks, input to SVSTAT
# t/stackedreads.txt.fastq	- stacked reads fastq
#                  derAB.*	- bwa, samtools, and fasta files for candidate junctions search
#            cj_report.txt	- result!
#
# versions:

# die if abused
if [[ ( -e $1 ) ]]
    then
	f=$(readlink -f $1)
    else
        echo "svstat ref.fa"
        exit 1
fi

#expects fastq input, presumably result of flowsim simulation of reads from a single sample
#flowsim_qc | tee pass_qc.fastq | FastqToTbl | awk '{print $1, length($2);}' > read_len.txt;
#make index of passed_qc reads for later
#cdbfasta -Q pass_qc.fastq;

#prepare resources
#cat /dev/stdin > in.fastq;
#cdbfasta -Q in.fastq;

#align to ref
#bwa bwasw $f in.fastq > 0.sam;
#remove 0-length alignments, if any
#cat <(grep "^@" 0.sam) <(grep -v "^@" 0.sam | awk '{if (NF >= 11) print $0;}') > 1.sam;
#create bam file
#samtools view -bt $f.fai -o 2.bam 1.sam;
#sort bam file
#samtools sort 2.bam 3;
#remove duplicate reads
#java -jar ~/src/picard-tools-1.40/MarkDuplicates.jar INPUT=3.bam OUTPUT=4.bam METRICS_FILE=4_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true;
#index the bam file for samtools
#samtools index 4.bam;

############### SVSTAT
### find stacks
#assumes these are the translocations you're looking for
#mkdir 4_11 12_21 9_22 1_19;
#find partial alignments sharing start or end coordinates (stacks)
join -j1 -o 2.1,2.2,1.2,2.3,2.4,2.5,2.6,2.7,2.8 <(sort +0 -1 read_len.txt) <(bamtobl.pl 4.bam | sort +0 -1) | awk 'BEGIN{plchldr_cols="1.2 1.4 1.5 1.6 1.7 1.8 1.9"}{ss=$6; se=$7; qs=$4; qe=$5; qlen=$3; ori="+"; if ($8!=1) {ori="-";tmp=qs; qs=qlen-qe+1; qe=qlen-tmp+1;} if (qs!=1) print $0,plchldr_cols,"start"; if (qe!=qlen) print $0,plchldr_cols,"end";}' | tee onehit_stackedreads.txt | awk '{if ($17=="start") print $2,$6,$17; else print $2,$7,$17;}' | sort | uniq -c | awk '{if ($1>1) print $2,$3,$4;}' | tee >(grep end$ > fstack_coord.txt) | grep start$ > rstack_coord.txt;
#join stacks with their corresponding reads, and output result to translocation-specific workspaces
cat <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat fstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep end$ onehit_stackedreads.txt | awk '{print $2$7$17,$0}' | sort +0 -1)) <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat rstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep start$ onehit_stackedreads.txt | awk '{print $2$6$17,$0}' | sort +0 -1)) | awk '{if (match($2,/chr/)==1) print $0; else {$2="chr"$2;print $0;}}' | awk '{ori="+"; if ($8!=1) ori="-"; print $1,$2,"-1","-1","-1","-1",$3,$4,$5,$6,$7,ori,"-1","-1","-1_-1",$17,$9,$10,$11,$12,$13,$14,$15,$16;}' > stackedreads.txt;
#grab fastq of stacked reads
cat stackedreads.txt | cut -d" " -f1 | sort | uniq | cdbyank in.fastq.cidx > stackedreads.fastq;

### align stacking reads to candidate junctions
#generate candidate junctions as probes for the array
#for d in 4_11 12_21 9_22 1_19; do 
#	c1=${d%_*}; 
#	c2=${d#*_}; 
#	cd $d; 
#	#report the type of translocation 
#        echo $d;
#	hypoDBgen hg18 chr$c1 chr$c2 1100 stackedreads.txt derAB.fasta;
#	#index the candidate junction probes
#	bwa index -a bwtsw derAB.fasta;
#	#index the candidates for samtools
#	samtools faidx derAB.fasta;
#	#align to candidate array
#	bwa bwasw derAB.fasta stackedreads.txt.fastq > 0.sam;
#	#create the bam
#	samtools view -bt derAB.fasta.fai -o 0.bam 0.sam
#	#sorting
#	samtools sort 0.bam 1;
#	#indexing
#	samtools index 1.bam;
#	#generate text file with full report, and show top hits in stdout
#	bamToBed -i 1.bam | tee 1.bed | awk '{cjID=$1; start=$2; end=$3; mq=$5; lb=501-start; rb=end-500; if (lb<rb) minb=lb; else minb=rb; if (minb < 0) minb = 0; print $1,$2,$3,$4,$5, minb, mq, minb*mq;}' | tee raw_cj_supp.txt | awk '{cjID=$1; b=$6; q=$7; bq=$8; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq;}END{for (cjID in bs){if (bqs[cjID]!=0) print cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10);};}' | sort +3 -4 -rn | tee cj_report.txt | head;
#	cd ..;
#done;
