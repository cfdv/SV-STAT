#!/bin/bash
# purpose: Generates stacked reads from input bam alignment. Stacked reads format is compatible with SVSTAT
# author: Caleb Forbes Davis V
# date: 10/06/11
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
# - you're looking for t(4;11), t(12;21), t(9;22), or t(1;19) canonical pediatric ALL markers
# - target coordinates:  hg18 { chr1:162870000-163070000
# 				chr11:117815000-117915000
# 				chr12:11871000-11971000
# 				chr19:1544000-1594000
# 				chr21:35135000-35385000
# 				chr22:21790000-21990000
# 				chr4:88095000-88295000
# 				chr9:132560000-132760000 }
# generates:
#  onehit_stackedreads.txt	- all partially aligned reads
#         rstack_coord.txt	- stacks with tails pointing upstream relative to reference
#         fstack_coord.txt	- stacks with tails pointing downstream
#  {4_11,9_22,12_21,1_19}/	- translocation-specific workspaces (t/)
#       t/stackedreads.txt	- read alignment annotations joined to their stacks, compatible with SVSTAT

# die if abused
if [[ ( -e $1 ) && ( -e $2 ) ]]
    then
	b=$(readlink -f $1)
	rl=$(readlink -f $2)
    else
        echo "bam2stacks bam read_lengths"
        exit 1
fi

#assumes these are the translocations you're looking for
mkdir 4_11 12_21 9_22 1_19;
#find partial alignments sharing start or end coordinates (stacks)
join -j1 -o 2.1,2.2,1.2,2.3,2.4,2.5,2.6,2.7,2.8 <(sort +0 -1 $rl) <(bamtobl.pl $b | sort +0 -1) | awk 'BEGIN{plchldr_cols="1.2 1.4 1.5 1.6 1.7 1.8 1.9"}{ss=$6; se=$7; qs=$4; qe=$5; qlen=$3; ori="+"; if ($8!=1) {ori="-";tmp=qs; qs=qlen-qe+1; qe=qlen-tmp+1;} if (qs!=1) print $0,plchldr_cols,"start"; if (qe!=qlen) print $0,plchldr_cols,"end";}' | tee onehit_stackedreads.txt | awk '{if ($17=="start") print $2,$6,$17; else print $2,$7,$17;}' | sort | uniq -c | awk '{if ($1>1) print $2,$3,$4;}' | tee >(grep end$ > fstack_coord.txt) | grep start$ > rstack_coord.txt;
#join stacks with their corresponding reads, and output result to translocation-specific workspaces
cat <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat fstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep end$ onehit_stackedreads.txt | awk '{print $2$7$17,$0}' | sort +0 -1)) <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat rstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep start$ onehit_stackedreads.txt | awk '{print $2$6$17,$0}' | sort +0 -1)) | awk '{if (match($2,/chr/)==1) print $0; else {$2="chr"$2;print $0;}}' | awk '{ori="+"; if ($8!=1) ori="-"; print $1,$2,"-1","-1","-1","-1",$3,$4,$5,$6,$7,ori,"-1","-1","-1_-1",$17,$9,$10,$11,$12,$13,$14,$15,$16;}' | tee >(awk '{if (($2=="chr11" && (($10>117815000 && $10<117915000) || ($11>117815000 && $11<117915000))) || ($2=="chr4" && (($10>88095000 && $10<88295000) || ($11>88095000 && $11<88295000)))) print $0;}' > 4_11/stackedreads.txt) |  tee >(awk '{if (($2=="chr21" && (($10>35135000 && $10<35385000) || ($11>35135000 && $11<35385000))) || ($2=="chr12" && (($10>11871000 && $10<11971000) || ($11>11871000 && $11<11971000)))) print $0;}' > 12_21/stackedreads.txt) | tee >(awk '{if (($2=="chr1" && (($10>162870000 && $10<163070000) || ($11>162870000 && $11<163070000))) || ($2=="chr19" && (($10>1544000 && $10<1594000) || ($11>1544000 && $11<1594000)))) print $0;}' > 1_19/stackedreads.txt) | awk '{if (($2=="chr9" && (($10>132560000 && $10<132760000) || ($11>132560000 && $11<132760000))) || ($2=="chr22" && (($10>21790000 && $10<21990000) || ($11>21790000 && $11<21990000)))) print $0;}' > 9_22/stackedreads.txt;
