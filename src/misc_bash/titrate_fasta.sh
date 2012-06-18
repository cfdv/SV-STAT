#!/bin/bash
# purpose: Titrates two genome fasta files in the desired proportion.
#          useful for diluting tumor in normal for simulations 
# author: Caleb Forbes Davis V
# date: 11/2/11
#
# assumes: 
#  - normal sample fasta file is haploid

# die if abused
if [[ ( $# -gt 0 ) && ( -e $1 ) && ( $2 -gt 0 ) && ( -e $3 ) && ( $4 -gt 0 ) ]]
    then
	samp1=$(readlink -f $1)
	samp2=$(readlink -f $3)
	num1=$2
	num2=$4
    else
        echo "titrate_fasta /path/to/tumor #tumor_copies /path/to/normal #normal_copies"
        exit 0
fi

for i in `seq 1 $num1`; do cat $samp1; done
#do 2x the #normal copies, assuming the existing normal fasta file is haploid
for i in `seq 1 $(($num2*2))`; do cat $samp2; done
