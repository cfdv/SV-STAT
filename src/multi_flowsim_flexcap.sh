#!/bin/bash
# Use multiple processors to simulate 454 reads

# die if abused
if [[ ( -e $1 ) && ( $2 -gt 0 ) && ( $3 -gt 0 ) ]]
    then
	f=$(readlink -f $1)
        let numreads=$2
        let numprocs=$3
    else
        echo "multi_flowsim fasta #reads #processors"
        exit 1
fi

let this_proc_numreads=$(( numreads / $(( numprocs )) ))
for (( i=1; i<=numprocs; i++ )); do mkdir $i; cd $i; { flexflowsim_flexcap $f $this_proc_numreads & }; cd ..; done
