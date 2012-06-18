#!/bin/bash
# Use multiple processors to simulate 454 reads

# die if abused
if [[ ( $# -gt 0 ) && ( $1 -gt 0 ) && ( $2 -gt 0 ) ]]
    then
        let numreads=$1
        let numprocs=$2
    else
        echo "multi_flowsim_cap_nokit #reads #processors"
        exit 0
fi

let this_proc_numreads=$(( numreads / $(( numprocs )) ))
for (( i=1; i<=numprocs; i++ )); do mkdir $i; cd $i; { flexflowsim_cap_nokit $this_proc_numreads & }; cd ..; done
