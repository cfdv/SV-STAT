#!/bin/bash
#Generates reads from three distributions, hoping to match empirical results
#accepts one parameter - the number of desired reads to simulate
#v1 - 3 distributions, comparison looks like H:\cfdavis\Lau\LeukSeqCap\results\96C_vs_flexflowsim_readlength_distributions.jpg
#v2 - 4 distributions, added a lognormal
#v3 - adds ... | mutator | ...

#die if abused
if [[ ( $# -gt 0 ) && ( $1 -gt 0 ) ]]
    then
        let numreads=$1
    else
        echo "flexflowsim number_of_reads"
        exit 1
fi

# set the weights of the 3 distributions. 
#1 6/1/3 gauss/uniform/lognorm
#2 6/3/1 gauss/uniform/lognorm
#3 55/35/10 gauss/uniform/lognorm
#4 55/30/15 gauss/uniform/lognorm
let lognorm_numreads=$(( numreads * 8 / 100 ))
let uniform_numreads=$(( numreads * 27 / 100 ))
let gauss_numreads=$(( numreads * 55 / 100 ))
let uniform_mid_numreads=$(( numreads * 10 / 100 ))	#another log-normal distribution for mid-range read lengths

#make the LogNorm component
mkdir lognorm
cd lognorm
##1 clonesim -c $lognorm_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="LogNormal 3.9 0.3" | kitsim | flowsim --generation=Titanium -o output.sff
##2 clonesim -c $lognorm_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="LogNormal 3.9 0.4" | kitsim | flowsim --generation=Titanium -o output.sff
##3 clonesim -c $lognorm_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="LogNormal 3.9 0.4" | kitsim | flowsim --generation=Titanium -o output.sff
##4 ""
clonesim -c $lognorm_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="LogNormal 3.9 0.2" | kitsim | mutator | flowsim --generation=Titanium -o output.sff
flower output.sff --fasta=output.fasta --fastq=output.fastq
grep ">" output.fasta | sed 's/>[^\.]*\.\.\(.*\)/\1/' > readlengths.txt

#make the uniform component
mkdir ../uniform
cd ../uniform
##1 clonesim -c $uniform_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Uniform 85 400" | kitsim | flowsim --generation=Titanium -o output.sff
##2 clonesim -c $uniform_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Uniform 85 400" | kitsim | flowsim --generation=Titanium -o output.sff
##3 clonesim -c $uniform_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Uniform 75 400" | kitsim | flowsim --generation=Titanium -o output.sff
##4 ""
clonesim -c $uniform_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Uniform 65 400" | kitsim | mutator | flowsim --generation=Titanium -o output.sff
flower output.sff --fasta=output.fasta --fastq=output.fastq
grep ">" output.fasta | sed 's/>[^\.]*\.\.\(.*\)/\1/' > readlengths.txt

#make the gaussian component
mkdir ../gauss
cd ../gauss
##1 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta | kitsim | flowsim --generation=Titanium --degradation="Normal 0.00275 0.001" -o output.sff
##2 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Normal 390 55" | kitsim | flowsim --generation=Titanium -o output.sff
##3 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Normal 390 55" | kitsim | flowsim --generation=Titanium -o output.sff
##4 ""
clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Normal 390 55" | kitsim | mutator | flowsim --generation=Titanium -o output.sff
flower output.sff --fasta=output.fasta --fastq=output.fastq
grep ">" output.fasta | sed 's/>[^\.]*\.\.\(.*\)/\1/' > readlengths.txt

#make the midrange lognormal component
mkdir ../midlognormal
cd ../midlognormal
##1 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta | kitsim | flowsim --generation=Titanium --degradation="Normal 0.00275 0.001" -o output.sff
##2 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Normal 390 55" | kitsim | flowsim --generation=Titanium -o output.sff
##3 clonesim -c $gauss_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="Normal 390 55" | kitsim | flowsim --generation=Titanium -o output.sff
##4 ""
clonesim -c $uniform_mid_numreads ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta --lengths="LogNormal 5.7 0.3" | kitsim | mutator | flowsim --generation=Titanium -o output.sff
flower output.sff --fasta=output.fasta --fastq=output.fastq
grep ">" output.fasta | sed 's/>[^\.]*\.\.\(.*\)/\1/' > readlengths.txt

#return to initial directory
cd ..
