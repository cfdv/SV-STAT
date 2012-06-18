#!/bin/bash
clonesim -c 350000 ~/goldenpath/hg18/leukseqcap/sample4_capture_target.fasta | kitsim | flowsim --generation=Titanium --degradation="Normal 0.00275 0.001" -o output_14e5_Titanium_deg_00275.sff
flower output_14e5_Titanium_deg_00275.sff --fasta=14e5_Titanium_deg_00275.fasta --fastq=14e5_Titanium_deg_00275.fastq
grep ">" 14e5_Titanium_deg_00275.fasta | sed 's/>[^\.]*\.\.\(.*\)/\1/' > 14e5_Titanium_deg_00275_readlengths.txt
