#!/bin/bash
# requires:
# - in2out
# - FastqToTbl

#INPUT: STDIN - expects raw, flowsim generated fastq. Assumes defline = "@readid qclip start..end"
#OUTPUT: STDOUT - expect adapters trimmed off, quality scores are rescaled, and duplicate start coordinate+orientation reads are removed.
#v2 - no sanity checking, too lossy

# trim adapters | rescale quality scores | remove duplicates | report.fastq
FastqToTbl /dev/stdin | awk '{split($5,x,/\.\./); l=x[2]-x[1]+1; print "@"$1"\n"substr($2,x[1],l)"\n+\n"substr($6,x[1],l);}' | in2out fastq-flowsim fastq | FastqToTbl | sort +0 -1 -u | awk '{print "@"$1"\n"$2"\n"$3"\n"$4;}'
