#!/bin/bash
#Structural variants (SVs) aberrantly connect genomic regions to each other. Existing tools use different approaches
# to solve problems related to observations of nucleotide sequences shorter than the repetitive regions connected by
# SVs. Assuming the user suspects fusions in specific orientations between specific regions, SV-STAT addresses these
# problems in three steps. 1) Create a synthetic reference of candidate SVs, 2) align misaligned reads to the synthetic
# reference, and 3) score support for SVs provided by what could now be correctly aligned reads.
#author - Caleb Davis
#date - 11/13/13
#purpose - test all the things

# die if abused
if [[ ( $# -ne 0 ) ]]
    then
        echo "test-svstat"
        exit 1
fi

MYEXITCODE=0

# eventually allow alternate test data. For now, assume "bam" and "metadata" and "result" folders are co-located
MYTESTPATH=`readlink -m $0 | xargs dirname`

# Some functionality depends on the versions of various tools
MYSAMTOOLSVER="(0.1.17+)"

# Is samtools in the user's path?
NOSAMTOOLS=`which samtools 2> /dev/stdout | grep "no samtools in" | wc -l`
if [ "$NOSAMTOOLS" -ne 0 ]; then echo "FATAL: install samtools $MYSAMTOOLSVER in \$PATH"; MYEXITCODE=1; fi

# Does samtools respond as expected?
SAMTOOLSRESPONSE=`samtools view 2> /dev/stdout | grep BAM | wc -l`
if [ "$SAMTOOLSRESPONSE" -eq 0 ]; then echo -e "FATAL: fix samtools installation"; MYEXITCODE=1; fi

# Is samtools loaded with bam2fq?
NOBAM2FQ=`samtools bam2fq $MYTESTPATH/bam/hs1011.bam 2>&1 | grep unrecognized | wc -l`
if [ "$NOBAM2FQ" -ne 0 ]; then echo "FATAL: upgrade samtools $MYSAMTOOLSVER"; MYEXITCODE=1; fi

exit $MYEXITCODE