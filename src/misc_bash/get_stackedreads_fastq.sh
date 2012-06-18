#!/bin/bash
# purpose: Generates fastq for stacked reads
# author: Caleb Forbes Davis V
# date: 10/06/11
#
# assumes:
# -folder and file structure like t/stackedreads.txt, where t=translocation model.
# -reads are stored in pass_qc.fastq.cidx in the current directory
# -readids in t/stackedreads.txt match the keys in the fastq index
#
# generates:
# t/stackedreads.fastq	- stacked reads fastq

#grab fastq of each stacked read with a unique combination of chromosome, start coordinate, and orientation
for file in */stackedreads.txt; do cat $file | cut -d" " -f1 | sort | uniq | cdbyank pass_qc.fastq.cidx > stackedreads.fastq; done;
