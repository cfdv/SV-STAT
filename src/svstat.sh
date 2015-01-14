#!/bin/bash
#Structural variants (SVs) aberrantly connect genomic regions to each other. Existing tools use different approaches
# to solve problems related to observations of nucleotide sequences shorter than the repetitive regions connected by
# SVs. Assuming the user suspects fusions in specific orientations between specific regions, SV-STAT addresses these
# problems in three steps. 1) Create a synthetic reference of candidate SVs, 2) align misaligned reads to the synthetic
# reference, and 3) score support for SVs provided by what could now be correctly aligned reads.
#author - Caleb Davis
#date - 2/28/11
#input - bam file location, and candidate SVs
#output format - candidate_junctionID readID chrom qstart qend sstart send ori BWAmapscore
#purpose - generate .bl-like alignment coordinates of reads mismapped to the reference but correctly mapped to junctions of SVs
# external requires:
# -samtools
# -java/picard-tools
# -cdbfasta
#
# internal requires:
# -FastqToTbl
#
# assumes:
# -up to ~ 1 Mb total uncertainty in positions of breakpoints for each SV
# -sorted bam
# -environmental variables
#  * CDBFASTA_BIN_DIR - path to cdbfasta
#  * SVSTAT_SRC_DIR   - path to svstat executables, e.g. FastqToTbl
#  * PICARD_JAR_DIR   - path to picard-tools jar files

# generates:
#               pass_qc.fastq	- adapters trimmed, base qualities scaled to match empirical distributions, duplicates (by ori and start coord) removed
#          pass_qc.fastq.cidx	- cdbfasta index of pass_qc.fastq
#                read_len.txt	- readid, read_length
#                       4.bam	- 3.bam filtered by picard to remove duplicates
#                   4.bam.bai	- samtools index of 4.bam
#     onehit_stackedreads.txt	- all partially aligned reads
#            rstack_coord.txt	- stacks with tails pointing upstream relative to reference
#            fstack_coord.txt	- stacks with tails pointing downstream
#        {SV_1,SV_2,...,SV_n}	- SV-specific workspaces (SVID/)
#       SVID/stackedreads.txt	- read alignment annotations joined to their stacks, input to SVSTAT
# SVID/stackedreads.txt.fastq	- stacked reads fastq
#      SVID/_FASTA/derAB_##.*	- bwa, samtools, and fasta files for candidate junctions search
#          SVID/cj_report.txt	- result!
#
# versions:
# -v1: accepts raw flowsim fastq, and reports correct junctions in sample 4 from among all possible translocation types
# -v2: works with latest version of hypoDBgen to use parallel processing for alignments of stacked reads to candidate junctions
# -v3: eliminate all the stacks caused by same-start+orientation reads
# -v4: read SV candidates from file, and test behavior for one well-understood non-B-ALL SV candidate
# -v5: TODO - allow arbitrary SV candidates

usage() {
    echo
    head -n2 $XDIR/../README | cut -c3-
    echo "
svstat.sh [-hma] [-d SCRATCHPATH] <bam> <candidate SVs>

  [ -h ]              print this help and exit
  [ -m ]              also provide output in a vcf-like format
  [ -a ]              include low-quality junctions in the vcf-like output
  [ -d SCRATCHPATH]   use a temporary workspace (default is `mktemp | xargs dirname`)
"
}

testme() {
    $XDIR/../test/test-svstat.sh
    EXITCODE=$?
    if [ "$TMPSVNAME" != 0 ]; then exit $EXITCODE; fi
}

metagen() {
    usage | awk 'NR==3{print "##program=svstat-" $2;}'
    echo "##META: path to SV-STAT=$XDIR"
    echo "##input_bam=$(readlink -m $MYBAM)"
    echo "##input_candidate_SVs=$(readlink -m $MYMODELS)"
    echo "##parameters=MAXCVG:${MAXCVG} SLOPE:${SLOPE} INTERCEPT:${INTERCEPT} MAXSUPPREADS:${MAXSUPPREADS}"
    echo "##META: all SV-STAT INFO lines will start with the letters VT"
    echo "##META: hg19 coordinates"
    echo "##INFO=<ID=\"SVPNAME\",Number=.,Type=String,Description=\"unique name given to this event (for uniting with viz info)\">"
    echo "##INFO=<ID=\"VCT\",Number=.,Type=String,Description=\"SV-Parlaiment:type of variant called. Value is Insertion, Deletion, or Mismatch.\">"
    echo "##INFO=<ID=\"VTO1\",Number=.,Type=String,Description=\"orientation 1\">"
    echo "##INFO=<ID=\"VTO2\",Number=.,Type=String,Description=\"orientation 2\">"
    echo "##INFO=<ID=\"VTBLSC\",Number=.,Type=Float,Description=\"log-support score of soft-clipped reads aligned to synthetic junction with BLAST\">"
    echo "##INFO=<ID=\"VTBLNR\",Number=.,Type=Integer,Description=\"number of reads aligned to synthetic junction with BLAST\">"
    echo "##INFO=<ID=\"VTBWSC\",Number=.,Type=Float,Description=\"log-support score of soft-clipped reads aligned to synthetic junction with BWA\">"
    echo "##INFO=<ID=\"VTBWNR\",Number=.,Type=Integer,Description=\"number of reads aligned to synthetic junction with BWA\">"
    echo "##INFO=<ID=\"VTTY\",Number=.,Type=String,Description=\"SV-STAT:type of SV\">"
    echo "##INFO=<ID=\"VTXP\",Number=.,Type=String,Description=\"inter-chromosomal translocation CHROM:END:VTO1_VTXP:VTO2\">"
    echo "##INFO=<ID=\"VTF\",Number=.,Type=String,Description=\"SV-STAT filter status. Value is PASS or FAIL.\">"
    echo -e "#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tREFSEQ\tVARSEQ\tINFO"
}

# internal variables
XDIR=`dirname $(readlink -m $0)`
ORIGDIR=$PWD
TMPSVNAME=""
MYTMPFILE=`mktemp`
MYTMPDIR=`dirname $(readlink -m $MYTMPFILE)`

while getopts "hmad:" OPTION; do
    case $OPTION in
        h) usage; exit 0;;
        m) META=True;;
        a) FILT=False;;
        d) SCRDIR=$OPTARG;;
    esac
done

# user-accessible parameters with defaults
META=${META:-False}
FILT=${FILT:-True}
SCRDIR=${SCRDIR:-$MYTMPDIR}

shift $((OPTIND-1))

# test install then exit with help if used without valid arguments
if [[ ( -e $1 ) && ( -e $2 ) ]]
    then
        MYBAM=$1
        MYMODELS=$2
    else
        usage
        testme
        exit 1
fi

#algorithm parameters
MAXCVG=50000
SLOPE=1.5149
INTERCEPT=-4.3414
MAXSUPPREADS=160

############### SVSTAT
### find stacks
#TODO: think about allowing an optional -sorted parameter, which may enable operations on a streaming BAM.
#process models file line by line
while read SVNAME FIRSTCHR FIRSTSTART FIRSTEND FIRSTORI SECONDCHR SECONDSTART SECONDEND SECONDORI RECIPROCAL; do
    if [ -d "$SCRDIR/$SVNAME" ]; then continue; fi
    # Make a note of the name of the first SV for this job. We'll use it in the name of a log file
    if [ "$TMPSVNAME" == "" ]; then FIRSTSVNAME=$SVNAME; TMPSVNAME=$SVNAME; fi
    echo "[SV-STAT:svstat.sh] Processing $SVNAME..."
    mkdir -p $SCRDIR/$SVNAME/fastq $SCRDIR/$SVNAME/bam/intermediates $SCRDIR/$SVNAME/bam/full $SCRDIR/$SVNAME/candidates
    MYREGIONS=`echo -en "$FIRSTCHR:$FIRSTSTART-$FIRSTEND $SECONDCHR:$SECONDSTART-$SECONDEND" | sed 's/chr//g'`
    samtools view -hb $MYBAM $MYREGIONS > $SCRDIR/$SVNAME/bam/intermediates/2.bam
    samtools sort $SCRDIR/$SVNAME/bam/intermediates/2.bam $SCRDIR/$SVNAME/bam/intermediates/3
    MYFASTQ=$SCRDIR/$SVNAME/fastq/samtools-bam2fq.fastq
    NUMREADS=`samtools view -hbu $SCRDIR/$SVNAME/bam/intermediates/3.bam | samtools bam2fq - | tee $MYFASTQ | wc -l`
    if (( "$NUMREADS" < "$MAXCVG" ))
        then
            $CDBFASTA_BIN_DIR/cdbfasta -Q $MYFASTQ
            $SVSTAT_SRC_DIR/FastqToTbl $MYFASTQ | awk '{print $1, length($2);}' > $SCRDIR/$SVNAME/read_len.txt
            java -jar $PICARD_JAR_DIR/MarkDuplicates.jar INPUT=$SCRDIR/$SVNAME/bam/intermediates/3.bam OUTPUT=$SCRDIR/$SVNAME/bam/full/4.bam METRICS_FILE=$SCRDIR/$SVNAME/bam/full/4_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
            samtools index $SCRDIR/$SVNAME/bam/full/4.bam
            join -j1 -o 2.1,2.2,1.2,2.3,2.4,2.5,2.6,2.7,2.8 <(sort +0 -1 $SCRDIR/$SVNAME/read_len.txt) <(perl $SVSTAT_SRC_DIR/bamtobl.pl $SCRDIR/$SVNAME/bam/full/4.bam | sort +0 -1) | awk 'BEGIN{plchldr_cols="1.2 1.4 1.5 1.6 1.7 1.8 1.9"}{ss=$6; se=$7; qs=$4; qe=$5; qlen=$3; ori="+"; if ($8!=1) {ori="-";tmp=qs; qs=qlen-qe+1; qe=qlen-tmp+1;} if (qs!=1) print $0,plchldr_cols,"start"; if (qe!=qlen) print $0,plchldr_cols,"end";}' | awk '$8<0{type=$17; if (type~/start/) $17="end"; if (type~/end/) $17="start"; print;}$8>0{print;}' | tee $SCRDIR/$SVNAME/onehit_stackedreads.txt | sort -k1,1 -k17 -u | awk '{if ($17=="start") print $2,$6,$17; else print $2,$7,$17;}' | sort | uniq -c | awk '{if ($1>1) print $2,$3,$4;}' | tee >(grep end$ > $SCRDIR/$SVNAME/fstack_coord.txt) | grep start$ > $SCRDIR/$SVNAME/rstack_coord.txt
            cat <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat $SCRDIR/$SVNAME/fstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep end$ $SCRDIR/$SVNAME/onehit_stackedreads.txt | awk '{print $2$7$17,$0}' | sort +0 -1)) <(join -j1 -o 2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18 <(cat $SCRDIR/$SVNAME/rstack_coord.txt | awk '{print $1$2$3,$0}' | sort +0 -1) <(grep start$ $SCRDIR/$SVNAME/onehit_stackedreads.txt | awk '{print $2$6$17,$0}' | sort +0 -1)) | awk '{if (match($2,/chr/)==1) print $0; else {$2="chr"$2;print $0;}}' | awk '{ori="+"; if ($8!=1) ori="-"; print $1,$2,"-1","-1","-1","-1",$3,$4,$5,$6,$7,ori,"-1","-1","-1_-1",$17,$9,$10,$11,$12,$13,$14,$15,$16;}' > $SCRDIR/$SVNAME/stackedreads.txt
            cut -d" " -f1 $SCRDIR/$SVNAME/stackedreads.txt | sort | uniq | $CDBFASTA_BIN_DIR/cdbyank $SCRDIR/$SVNAME/fastq/samtools-bam2fq.fastq.cidx > $SCRDIR/$SVNAME/stackedreads.txt.fastq
            $SVSTAT_SRC_DIR/FastqToTbl $SCRDIR/$SVNAME/stackedreads.txt.fastq | awk '{print ">"$1"\n"$2;}' > $SCRDIR/$SVNAME/stackedreads.txt.fasta
            #TODO FIX THIS TO choose the correct set of models given AORI and BORI
            cd $SCRDIR/$SVNAME/candidates
            perl $SVSTAT_SRC_DIR/hypoDBgen.pl $REF_DIR $FIRSTCHR $FIRSTORI $SECONDCHR $SECONDORI $RECIPROCAL ../stackedreads.txt derAB.fasta
            #TODO parameterize number of processors for multi_aln_sr2cj.sh -- currently 1
            #$SVSTAT_SRC_DIR/multi_aln_sr2cj.sh 1 ../stackedreads.txt.fastq
            for d in *; do
                if [ -d "$d" ]
                    then
                        cd $d
                        $SVSTAT_SRC_DIR/aln_sr2cj.sh ../../stackedreads.txt.fastq
                        cd ..
                    else
                        echo "[SV-STAT:svstat.sh] no candidates for $SVNAME"
                fi
            done
        else
            echo "[SV-STAT:svstat.sh] $SVNAME exceeded coverage threshold. Exiting..."
    fi
    cd $ORIGDIR
    # reduce disk usage by consolidating/removing intermediate files
    # TODO allow debugging mode to keep these upon request
    rm $SCRDIR/$SVNAME/*.txt*
    cat $SCRDIR/$SVNAME/candidates/*/raw_cj_supp.txt > $SCRDIR/$SVNAME/raw_cj_supp.txt
    cat $SCRDIR/$SVNAME/candidates/*/cj_report.txt > $SCRDIR/$SVNAME/cj_report.txt
    cat $SCRDIR/$SVNAME/candidates/*/derAB.fasta > $SCRDIR/$SVNAME/derAB.fasta
    rm -rf $SCRDIR/$SVNAME/fastq $SCRDIR/$SVNAME/bam $SCRDIR/$SVNAME/candidates
    # Make a note of the name of the last SV in the job. We'll use it in the name of a log file
    LASTSVNAME=$SVNAME
    # Keep track of the SVs processed in this job
    echo "$SVNAME" >> $MYTMPFILE
done < $MYMODELS

# summarize results and clean up
MYRUN="${FIRSTSVNAME}.${LASTSVNAME}"
cat $SCRDIR/*/cj_report.txt > ${MYRUN}.cjrep.txt
cat $SCRDIR/*/raw_cj_supp.txt > ${MYRUN}.rawcjsupp.txt
cat $SCRDIR/*/derAB.fasta > ${MYRUN}.derAB.fasta
while read SVNAME; do rm -rf $SCRDIR/$SVNAME; done < $MYTMPFILE
cp -p $MYTMPFILE ${MYRUN}.svlist.txt
rm $MYTMPFILE
# attach no. of supporting reads to each junction. note: a junction appearing multiple times in the run will have artificially inflated number of supporting reads
MYALN="BWA"
join -j1 <(sort -u -k1,1 ${MYRUN}.cjrep.txt) <(cut -d" " -f1 ${MYRUN}.rawcjsupp.txt | sort | uniq -c | awk '{print $2,$1;}') > ${MYRUN}_${MYALN}_tmpreport.txt
# attach no. of times junction appears in run, and normalize the number of supporting reads by the number of times the junction appears in the run
join -j1 <(sort -k1,1 ${MYRUN}_${MYALN}_tmpreport.txt) <(cut -d" " -f1 ${MYRUN}.cjrep.txt | sort | uniq -c | awk '{print $2,$1;}') | awk '{nsupp=$6; nhit=$7; print $0,nsupp/nhit;}' > ${MYRUN}_${MYALN}_out.txt
# remove redundant files
rm ${MYRUN}_${MYALN}_tmpreport.txt
rm ${MYRUN}.cjrep.txt

# report SVs in vcf-like format upon request
if [ "$META" == "True" ]
    then
        MYVCFOUT="${MYRUN}.svp"
        # print headers
        metagen > $MYVCFOUT
        # process each reported SV
        cat ${MYRUN}_${MYALN}_out.txt | tr ":" " " | tr "_" " " | tr " " "\t" |
        awk ' \
          BEGIN{i=0; OFS="\t";} \
          { \
            fchr=$1; fcoord=$2; fori=$3; schr=$4; scoord=$5; sori=$6; logsupp=$10; normnumreads=$13; \
            start=fcoord-1; end=scoord; svpt="DEL"; svt="DEL"; svx="."; svo1="+"; svo2="+"; svf="FAIL"; \
            svblsc="."; svblnr="."; tsize=scoord-fcoord; svs=sqrt(tsize*tsize); printme="True"; \
            if (fchr==schr) \
              { \
                if (fori==sori) \
                  { \
                    if (fori=="+") \
                      { if (fcoord>scoord) {svt="DUP"; svpt="INS"; start=scoord-1; end=scoord;} } \
                    else \
                      { if (fcoord<scoord) {svt="DUP"; svpt="INS"; end=fcoord;} } \
                  } \
                else \
                  { \
                    svt="INV"; svpt="MIS"; \
                    if (fcoord>scoord) { start=scoord-1; end=fcoord; svo1=sori; svo2=fori; } \
                  }; \
              } \
            else \
              { \
                svt="CTX"; svpt="MIS"; svo1=fori; svo2=sori; svx=schr":"scoord; svs="."; end=fcoord;\
              } \
            \
            if ( '$SLOPE'*logsupp+'$INTERCEPT'-log(normnumreads)/log(10)>=0 && normnumreads<='$MAXSUPPREADS') \
              { \
                svf="PASS"; \
              } \
            else \
              { \
                if ("'$FILT'"=="True") printme="False"; \
              } \
            \
            if (printme=="True") \
              { \
                i++; \
                thisinfo = "SVPNAME=VT"i; \
                thisinfo = thisinfo ";VCT="svpt; \
                thisinfo = thisinfo ";VTO1="svo1; \
                thisinfo = thisinfo ";VTO2="svo2; \
                thisinfo = thisinfo ";VTBLSC="svblsc; \
                thisinfo = thisinfo ";VTBLNR="svblnr; \
                thisinfo = thisinfo ";VTBWSC="logsupp; \
                thisinfo = thisinfo ";VTBWNR="normnumreads; \
                thisinfo = thisinfo ";VTTY="svt; \
                thisinfo = thisinfo ";VTXP="svx; \
                thisinfo = thisinfo ";VTF="svf; \
                print fchr,".",start,".",".",end,".",svpt,svs,".",".",thisinfo; \
              } \
          } \
        ' >> $MYVCFOUT # 0-based, half-open coordinates, e.g. substring[2,5) of "ATCGCTAT" returns "CGC"
fi