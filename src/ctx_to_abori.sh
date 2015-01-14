#!/bin/bash
#Structural variants (SVs) aberrantly connect genomic regions to each other. Existing tools use different approaches
# to solve problems related to observations of nucleotide sequences shorter than the repetitive regions connected by
# SVs. Assuming the user suspects fusions in specific orientations between specific regions, SV-STAT addresses these
# problems in three steps. 1) Create a synthetic reference of candidate SVs, 2) align misaligned reads to the synthetic
# reference, and 3) score support for SVs provided by what could now be correctly aligned reads.
#author - Caleb Davis
#date - 8/15/13
#input - bam file location, and candidate SVs ( in breakdancer's ctx format)
#output format - SVNAME FIRSTCHR FIRSTSTART FIRSTEND FIRSTORI SECONDCHR SECONDSTART SECONDEND SECONDORI RECIPROCAL
#purpose - Breakdancer's ctx format tells you which regions are connected, and roughly how many reads on either side of the 
# breakpoint are mis-mated and support the junction. However, the format does not explicitly specify the order and orientations
# of the regions in the SV. Given a candidate SV in breakdancer format, this code will reach into the bam and determine
# the order and orientations of regions for all junction(s) supported by reads flanking the breakpoints
# external requires:
# -samtools

# die if abused
if [[ ( -e $1 ) && ( -e $2 ) ]]
    then
        MYBAM=$1
        MYCTX=$2
    else
        echo "ctx_to_abori aln.bam candidate_SVs.txt"
        exit 1
fi

#process candidate SVs in breakdancer format line by line
#set window size around breakpoint around which to search for discordant pairs
BUFFERBP=1000
#set minimum distance between two breakpoints on same chromosome
SAMECHRMINDIST=50
awk '$1!~/GL000/ && $4!~/GL000/ && $1!~/MT/ && $4!~/MT/' $MYCTX | #get rid of SVs in GL000/MT
awk '$1==$4{diff=$5-$2; if (diff >= '$SAMECHRMINDIST') print $0;}$1!=$4{print $0;}' | #do not consider small SVs/indels
cut -f1-2,4-5 | while read FIRSTCHR FIRSTCOORD SECONDCHR SECONDCOORD; do
  FIRSTSTART=`echo -e "$(( $FIRSTCOORD-$BUFFERBP ))"`; FIRSTEND=`echo -e "$(( $FIRSTCOORD+$BUFFERBP ))"`
  SECONDSTART=`echo -e "$(( $SECONDCOORD-$BUFFERBP ))"`; SECONDEND=`echo -e "$(( $SECONDCOORD+$BUFFERBP ))"`
  FIRSTRANGE=`echo -e "$FIRSTCHR:$FIRSTSTART-$FIRSTEND"`; SECONDRANGE=`echo -e "$SECONDCHR:$SECONDSTART-$SECONDEND"`
  samtools view -F14 $MYBAM $FIRSTRANGE $SECONDRANGE | 
  sort -u -k1,1 |
  awk 'BEGIN{OFS="\t";}{if ($7=="=")$7="'$FIRSTCHR'"; print $0;}' | #explicitly annotate chromosome of mate
  awk '\
    BEGIN{\
      f['$FIRSTCHR',"start"]='$FIRSTSTART'; f['$FIRSTCHR',"end"]='$FIRSTEND'; \
      s['$SECONDCHR',"start"]='$SECONDSTART'; s['$SECONDCHR',"end"]='$SECONDEND'; \
    } \
    ($4>f[$3,"start"] && $4<f[$3,"end"] && $8>s[$7,"start"] && $8<s[$7,"end"]) || \
      ($4>s[$3,"start"] && $4<s[$3,"end"] && $8>f[$7,"start"] && $8<f[$7,"end"])\
  ' |                         # only passes reads along if they fall within a 2kb window centered at the breakpoints
  cut -f2-3,7 |
  while read DECFLAG READCHR MATECHR; do
    BINFLAG=`echo -e "obase=2;$DECFLAG" | bc | rev`
    echo $DECFLAG $BINFLAG $READCHR $MATECHR
  done |                      # flags converted to binary
  awk '\
    BEGIN{\
      a['$FIRSTCHR',"+",'$SECONDCHR',"+"]=0; a['$FIRSTCHR',"+",'$SECONDCHR',"-"]=0; \
      a['$FIRSTCHR',"-",'$SECONDCHR',"+"]=0; a['$FIRSTCHR',"-",'$SECONDCHR',"-"]=0; \
    } \
    { \
      readori="+"; mateori="+"; \
      if (substr($2,5,1)==1) readori="-"; \
      if (substr($2,6,1)==1) mateori="-"; \
      a[$3,readori,$4,mateori]++; \
    } \
    END{ \
      print a['$FIRSTCHR',"+",'$SECONDCHR',"+"] "\t" a['$FIRSTCHR',"+",'$SECONDCHR',"-"] "\t" \
            a['$FIRSTCHR',"-",'$SECONDCHR',"+"] "\t" a['$FIRSTCHR',"-",'$SECONDCHR',"-"]; \
    } \
  ' |                         # tally of reads supporting the four possible orientations 
  awk '\
    { \
      for (i=1; i<=4; i++) { \
        a[i,"r"]=0; \
        if ($i>0) a[i,"w"]=1; else a[i,"w"]=0; \
      } \
      if ($1>0 && $4>0) { \
        a[1,"r"]=1; a[4,"w"]=0; \
      } \
      if ($2>0 && $3>0) { \
        a[2,"r"]=1; a[3,"w"]=0; \
      } \
      for (i=1; i<=4; i++) { \
        if (i<=2) firstori="+"; else firstori="-"; \
        if ((i%2)==0) secondori="+"; else secondori="-"; \
        svname="'$FIRSTCHR'_'$FIRSTCOORD'_" firstori "_" "'$SECONDCHR'_'$SECONDCOORD'_" secondori; \
        if (a[i,"w"]==1) \
          print svname "\t" \
                "chr'$FIRSTCHR'" "\t" "'$FIRSTSTART'" "\t" "'$FIRSTEND'" "\t" firstori "\t" \
                "chr'$SECONDCHR'" "\t" "'$SECONDSTART'" "\t" "'$SECONDEND'" "\t" secondori "\t" \
                a[i,"r"]; \
      } \
    } \
  '                           # switches: "r" for reciprocal, and "w" for write; four orientations
done
