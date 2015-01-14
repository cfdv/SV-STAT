#!/bin/bash
# purpose: Bring in all partially aligned reads around top-scoring candidate junctions to refine the support scores
# author: Caleb Forbes Davis V
# date: 11/8/11
# external requires:
#
# internal requires:
#
# assumes:
#
# generates:
#
# versions:

# die if abused
if [[ ( $# -gt 0 ) ]]
    then
	s=$1
    else
        echo "svstat_refine sampleID"
        exit 1
fi

#go to this sample's directory
cd $s/SVSTAT/
#annotate raw candidate junction support files with translocation type and sample
for t in 1_19 12_21 4_11 9_22; do cd $t/_FASTA; for p in */; do cd $p; page=${p%\/*}; awk '{print $0,"'$s'","'$t'","'$p'";}' raw_cj_supp.txt > raw_cj_supp_v2.txt; cd ../; done; cd ../..; done
#compile all supporting reads, and pick just the ones with the top alignment scores
cat */_FASTA/*/raw_cj_supp_v2.txt | sort +3 -4 +7 -8 -nr | sort -k4,4 -u > bwa_best_hit_raw_cj_supp.txt
#summarize support for candidate junctions
awk '{cjID=$1; b=$6; q=$7; bq=$8; jtt=$10; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq; nr[cjID]+=1; tt[cjID]=jtt;}END{for (cjID in bs){if (bqs[cjID]!=0) print tt[cjID], cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10), nr[cjID];};}' bwa_best_hit_raw_cj_supp.txt | awk '{print "'$s'",$0;}' > bwa_best_hit_cj_report.txt
#determine pairwise distances between each candidate junction
awk '{cjID=$3; split(cjID,c,":"); log10_supp=$7; split(c[3],b2,"_"); ori=b2[1]; chrA=c[1]; Acoord=c[2]; chrB=b2[2]; Bcoord=c[4]; cj[NR]=chrA","Acoord","chrB","Bcoord","log10_supp","cjID;} END{for (i=1;i<=length(cj)-1;i++){ for (j=i+1;j<=length(cj);j++){ split(cj[i],ibkpnts,","); split(cj[j],jbkpnts,","); ichrA=ibkpnts[1];iAcoord=ibkpnts[2];ichrB=ibkpnts[3];iBcoord=ibkpnts[4]; ilog10_supp=ibkpnts[5]; icjID=ibkpnts[6]; jchrA=jbkpnts[1];jAcoord=jbkpnts[2];jchrB=jbkpnts[3];jBcoord=jbkpnts[4]; jlog10_supp=jbkpnts[5]; jcjID=jbkpnts[6]; dist=sqrt((iAcoord-jAcoord)^2+(iBcoord-jBcoord)^2); pw_supp=ilog10_supp+jlog10_supp; if(ichrA==jchrA && dist!=0) print cj[i],cj[j],dist,pw_supp;} }}' bwa_best_hit_cj_report.txt > cj_junction_bkpnt_dist_v3.txt
#consolidate any highly supported, similar junctions
cp bwa_best_hit_raw_cj_supp.txt cons_bwa_best_hit_raw_cj_supp.txt; awk '{zsupp_thresh=2.27; d_thresh=20; med=2.79239; mad=.7706258; supp=$4; d=$3; zsupp=(supp-med)/mad; if (zsupp>=zsupp_thresh && d<=d_thresh) print $0;}' cj_junction_bkpnt_dist_v3.txt | sort +3 -4 -n | awk '{split($1,icj,","); split($2,jcj,","); min_cjID=icj[6]; maj_cjID=jcj[6]; iscore=icj[5]; jscore=jcj[5]; if (iscore>jscore) {min_cjID=jcj[6]; maj_cjID=icj[6];} print min_cjID,maj_cjID;}' | while IFS=' ' read -ra CJS; do sed -i 's/'${CJS[0]}'/'${CJS[1]}'/g' cons_bwa_best_hit_raw_cj_supp.txt; done
#tally support across all reads for consolidated junctions
awk '{cjID=$1; b=$6; q=$7; bq=$8; jtt=$10; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq; nr[cjID]+=1; tt[cjID]=jtt;}END{for (cjID in bs){if (bqs[cjID]!=0) print tt[cjID], cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10), nr[cjID];};}' cons_bwa_best_hit_raw_cj_supp.txt | awk '{print "'$s'",$0;}' > cons_bwa_best_hit_cj_report.txt

#refine the support by bringing in all partially aligned reads flanking top-scoring candidate junctions
mkdir finish
cd finish
#grab the fasta for the top-scoring junctions
awk '{zsupp_med=1.44716; zsupp_mad=.5899858; zsupp=($7-zsupp_med)/zsupp_mad; if (zsupp>2.1) print $0,zsupp;}' ../cons_bwa_best_hit_cj_report.txt | sort +6 -7 -nr | cut -d" " -f3 | grep -A1 -f - ../*/_FASTA/*/derAB.fasta | sed -e 's/..\/[0-9]*_[0-9]*\/_FASTA\/[0-9]*\/derAB.fasta.//g' | grep -v "^--$" > best_derAB.fasta
#index the top-scoring junctions for bwa
bwa index -a is best_derAB.fasta
#index with samtools
samtools faidx best_derAB.fasta
#express local coordinates along + strand of query
bamtobl.pl ../4.bam | awk '$7<0{l=$9; qs=$3; qe=$4; $3=l-qe+1; $4=l-qs+1; print $0}$7>=0{print $0;}' > 6.bl
#retain only the partially aligned reads
awk '{if ($3!=1 || $4!=$9) print $0;}' 6.bl > 7.bl

#grab only those reads falling within the regions flanking the best candidate junctions. Project the reads onto the reference so we know what the read would have looked like if it was from the reference
FastaToTbl best_derAB.fasta | awk '{split($1,c,":"); split(c[3],b2,"_"); chrA=c[1]; Acoord=c[2]; Aori=b2[1]; chrB=b2[2]; Bcoord=c[4]; Astart=Acoord-1000; Aend=Acoord; Bstart=Bcoord; Bend=Bcoord+1000; Bori="-"; if (Aori=="-") {Astart=Acoord; Aend=Acoord+1000;} if (c[5]=="-") {Bstart=Bcoord-1000; Bend=Bcoord; Bori="+";} print chrA,Astart,Aend,Aori,$1"\n"chrB,Bstart,Bend,Bori,$1;}' | while IFS=" " read -ra COORDS; do awk '{chr=$2; if ($7>0) ori="+"; else ori="-"; tail5=$3-1; tail3=$9-$4; if (ori=="+") {start=$5-tail5; end=$6+tail3; coord=start;} else {start=$5-tail3; end=$6+tail5; coord=end;} js='${COORDS[1]}'; je='${COORDS[2]}'; if (chr=="'${COORDS[0]}'" && coord>js && coord<je && ori=="'${COORDS[3]}'") print $2,$1,"'${COORDS[4]}'",start,end,ori;}' 7.bl; done | sort -k1,2 -u > reads_close_and_ori_to_breakpoints.txt
#grab the sequence and quality of the reads
paste -d" " reads_close_and_ori_to_breakpoints.txt <(awk '{print $1":"$4"-"$5;}' reads_close_and_ori_to_breakpoints.txt  | samtools faidx ~/goldenpath/hg18/hg18.fasta `cat /dev/stdin` | FastaToTbl | cut -d" " -f2) <(cut -d" " -f2 reads_close_and_ori_to_breakpoints.txt | cdbyank ../pass_qc.fastq.cidx | FastqToTbl | cut -d" " -f2,4) > seq_reads_close_ori_breakpoints.txt
#grab the reference sequence from the same region
paste -d" " seq_reads_close_ori_breakpoints.txt <(awk '{print ">"$1":"$4"-"$5"\n"toupper($7);}' seq_reads_close_ori_breakpoints.txt | fastx_reverse_complement | FastaToTbl | cut -d" " -f2) > rc_seq_reads_close_ori_breakpoints.txt
#accumulate test and reference reads flanking top-scoring junctions
awk '{ref_seq=$7; if ($6=="-") ref_seq=$10; r_t_len_diff=length($8)-length(ref_seq); if (r_t_len_diff > -30 && r_t_len_diff < 30) {if (length($8)<=length(ref_seq)) clip_len=length($8); else clip_len=length(ref_seq); print "@"$2"_"$1"_test\n"$8"\n+\n"$9"\n@"$2"_"$1"_ref\n"substr(ref_seq,1,clip_len)"\n+\n"substr($9,1,clip_len);}}' rc_seq_reads_close_ori_breakpoints.txt > flanking_seq_test_ref.fastq

#align the flanking reads to the top-scoring junctions
bwa bwasw best_derAB.fasta flanking_seq_test_ref.fastq > 0.sam
#index, sort, etc
samtools view -bt derAB.fasta.fai -o 0.bam 0.sam; samtools sort 0.bam 1; samtools index 1.bam
#generate junction summary report
bamToBed -i 1.bam | tee -i 1.bed | awk '{cjID=$1; start=$2; end=$3; mq=$5; lb=501-start; rb=end-500; if (lb<rb) minb=lb; else minb=rb; if (minb < 0) minb = 0; print $1,$2,$3,$4,$5, minb, mq, minb*mq;}' | tee -i raw_cj_supp.txt | awk '{cjID=$1; b=$6; q=$7; bq=$8; bs[cjID]+=b; qs[cjID]+=q; bqs[cjID]+=bq;}END{for (cjID in bs){if (bqs[cjID]!=0) print cjID, bs[cjID], qs[cjID], bqs[cjID], log(bqs[cjID])/log(10);};}' | sort +3 -4 -rn > cj_report.txt
#subtract support from normal reads
awk '{mult=1; cjID=$1; b=$6; q=$7; bq=$8; split($4,readID,"_"); origin=readID[length(readID)]; if (origin=="ref") mult=-1; bs[cjID]+=(b*mult); qs[cjID]+=(q*mult); bqs[cjID]+=(bq*mult);}END{for (cjID in bs) {if (bqs[cjID]<=0) log10bqs=-1; else log10bqs=log(bqs[cjID])/log(10); print cjID, bs[cjID], qs[cjID], bqs[cjID], log10bqs;};}' raw_cj_supp.txt > normalized_cj_report.txt
#go back to the original directory (original_directory/sample/SVSTAT/finish)
cd ../../..
