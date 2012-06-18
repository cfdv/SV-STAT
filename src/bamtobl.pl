#!/usr/bin/perl
#SAM/BAM-formatted data are the preferred method of packaging NGS alignment results. I need to adapt the stack-n-tail method to analyze this format. The first step is to export the query start, end, and subject start & end coordinates from the bam file. Lincoln Stein has some code in Bio::DB::Sam that does this readily. Here I modify it slightly to also export the readids. Read mapping scores are also obtained.
#This code is modified from the Bio::DB::Sam documentation (http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm). 
#author - Caleb Davis
#date - 2/28/11
#input - bam file location. use full path
#output format - readid chrom qstart qend sstart send ori BWAmapscore
#purpose - generate .bl-like alignment coordinates from bam file.

 use Bio::DB::Sam;

 my $usage = "perl bamtobl.pl path_to_bam";
 my $bampath = shift or die $usage;

 # high level API
 my $sam = Bio::DB::Sam->new(-bam => $bampath);

 my @targets    = $sam->seq_ids;
 foreach $target (@targets) {
    my $length = $sam->length($target);
    my $segment = $sam->segment(-seq_id => $target); 
    my @alignments = $segment->features;
    for my $a (@alignments) {

      # where does the alignment start in the reference sequence
      my $seqid  = $a->display_name;
      my $start  = $a->start;
      my $end    = $a->end;
      my $strand = $a->strand;

      # where does the alignment start in the query sequence
      my $query_start = $a->query->start;     
      my $query_end   = $a->query->end;
      my $query_seq   = $a->query->seq->seq;
      my $query_len   = length($query_seq);

      my $match_qual= $a->qual;       # quality of the match
      print "$seqid $target $query_start $query_end $start $end $strand $match_qual $query_len\n";
    }

 }
