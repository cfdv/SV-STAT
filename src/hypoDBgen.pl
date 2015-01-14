#!/usr/bin/perl
#this code converts stacked reads and a reference sequence into a fasta file containing all the possible combinations of stacks. How the stacks are joined together depends on the model you want to test. You could be trying to find deletions, duplications, inversions, translocations, etc. Ideally, it would let me specify whether I wanted derA (for A=chr4, B=chr11, this would be der4. The code expects qarm-qarm junctions), derB, invA (two forward stacks or two reverse stacks from the same chromosomes are joined), invB, etc. or some combination, but for now I'll write it to generate a single fasta file with derA, derB, invA, and invB where:
#-derA = chrA:fstack-chrB:rstack
#-derB = chrB:fstack-chrA:rstack
#-invA = chrA:fstack-(chrA:fstack)	#(breakpoint) means reverse complement before joining with other breakpoint's to form the candidate fusion
#-invB = chrB:fstack-(chrB:fstack)
#v8 - removed path dependencies, fixed paths, starting to elaborate junction connection types
#v9 - implements output paging to avoid hitting the 2^32 limit in BWA indices
#v10 - removes dependence on /_FASTA directory structure.

#use File::Copy;
use Data::Dumper;
use IO::String;
use Bio::PrimarySeq;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub local_fetchrefdna{
	#this code retrieves a range of dna given assembly chromosome, start, and stop
	#assumes .ol format, where fasta header is stripped, and entire sequence is stored on the first line of the file.
	#usage: my $seq = local_fetchrefdna(assembly,chr,start,stop)
	my $assembly = shift;  #path to reference genome; e.g. /path/to/hg19 with $assembly/chromosomes/$chr.ol
	my $chr = shift;
	my $start = shift;
	my $stop = shift;

	#define path to the chromosome files in the assembly directory
	#TODO: use samtools instead of .ol
	my $filename = $assembly . "/chromosomes/$chr" . ".ol";         #e.g. /path/to/hg19/chromosomes/chr22.ol
	open(SEQ, "<$filename") || die("Cannot open $filename");
	
	my $seqout = '';
	my $length = $stop-$start+1;
	my $offset = $start - 1;

	#print "length = $length; offset = $offset\n";
	seek SEQ,$offset,0;	#move pointer to start coordinate
	read SEQ,$seqout,$length;	#extract sequence

	return $seqout;
	close(SEQ);
}

sub make_junction{
	#input: two reference sequence region hashes, A and B. The hashes of the form $stacks{$chr}{$stacktype}{$coord} containing "seq", "sum_tail_lengths", and "has_significant_tail" keys.
	#output: junction hash. It contains an identifier ['Achr':'(1-based) Abreakpoint coordinate':'Aori'_'Bchr':'(1-based) Bbreakpoint coordinate':'Bori'] and sequence.
	#	regardless of orientation, the start coordinate is always less than or equal to the stop coordinate for A,B,and junction hashes. Bchr and Achr can be the same or different.
	#usage: my $jun = make_junction($outfile,\%Astack, $Achr, $Acoord, $Aori, \%Bstack, $Bchr, $Bcoord, $Bori)
	my ($outfile, $Astack, $Achr, $Acoord, $Aori, $Bstack, $Bchr, $Bcoord, $Bori) = @_;
	
	#perform qc on the junction and print
	if (junction_pass_qc($Astack,$Bstack)) {
		my $Aid = $Achr . ":" . $Acoord . ":" . $Aori;
		my $Bid = $Bchr . ":" . $Bcoord . ":" . $Bori;
		my $A_seqobj = Bio::PrimarySeq->new ( -seq => $Astack->{"seq"},
																		 -id  => $Aid);
		my $B_seqobj = Bio::PrimarySeq->new ( -seq => $Bstack->{"seq"},
																		 -id  => $Bid);
		$A_seqobj = $A_seqobj->revcom unless ($Aori eq "+");
		$B_seqobj = $B_seqobj->revcom unless ($Bori eq "+");

		$num_cj++;
		if (($num_cj % $cj_per_page) == 1){
			close(OUT) unless ($current_page==0);
			$current_page++;
			$page_file_name = "$current_page/$outfile";
			if (! -e "$current_page") {
				mkdir("$current_page/") or die "Can't create $current_page:$!\n";
			}
			open(OUT, ">$page_file_name") || die("Cannot open $page_file_name");
		}
		
		print OUT ">" . $Aid . "_" . $Bid . "\n" . $A_seqobj->seq . $B_seqobj->seq . "\n";
	}
}

sub make_translocations{
	#input: two hashes - one per chromosome - of the form $stacks{$chr}{$stacktype}{$coord} containing all forward and reverse stacks on two different chromosomes
	#output: junction hash with all candidate junctions of derA. It contains an identifier ['Achr':'(1-based) Abreakpoint coordinate'-'Aori'_'Bchr':'(1-based) Bbreakpoint coordinate'-'Bori'] and sequence for each candidate translocation junction.
	#	regardless of orientation, the start coordinate is always less than or equal to the stop coordinate for A,B,and junction hashes. 
	#usage: my $juns = make_translocations($outfile,\%Astacks,\%Bstacks)
	my ($outfile,$first_stacks,$second_stacks,$first_ori,$second_ori,$first_type,$second_type) = @_;
	
	#which chromosomes are we dealing with? A and Bstacks should have one chromosome represented, each, and they should not be the same chromosome
	my @first_chr = keys(%{$first_stacks});
	my @second_chr = keys(%{$second_stacks});
  #TODO: fix the name of this function to make_SV
	#( ( (@second_chr == 1) && (@first_chr == 1) ) && ( ($first_chr[0] ne $second_chr[0]) ) ) || die("make_translocations: There are too many chromosomes represented in your stacks, or you are trying to create a translocation between two regions in the same chromosome.\n");

	#configure junction
	my $first_reg = {stacks=>$first_stacks, chrm=>$first_chr[0], type=>$first_type, ori=>$first_ori};
	my $second_reg = {stacks=>$second_stacks, chrm=>$second_chr[0], type=>$second_type, ori=>$second_ori};
	my %jun = ("A", $first_reg,
		   "B", $second_reg);

	#my $Astacktype = $Atype->[0];	#"start" and "end" refer to which side of the stack the tail is on. Therefore "end" stacks are forward stacks, and "start" stacks are reverse stacks. These terms are used interchangeably
	my $Astacks = $jun{"A"}->{"stacks"};
	my $Bstacks = $jun{"B"}->{"stacks"};
	my $Achrm = $jun{"A"}->{"chrm"};
	my $Bchrm = $jun{"B"}->{"chrm"};
	my $Astacktype = $jun{"A"}->{"type"};
	my $Aori = $jun{"A"}->{"ori"};	#"+" and "-" refer to which strand the stack should be read from in order to form the junction.
	my $Bstacktype = $jun{"B"}->{"type"};
	my $Bori = $jun{"B"}->{"ori"};
	foreach my $Acoord (sort { $a <=> $b } keys(%{$Astacks->{$Achrm}{$Astacktype}}) ) {	#now iterate through the A chromosome stacks, using the appropriate stack type
		foreach my $Bcoord (sort { $a <=> $b } keys(%{$Bstacks->{$Bchrm}{$Bstacktype}}) )	 { #now iterate through the B chromosome stacks, using the appropriate stack type
			make_junction($outfile,$Astacks->{$Achrm}{$Astacktype}{$Acoord},$Achrm,$Acoord,$Aori,$Bstacks->{$Bchrm}{$Bstacktype}{$Bcoord},$Bchrm,$Bcoord,$Bori);
		}
	}
}

sub make_inversions{
  #TODO: allow all intrachromosomal SV -- consider deleting this function if "make_translocations" also works for intrachromosomal SV
	#input: one hash of the form $stacks{$chr}{$stacktype}{$coord} with "seq", "sum_tail_lengths", "has_significant_tail", and "reads" keys
	#output: hash of junctions with all candidate junctions for invA. It contains an identifier ['Achr':'(1-based) Abreakpoint1 coordinate'-'Aori'_'Achr':'(1-based) Abreakpoint2 coordinate'-'Aori'] and sequence for each candidate inversion junction.
	#usage: my $juns = make_inversions($outfile,\%Astacks)
	my ($outfile,$Astacks) = @_;
	my $Bstacks = $Astacks;

	my (@Achr) = keys(%{$Astacks});
	(@Achr == 1) || die("make_inversions: There are too many chromosomes represented in your stacks. An inversion can only involve one chromosome.\n");
	my @Bchr = @Achr; #for an inversion, chrA and chrB are the same

	#There are only two ways to make inversions given a set of stacks from a chromosome. Either join the forward or the reverse stacks together. Pick one stack and flip it. Iterate 1:n-1 outer, i+1:n inner
	my @invtypes = (["end",["+","-"]],["start",["-","+"]]);	#join forward stacks AF+_AF-, then reverse stacks AR-_AR+
	foreach my $invtype (@invtypes){
		#this inner loop is called twice, once for each inversion type. 
		my $Astacktype = $invtype->[0];	#forward stacks the first time through, and reverse stacks in the second iteration
		my $Bstacktype = $Astacktype;
		my $Aori = $invtype->[1][0];
		my $Bori = $invtype->[1][1];
		my @stackcoords =();
		if ($Astacktype eq "end"){
			@stackcoords = ( sort { $a <=> $b } keys(%{$Astacks->{$Achr[0]}{$Astacktype}}) );	#sort forward stacks by coordinate, ascending
		} else {
			@stackcoords = ( sort { $b <=> $a } keys(%{$Astacks->{$Achr[0]}{$Astacktype}}) );	#sort reverse stacks by coordinate, descending
		}
		for (my $i=0; $i<= ($#stackcoords - 1); $i++){	#Go through stack combinations pairwise. Order is not important because AF1+_AF2- === AF2+_AF1-
			for (my $j=$i+1; $j<=$#stackcoords; $j++){	#$i and $j cooperate to iterate across a triangular portion of the pairwise matrix excluding the diagonal ($i=$j would be a duplication)
				my $Acoord=$stackcoords[$i];
				my $Bcoord=$stackcoords[$j];
				#construct the junction using the content, chr, coord, and orientation of both stacks, make_junction performs qc and prints
				make_junction($outfile,$Astacks->{$Achr[0]}{$Astacktype}{$Acoord},$Achr[0],$Acoord,$Aori,$Bstacks->{$Bchr[0]}{$Bstacktype}{$Bcoord},$Bchr[0],$Bcoord,$Bori);
			}
		}
	}
}

sub junction_pass_qc{
	#input: two stacks forming a junction together. Expecting a hash with "sum_tail_lengths", and "has_significant_tail" keys
	#output: 1 for pass, 0 for fail
	#usage: junction_pass_qc($Astack,$Bstack)
	my ($Astack, $Bstack) = @_;
	my $junction_qc = 0;
	
	#equals 1 when there are enough tail bases in the candidate juction, and at least one tail in either stack is long enough
	if (($Astack->{"sum_tail_lengths"} + $Bstack->{"sum_tail_lengths"} > $significant_stack_sum_tail_lengths) && 
							($Astack->{"has_significant_tail"} || $Bstack->{"has_significant_tail"})){
		$junction_qc = 1;	
	}
	return $junction_qc;
}

sub flip_stacktype{
	#input: type of stack, either "end" or "start". These correspond to breakpoints with tails (unaligned portion of reads)
  # pointing down- or up-stream, respectively
	#output: "end" if input is "start", or "start" if input is "end"
	#usage: flip_stacktype($my_stacktype)
  #TODO: error check
	my ($input_stacktype) = @_;
  if ($input_stacktype eq "start") {return "end";} elsif ($input_stacktype eq "end") {return "start";}
}

sub flip_ori{
	#input: type of orientation, either "+" or "-". These correspond to the forward or reverse strand in the reference
	#output: "+" if input is "-", or "-" if input is "+"
	#usage: flip_ori($my_ori)
  #TODO: error check
	my ($input_ori) = @_;
  if ($input_ori eq "+") {return "-";} elsif ($input_ori eq "-") {return "+";}
}

#####################           MAIN              ###################################

my $usage = "perl hypoDBgen.pl assembly firstchr firstori secondchr secondori reciprocal stackedreadsfile outfile\n";
my $assembly= shift or die $usage;
my $firstchr = shift or die $usage;
my $firstori = shift or die $usage;
my $secondchr = shift or die $usage;
my $secondori = shift or die $usage;
#sometimes $reciprocal = 0, in which case "$reciprocal = shift" evaluates as FALSE...
my $reciprocal = shift; die $usage unless (defined $reciprocal && $reciprocal ne '');
my $stackedreadsfile = shift or die $usage;
my $outfile = shift or die $usage;	#fasta file of hypothetically rearranged sequences

#define file IO
open(STKRDS, "<$stackedreadsfile") || die("Cannot open $stackedreadsfile");

#algorithm parameters to limit the number of junctions considered
$significant_read_tail_length = 4;	#a junction is not considered if neither one has a read of this length or greater
$significant_stack_sum_tail_lengths = 8;	#a junction is not considered if the sum of the tail lengths in both stacks is this number or less
$significant_deletion_length = 2;	#forward and reverse stack coordinates have to be separated by a gap of this size or larger

#how many candidate junctions per file?
$cj_per_page = 3000000;
$num_cj = 0;
$current_page = 0;

#import the stacked reads file.  form a hash
#the code expects the following column order:
#QUERYID	SUBJECTID	PERC_IDENT	ALIGN_LENGTH	MISMATCHES	GAPOPENINGS	QLENGTH	QSTART	QEND	SSTART	SEND	ORI	EVALUE	BITSCORE	BESTHIT	STACKTYPE
#generate a hash of the form: stacktype->coord->{refseq,readid->seq}
#my $j = 0;
my %stkrds = ();
while (my $line = <STKRDS>){
	@cols = map { trim($_) } split(/ /, $line);	#read in column headers
	my $start = 0;
	my $stop = 0;
	my $chr = $cols[1];
	my $stacktype = $cols[15];
	#my $readid = $cols[0];	
	my @readannots = ();
	push(@readannots,($cols[0],$cols[17],$cols[6]));
	for (my $k=18; $k<=23; $k++){
		push(@readannots,$cols[$k]);
	}
	my $readid = join("_",@readannots);
	my $tail_length = 0;
	if ($stacktype eq 'start'){
		$coord = $cols[9];	#sstart
		$start = $coord;	
		$stop = $start + 500;
		$tail_length = $cols[7] - 1;	#qstart - 1
	} else {
		$coord = $cols[10];	#send
		$start = $coord - 500;
		$stop = $coord;	
		$tail_length = $cols[6] - $cols[8];	#qlength - qend
	}
	#if there is not already sequence associated with the stack then go get it
	if (!exists $stkrds{$chr}{$stacktype}{$coord}{"seq"}){
		$stkrds{$chr}{$stacktype}{$coord}{"seq"} = local_fetchrefdna($assembly, $chr, $start, $stop);
		$stkrds{$chr}{$stacktype}{$coord}{"reads"} = ();
		$stkrds{$chr}{$stacktype}{$coord}{"sum_tail_lengths"} = 0;
		$stkrds{$chr}{$stacktype}{$coord}{"has_significant_tail"} = 0;
		#print "stack sequence retrieved for $stacktype:$coord . . . ";
	}
	push @{$stkrds{$chr}{$stacktype}{$coord}{"reads"}}, $readid;
	$stkrds{$chr}{$stacktype}{$coord}{"sum_tail_lengths"} += $tail_length;
	$stkrds{$chr}{$stacktype}{$coord}{"has_significant_tail"} = 1 unless $tail_length < $significant_read_tail_length;
	#print "$readid added to $stacktype:$coord\n";
	
	#$j++;
	#last if $j == 20;
}

#create two stacks, one per partner chromosome
my %firststacks = ();
my %secondstacks = ();
$firststacks{$firstchr} = $stkrds{$firstchr};
$secondstacks{$secondchr} = $stkrds{$secondchr};

#now collect candidate junctions according to the junction type include toggle settings (I-IV)
#TODO: change names of functions "make_translocations" and "make_inversions" to
#TODO(cont.):  indicate either inter- and intra-chromosomal SV. Refactoring a bit here.
#if ($firstchr ne $secondchr){
  if ($firstori ne $secondori){
    my $stacktype = "end";                                                                                   # t(A;B)(q;p)
    $stacktype = "start" unless ($firstori eq "+");
    make_translocations($outfile,\%firststacks,\%secondstacks,$firstori,$secondori,$stacktype,$stacktype);                          # derA
    if ($reciprocal) {
      $stacktype = flip_stacktype($stacktype);
      $firstori = flip_ori($firstori);
      $secondori = flip_ori($secondori);
      make_translocations($outfile,\%firststacks,\%secondstacks,$firstori,$secondori,$stacktype,$stacktype); # derB
    }
  }
  else {                                                                                     # t(A;B)(q;q)
    make_translocations($outfile,\%firststacks,\%secondstacks,"+","+","end","start") if ($firstori eq "+");                        # derA
    make_translocations($outfile,\%secondstacks,\%firststacks,"+","+","end","start") if ($reciprocal || ($firstori eq "-"));   # derB
  }
#}
#else {
  #TODO: allow intra-chromosomal SV
  #make_inversions($outfile,\%Astacks) if ($inc_I);
  #make_inversions($outfile,\%Bstacks) if ($inc_II);
  #make_inversions($outfile,\%Astacks) if ($inc_III);
  #make_inversions($outfile,\%Bstacks) if ($inc_IV);
#}


#possibly should report the number of stacks, number of stacked reads, number of hypothetical rearrangements

close(STKRDS);
