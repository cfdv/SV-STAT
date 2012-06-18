#!/usr/bin/perl -w
# generic bioperl format converter
# source -> http://www.bioperl.org/wiki/HOWTO:SeqIO

use Bio::SeqIO;
# get command-line arguments, or die with a usage statement
my $usage     = "in2out.pl informat outformat\n";
my $informat  = shift or die $usage;
my $outformat = shift or die $usage;

#### author = Caleb Davis
if ($informat eq "fastq-flowsim" && $outformat eq "fastq"){
        # flowsim's quality scores need to be rescaled
	# create one SeqIO object to read in
	my $seqin = Bio::SeqIO->new(
		                    -fh     => \*STDIN,
		                    -format => 'fastq',
		                    );

        my $outseq = Bio::SeqIO->new(
                             -fh     => \*STDOUT,
                             -format => 'fastq',
                             );
	 
	# write each entry in the input to the output
	while (my $inseq = $seqin->next_seq) {
            # rescale quality scores from 28-60 -> 28-40
            my @phreds = ();
            $str_qual = $inseq->qual_text;
            @quals = split(/ /, $str_qual);
	    foreach my $p (@quals){
                if ($p <= 28) { push(@phreds, $p); }
                else { push ( @phreds, int(28 + 12*(($p-28)/32)) ); }
            }
            #well, easiest to ignore the description for now
            #my $rescaled_inseq = Bio::Seq::Quality->new( -display_id => $inseq->display_id . " " . $inseq->desc,
            my $rescaled_inseq = Bio::Seq::Quality->new( -display_id => $inseq->display_id,
                                                         -seq => $inseq->seq,
                                                         -qual => join(" ", @phreds));
            $outseq->write_seq($rescaled_inseq);
	}
}
#### sincerest apologies to bioperl for this outrageous hack

else{ 
    # create one SeqIO object to read in, and another to write out
    my $seqin = Bio::SeqIO->new(
                            -fh     => \*STDIN,
                            -format => $informat,
                            );
    my $outseq = Bio::SeqIO->new(
                             -fh     => \*STDOUT,
                             -format => $outformat,
                             );
 
    # write each entry in the input to the output
    while (my $inseq = $seqin->next_seq) {
        $outseq->write_seq($inseq);
    }
}
