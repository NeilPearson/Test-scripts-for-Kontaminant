#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

# Takes a reads file, a contaminant reads file and an integer as input.
# Makes the supplied number of reads chimeric, shuffles them up so they're in a different order, and
# dumps everything out to stdout (from which it can be redirected to a file).
# Should also specify a minimum length of chimeric sequence I suppose.

my $reads_file = $ARGV[0];
my $contaminant_file = $ARGV[1];
my $use_n_contaminant_reads = $ARGV[2];
my $min_chimeric_length = $ARGV[3];
if (!$min_chimeric_length) {
    $min_chimeric_length = 20;
}


open(READS, "<", $reads_file) or die "Cannot open reads file $reads_file: $!\n";
my @readslines = <READS>;   close READS;

my @reads = ();
my @buffer = ();
foreach my $line (@readslines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        push @reads, $read;
        @buffer = ();
    }
}
@readslines = ();

open(CONTAMINANTS, "<", $contaminant_file) or die "Cannot open contaminant file $contaminant_file: $!\n";
my @contaminantlines = <CONTAMINANTS>;

my @contaminant_reads = ();
@buffer = ();
foreach my $line (@contaminantlines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        push @contaminant_reads, $read;
        @buffer = ();
    }
}
@contaminantlines = ();

# When simply contaminating a dataset, we can just push in a set number of contaminant reads, but when making a chimeric set it's a bit more complex (albeit not that much more).
# I'll be making heavy use of splice here. Modify the first n reads - they'll be shuffled shortly.
# So we get the basic read and pick a contaminant read at random. We get the sequences out of both, and split them up into an array.
# Then use splice to replace the character scores.
# I could calculate the right numbers to make sure the read stays the same, but it's easier and equally effective to just shorten it to the right length before re-joining it.
for my $i (1..$use_n_contaminant_reads) {
    # push @reads, $contaminant_reads[$i-1];
    # Get the read to modify
    my $read = $reads[$i-1];
    my $contaminant_read = $contaminant_reads[rand @contaminant_reads];
    
    my @sp1 = split /\n/, $read;    my @sp2 = split /\n/, $contaminant_read;
    my $seq = $sp1[1];  my $contaminant_seq = $sp2[1];
    chomp $seq; chomp $contaminant_seq;
    my $seqlen = length $seq;
    my @seqsp = split '', $seq;
    my @consp = split '', $contaminant_seq;
    
    @seqsp = splice @seqsp, int(rand($seqlen-$min_chimeric_length)), 0, @consp;
    @seqsp = splice @seqsp, 0, $seqlen;
    $seq = join '', @seqsp;
    $sp1[1] = "$seq\n";
    
    # Oh yeah, and just to be sure, I probably want to tag the thing as chimeric so I can definitively identify it later.
    my $tag = $sp1[0];
    chomp $tag;
    $tag .= "_CHIMERIC\n";
    $sp1[0] = $tag;
    
    $read = join '', @sp1;
    $reads[$i-1] = $read;
}

@reads = shuffle(@reads);

foreach my $read (@reads) {
    print $read;
}
