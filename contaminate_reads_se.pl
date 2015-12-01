#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

# Takes a reads file, a contaminant reads file and an integer as input.
# Mixes the supplied number of reads from the contaminant in with the other reads, shuffles them up so they're in a different order, and
# dumps everything out to stdout (from which it can be redirected to a file).

my $reads_file = $ARGV[0];
my $contaminant_file = $ARGV[1];
my $use_n_contaminant_reads = $ARGV[2];

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

# OK, pick some contaminant reads and put them in with the rest.
for my $i (1..$use_n_contaminant_reads) {
    # Tag the read as a contaminant so I can identify it later.
    my $read = $contaminant_reads[$i-1];
    my @sp1 = split /\n/, $read;
    my $tag = $sp1[0];
    chomp $tag;
    $tag .= "_CONTAMINANT\n";
    $read = join '', @sp1;
    
    push @reads, $read;
}

@reads = shuffle(@reads);

foreach my $read (@reads) {
    print $read;
}