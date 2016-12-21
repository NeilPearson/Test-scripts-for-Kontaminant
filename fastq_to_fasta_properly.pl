#!/usr/bin/perl -w
use strict;

# The quickest, easiest, most no-nonsense way to convert fastq->fasta.

my $f = $ARGV[0];

open INFILE, "<", $f or die "ERROR: Could not open fastq file $f: $!\n";
my @buffer = ();
while (my $line = <INFILE>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        my $seq = $buffer[1];
        $header =~ s/^\@/\>/;
        print "$header\n$seq\n"; 
        @buffer = ();
    }
}
close INFILE;
