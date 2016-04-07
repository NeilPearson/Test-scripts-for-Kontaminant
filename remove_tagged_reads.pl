#!/usr/bin/perl -w
use strict;

# This removes reads that have the 'CONTAMINANT' tag that I added in another script, when introducing artificial contamination.
# Useful after you've done filtering, and want to look only at reads that are in the original dataset.

my $r1 = $ARGV[0];

open R1, "<", $r1 or die "ERROR: Could not open fastq file $r1: $!\n";

my @buffer = ();
my @r1 = ();
while (my $line = <R1>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        my $seq = $buffer[1];
        unless ($header =~ / CONTAMINANT /) {
            foreach my $i (@buffer) { print "$i\n"; }
        }
        @buffer = ();
    }
}
close R1;

