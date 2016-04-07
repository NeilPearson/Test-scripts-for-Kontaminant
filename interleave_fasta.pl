#!/usr/bin/perl -w
use strict;

# This deals with supplying paired-end data to BLASTn, in a way that I can later screen it in a more informative manner.
#ÊI may as well wrap fastq->fasta up in here as well, to cut down the number of commands a bit too.

my $r1 = $ARGV[0];
my $r2 = $ARGV[1];

open R1, "<", $r1 or die "ERROR: Could not open fastq file $r1: $!\n";
open R2, "<", $r2 or die "ERROR: Could not open fastq file $r2: $!\n";

my @buffer = ();
my @r1 = ();
while (my $line = <R1>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        my $seq = $buffer[1];
        $header =~ s/^\@/\>/;
        push @r1, "$header\n$seq\n"; 
        @buffer = ();
    }
}
close R1;

@buffer = ();
my @r2 = ();
while (my $line = <R2>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        my $seq = $buffer[1];
        $header =~ s/^\@/\>/;
        push @r2, "$header\n$seq\n"; 
        @buffer = ();
    }
}
close R2;

foreach my $i (1..@r1) {
    $i --;
    print $r1[$i];
    print $r2[$i];
}
