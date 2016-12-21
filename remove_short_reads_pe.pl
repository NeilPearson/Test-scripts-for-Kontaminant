#!/usr/bin/perl -w
use strict;

# This removes reads that have the 'CONTAMINANT' tag that I added in another script, when introducing artificial contamination.
# Useful after you've done filtering, and want to look only at reads that are in the original dataset.

# In order to make BWA work, this needs to have a paired-end version. This is it.


my $r1 = $ARGV[0];
my $r2 = $ARGV[1];
my $minlength = $ARGV[2];
my $out1 = $ARGV[3];
my $out2 = $ARGV[4];

open R1, "<", $r1 or die "ERROR: Could not open fastq file $r1: $!\n";
open O1, ">", $out1 or die "ERROR: Could not open output fastq $out1: $!\n";

my @buffer = ();
my %r1_buffer = ();
while (my $line = <R1>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        # Get the (unique) first word from the header
        my @tagsplit = split /\s/, $header;
        $tagsplit[0] =~ s/[12]$//g;
        
        my $seq = $buffer[1];
        if (length($seq) >= $minlength) {
            $r1_buffer{$tagsplit[0]} = \@buffer;
        }
        @buffer = ();
    }
}
close R1;

# Read the second file; if the read here is valid, and there's a corresponding read in %r1_buffer, write both. Neat.
open R2, "<", $r2 or die "ERROR: Could not open fastq file $r2: $!\n";
open O1, ">", $out1 or die "ERROR: Could not open output fastq $out1: $!\n";
open O2, ">", $out2 or die "ERROR: Could not open output fastq $out2: $!\n";

@buffer = ();
while (my $line = <R2>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer == 4) {
        my $header = $buffer[0];
        # Get the (unique) first word from the header
        my @tagsplit = split /\s/, $header;
        $tagsplit[0] =~ s/[12]$//g;
        
        my $seq = $buffer[1];
        if ((length($seq) >= $minlength) && ($r1_buffer{$tagsplit[0]})) {
            foreach my $i (@{$r1_buffer{$tagsplit[0]}}) { print O1 "$i\n"; }
            foreach my $i (@buffer) { print O2 "$i\n"; }
        }
        @buffer = ();
    }
}

close R2;
close O1;
close O2;
