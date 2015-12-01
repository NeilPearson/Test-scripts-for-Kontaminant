#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

# Takes a reads file, a contaminant reads file and an integer as input.
# Mixes the supplied number of reads from the contaminant in with the other reads, shuffles them up so they're in a different order, and
# dumps everything out to stdout (from which it can be redirected to a file).
# Process is a tad more complex for paired end reads, but nothing we can't handle.

my $r1_file = $ARGV[0];
my $r2_file = $ARGV[1];
my $contaminant_file = $ARGV[2];
my $use_n_contaminant_reads = $ARGV[3];

open(READS, "<", $r1_file) or die "Cannot open reads file $r1_file: $!\n";
my @r1lines = <READS>;   close READS;

my @r1s = ();
my @buffer = ();
foreach my $line (@r1lines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        push @r1s, $read;
        @buffer = ();
    }
}
@r1lines = ();

open(READS, "<", $r2_file) or die "Cannot open reads file $r2_file: $!\n";
my @r2lines = <READS>;   close READS;

my @r2s = ();
@buffer = ();
foreach my $line (@r2lines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        push @r2s, $read;
        @buffer = ();
    }
}
@r2lines = ();

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
    my $read = $contaminant_reads[$i-1];
    # Tag the read as a contaminant so I can identify it later.
    my @sp1 = split /\n/, $read;
    my $tag = $sp1[0];
    chomp $tag;
    $tag .= "_CONTAMINANT\n";
    $read = join '', @sp1;
    
    my $r1 = $read;
    my $r2 = $read;
    
    @sp1 = split /\n/, $r2;
    my $seq = $sp1[1];
    chomp $seq;
    $seq = reverse_complement($seq);
    $sp1[1] = "$seq\n";
    $r2 = join '', @sp1;
    
    push @r1s, $r1;
    push @r2s, $r2;
}

my @indexes = (1..@r1s);
@indexes = shuffle(@indexes);

# Oh man, I can't just dump output to stdout any more. What a pain.
# I'll have to modify the input filenames and set up those as output.
my $r1out = $r1_file;
my $r2out = $r2_file;
$r1out =~ s/.fastq/_CONTAMINATED.fastq/;
$r2out =~ s/.fastq/_CONTAMINATED.fastq/;
open(R1OUT, ">", $r1out) or die "Cannot open output $r1out: $!\n";
open(R2OUT, ">", $r2out) or die "Cannot open output $r2out: $!\n";
foreach my $i (@indexes) {
    my $r1 = $r1s[$i-1];
    my $r2 = $r2s[$i-1];
    print R1OUT $r1;
    print R2OUT $r2;
}
close R1OUT;
close R2OUT;

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}