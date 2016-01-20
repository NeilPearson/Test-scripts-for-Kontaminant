#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);
use File::Basename;

# Takes a reads file, a contaminant reads file and an integer as input.
# Mixes the supplied number of reads from the contaminant in with the other reads, shuffles them up so they're in a different order, and
# dumps everything out to stdout (from which it can be redirected to a file).
# Process is a tad more complex for paired end reads, but nothing we can't handle.

my $r1_file = $ARGV[0];
my $r2_file = $ARGV[1];
my $contaminant_file = $ARGV[2];
my $use_n_contaminant_reads = $ARGV[3];
my $contaminant = basename($contaminant_file);

open(READS, "<", $r1_file) or die "Cannot open reads file $r1_file: $!\n";
my @r1lines = <READS>;   close READS;

my @r1s = ();
my @buffer = ();
my $maxread = 0;
foreach my $line (@r1lines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        chomp $read;
        push @r1s, $read;
        if (length $buffer[1] > $maxread) { $maxread = length $buffer[1]; }
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
        chomp $read;
        push @r2s, $read;
        if (length $buffer[1] > $maxread) { $maxread = length $buffer[1]; }
        @buffer = ();
    }
}
@r2lines = ();

open(CONTAMINANTS, "<", $contaminant_file) or die "Cannot open contaminant file $contaminant_file: $!\n";
my @contaminantlines = <CONTAMINANTS>;
# Assume contaminant is in FASTA format
# (On the grounds that this program is meant to generate fake reads from a reference genome).
# However, contaminant reads need to be output in FASTQ format, which means adding some extra stuff.
# (Assume quality scores don't matter, so they can be anything).
# Oh right, and I'll need to match the length of the contaminant reads to those in the existing dataset.
# That will mean getting the max length of line in that data (see above).
# (Assume Illumina-type data with a fixed max read length).
# Oh, and one more thing: the contaminant ref data will be quite long, mostly, and in one or a few contigs.
# I'll have to handle that in an appropriate way. 

my @contaminant_contigs = ();
@buffer = ();
foreach my $line (@contaminantlines) {
    if (($line =~ /^>/) && (@buffer)) {
        # FASTA read may be spread over multiple lines. Header will always be in the first line however.
        my $header = shift @buffer;
        my $seq = join '', @buffer;
        my $contig = "$header\n$seq";
        # A sanity check: only store this contig if the seq is longer than the length of fake read we want.
        if (length $seq > $maxread) { push @contaminant_contigs, $contig; }
        @buffer = ();
    }
    chomp $line;
    push @buffer, $line;
}
# Make sure to get the last contig out of the buffer. This is important because there may only be one sequence in the file.
if (@buffer) {
    # FASTA read may be spread over multiple lines. Header will always be in the first line however.
    my $header = shift @buffer;
    my $seq = join '', @buffer;
    my $contig = "$header\n$seq";
    # A sanity check: only store this contig if the seq is longer than the length of fake read we want.
    if (length $seq > $maxread) { push @contaminant_contigs, $contig; }
    @buffer = ();
}
@contaminantlines = ();

# Now, it's possible to generate a number of fake FASTQ reads and put them in with the rest.
my @qscores = qw/A B C D E F/;
for my $i (1..$use_n_contaminant_reads) {
    # Select a contig at random. (Yeah, this will give a bias towards the shorter contigs, if there's more than one. I'll do a weighted random choice in a later version).
    my $contig = $contaminant_contigs[rand @contaminant_contigs];
    # Tag the read as a contaminant so I can identify it later.
    # This tag also needs to be unique, so I'll probably need to add the number stored in $i or something.
    my @sp1 = split /\n/, $contig;
    my $tag = $sp1[0];
    chomp $tag;
    $tag .= " CONTAMINANT $contaminant uid $i";
    $tag =~ s/\>/\@/;
    $sp1[0] = $tag;
    # Get a subset of the sequence with an appropriate length
    my $seq = substr($sp1[1], rand(length($sp1[1]) - $maxread - 2), $maxread - 1);
    $sp1[1] = $seq;
    $sp1[2] = "+";
    my $qstr = ();
    for (1..length($sp1[1])) { $qstr .= $qscores[rand @qscores]; }
    $sp1[3] = $qstr;
    
    my $fakeread = join "\n", @sp1;
    
    my $r1 = $fakeread;
    my $r2 = $fakeread;
    
    # Read 2 needs to be the reverse complement of read 1
    # (Yes, I know that due to the way paired ends work, I should pick a sequence a bit further along the contig too.
    # I'll do that in the future).
    my @sp2 = split /\n/, $r2;
    $seq = $sp2[1];
    chomp $seq;
    $seq = reverse_complement($seq);
    $sp2[1] = $seq;
    $r2 = join "\n", @sp2;
    
    push @r1s, $r1;
    push @r2s, $r2;
}

my @indexes = (1..@r1s);
@indexes = shuffle(@indexes);

# Oh man, I can't just dump output to stdout any more. What a pain.
# I'll have to modify the input filenames and set up those as output.

$contaminant =~ s/.fasta//g;
my $r1out = $r1_file;
my $r2out = $r2_file;
$r1out =~ s/.fastq/_CONTAMINATED_$contaminant.fastq/;
$r2out =~ s/.fastq/_CONTAMINATED_$contaminant.fastq/;
open(R1OUT, ">", $r1out) or die "Cannot open output $r1out: $!\n";
open(R2OUT, ">", $r2out) or die "Cannot open output $r2out: $!\n";
my $c = 0;
foreach my $i (@indexes) {
    $c ++;
    my $r1 = $r1s[$i-1];
    my $r2 = $r2s[$i-1];
    print R1OUT $r1;
    print R2OUT $r2;
    print R1OUT "\n"; print R2OUT "\n";
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