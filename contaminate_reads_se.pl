#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

# Takes a reads file, a contaminant reads file and an integer as input.
# Mixes the supplied number of reads from the contaminant in with the other reads, shuffles them up so they're in a different order, and
# prints the resultant dataset to a new file.

my $reads_file = $ARGV[0];
my $contaminant_file = $ARGV[1];
my $use_n_contaminant_reads = $ARGV[2];
my $contaminant = basename($contaminant_file);

open(READS, "<", $reads_file) or die "Cannot open reads file $reads_file: $!\n";
my @readslines = <READS>;   close READS;

my @reads = ();
my @buffer = ();
my $maxread = 0;
foreach my $line (@readslines) {
    push @buffer, $line;
    if (@buffer >= 4) {
        my $read = join '', @buffer;
        push @reads, $read;
        if (length $buffer[1] > $maxread) { $maxread = length $buffer[1]; }
        @buffer = ();
    }
}
@readslines = ();

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
    my $seq = substr($sp1[1], rand(length($sp1[1]) - $maxread - 1), $maxread);
    $sp1[1] = $seq;
    $sp1[2] = "+";
    my $qstr = ();
    for (1..length($sp1[1])) { $qstr .= $qscores[rand @qscores]; }
    $sp1[3] = $qstr;
    
    my $fakeread = join "\n", @sp1;
    push @reads, $fakeread;
}

@reads = shuffle(@reads);

# For the sake of consistency with the PE version, I may as well have this output a file.
$contaminant =~ s/.fasta//g;
my $out = $reads_file;
$out =~ s/.fastq/_CONTAMINATED_$contaminant.fastq/;
open(OUT, ">", $out) or die "Cannot open output $out: $!\n";
foreach my $read (@reads) {
    print OUT $read;
}
close OUT;
