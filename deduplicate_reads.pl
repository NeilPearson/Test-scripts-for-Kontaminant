#!/usr/bin/perl -w
use strict;
use File::Basename;

# A script to remove PCR duplicates from a reads file.
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $check_n_characters = $ARGV[2];

if ((!$infile) || (!$outfile)) {
    die "Usage:
    perl deduplicate_reads.pl <input_file> <output_file> [<up_to_n_bases>]\n";
}

my $lines_per_read = 4;
my ($basename, $parentdir, $extension) = fileparse($infile, qr/\.[^.]*$/);
if ($extension =~ /\.fasta$/) { $lines_per_read = 2; }

open INFILE, "<", $infile or die "ERROR: Could not open file $infile: $!\n";
my @buffer = ();    my %seen_reads = ();    my @output_reads = ();
while (my $line = <INFILE>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer >= $lines_per_read) {
        # Check if a hash exists for the sequence line
        # If not, write it, and set a value in the corresponding hash.
        my $seq = $buffer[1];
        my $checkseq = $seq;
        # In cases where we know the read quality drops towards the end of the reads, we may want to work out PCR duplicates
        # based on only the first n reads. There is a parameter in the config file for that. (Leave it blank to check everything).
        if ($check_n_characters) { $checkseq = substr($seq, 0, $check_n_characters); }
        
        if (!$seen_reads{$checkseq}) {
            push @output_reads, @buffer;
            $seen_reads{$checkseq} = 1;
        }
        @buffer = ();
    }
}
close INFILE;

open (OUT, ">", $outfile) or die "ERROR: Cannot open output file $outfile\n"; 
foreach my $l (@output_reads) { print OUT "$l\n"; }
close OUT;

# Remove trailing newlines from the end of the file
my $gibberish = '/./,$!d;/^\n*$/{$d;N;};/\n$/ba';
`sed -e :a -e '$gibberish' $outfile`;








