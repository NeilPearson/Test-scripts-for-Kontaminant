#!/usr/bin/perl -w
use strict;
$|++;
my $reads_file = $ARGV[0];
my $substitution_freq = $ARGV[1];
my $indel_freq = $ARGV[2];
my $max_indel_length = $ARGV[3];

# This script introduces random mutations into FASTQ NGS reads, and prints 'em back out.
# Both substitutions and indels are possible.
# Frequencies are given as percentages.
# When an insertion is made, bases are picked at random and inserted at a random position, and then the read is shortened back to its original length.
# When a deletion is made, bases are added at random at the end of the read to make it back up to the original length.
# We'll ignore qscores completely, because none of the software used in these experiments actually makes use of it.

# Since freqs are expected as percentages, and I'll be working with random numbers in the 0-1 range, I want to divide the freq figures by 100.
$substitution_freq = $substitution_freq / 100;
$indel_freq = $indel_freq / 100;

open(READS, "<", $reads_file) or die "Cannot open reads file $reads_file: $!\n";
my @readslines = <READS>;   close READS;

my @buffer = ();
my @bases = qw/A C G T/;
my $maxread = 0;
foreach my $line (@readslines) {
    chomp $line;
    push @buffer, $line;
    if (@buffer >= 4) {
        # This is a single FASTQ read. I can now introduce random mutations to its sequence.
        my @seq = split '', $buffer[1];
        my $orig_length = @seq;
        for my $i (1..@seq) {
            $i --;
            my $b = $seq[$i];
            
            # Work out if there's going to be a substitution mutation on this base.
            if (rand() <= $substitution_freq) {
                # The substitution can't be the existing base; that would reduce the real sub rate by 1/4. 
                my @rbases = not_this_base($b);
                $seq[$i] = $rbases[rand 3];
            }
            
            # Then work out if there's going to be an indel starting here.
            # (To keep conflict with the sub rate to a minimum, I'll make indels begin from the following base).
            if (rand() <= $indel_freq) {
                # Indels can be achieved using splice.
                my $len = int rand($max_indel_length);
                # make the decision between insertion and deletion (equal weight for now; will change later if necessary).
                if (rand() < 0.5) {
                    # Insertion - pick random bases and insert them
                    my @insert = ();
                    for (my $i=1; $i <= $len; $i++) {
                        push @insert, $bases[rand 4];
                    }
                    splice(@seq, $i+1, 0, @insert);
                }
                else {
                    # Deletion - just remove some bases
                    splice(@seq, $i+1, $len);
                }
            }
            
            # Now check that the sequence matches the original length, and act appropriately if it doesn't.
            if (@seq > $orig_length) {
                # Remove some bases because it's too long
                while (@seq > $orig_length) {
                    pop @seq;
                }
            }
            elsif (@seq < $orig_length) {
                # Add some bases because it's too short
                my @insert = ();
                for (my $i=1; $i <= ($orig_length - @seq); $i++) {
                    push @insert, $bases[rand 4];
                }
                push @seq, @insert;
            }
        }
        
        # Now join the sequence back together and print the read data out.
        $buffer[1] = join '', @seq;
        foreach my $outline (@buffer) {
            print "$outline\n";
        }
        
        # CLEAR THE BUFFER
        @buffer = ();
    }
}
@readslines = ();

sub not_this_base {
    # Returns a list of bases that are NOT the one supplied.
    # We use this in order to avoid reducing the effective substitution rate by 1/4, by substituting in the same base.
    my $base = shift;
    
    my @outbases = ();
    for my $b (@bases) {
        unless ($base eq $b) { push @outbases, $b; }
    }
    return @outbases;
}

