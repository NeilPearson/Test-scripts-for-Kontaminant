#!/usr/bin/perl -w
use strict;
my $file = $ARGV[0];

# Count up the number of BLASTn alignments fitting the parameters we're looking for.
# Those will be hard-coded here, for the sake of convenience.
# I'd like to add capability to check for cases where both pairs of a read are represented, but I'm not too sure how to do that yet.
# FASTA files should now have a pair tag in their ID string which will be visible even after alignment. They should also be right next to each other.
# I'll need a results file to verify that though. Go and generate one. (I'll need the thing to interleave FASTA first...)
# If that does hold true, I'll need to make a wrapper loop on all the stuff I have here already. Or something to that effect, anyway.

#ÊThis modification is actually pretty easy. In my chosen dataset, reads are marked by a /1 and /2 tag. I've made my contaminant reads have those same tags.
# I've also put the reads in such an order (interleaved) that we'll get both reads of a pair as a contiguous block of lines in the output.
# That means I just need to add a function to strip the /1 and /2 parts from the query ID, in order to effect that.

# However, things get a little more complicated after that. I'll need to count both R1 and R2 results separately, and then make a final decision on whether to count
# a pair as mapped once I have both numbers. Still pretty easy stuff.

open(BLAST, "<", $file) or die "ERROR: Unable to open file $file\n$!\n";
my @blast = <BLAST>;
close BLAST;

# Columns in outfmt 6 are as follows:
# 0      1      2      3      4        5       6      7    8      9    10     11
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# We need to get the best (primary) alignment for each query sequence. Use my buffer trick.
my $count = 0;  my $r1only = 0; my $r2only = 0;
my $prev_query = "Nothing yet!";
my $query = "Nothing yet!";
my @buffer = ();
my $i = 0;
my $lc = `wc -l $file`;
my @lc = split /\s/, $lc;
$lc = $lc[0];
chomp $lc;
chomp $file;
open (INFILE, "<", "$file") or die "Could not open $file: $!";
LINE: while(my $line = <INFILE>) {
    $i ++;
    chomp $line;
    if (($line) || ($i == $lc)) {
        my @sp = split /\t/, $line;
        $query = $sp[0];
        # Do the removal of /1 and /2 tags, as described above.
        my @qry = split /\//, $query;
        $query = $qry[0];
        
        if ($prev_query eq 'Nothing yet!') { $prev_query = $query; }
        
        # This is where we check the query buffer.
        # If the previous query is not the same as the current one, it means we've moved to another one, which means we've
        # successfully picked out all the records (if any) that match the previous one.
        # Take a look, pick the best one, and add it to the results thus far.
        # (I'm not completely sure, but I think 'best' in this case can be most effectively quantified by lowest e-value)
        if (($query ne $prev_query) || ($i == $lc)) {
            if (@buffer) {
                # Modified to have 2 concurrent counts, as described above!
                my $bestline_r1 = ();   my $bestline_r2 = ();
                my $bestscore_r1 = 0;   my $bestscore_r2 = 0;
                
                foreach my $passedline (@buffer) {
                    my $p_log_evalue = $passedline->[10];
                    my $mapped_query = $passedline->[0];
                    my @qry = split /\//, $mapped_query;
                    my $read = $qry[1];
                    if ($read == 1) {
                        if (($p_log_evalue < $bestscore_r1) || (!$bestscore_r1)) {
                            $bestline_r1 = $passedline;
                            $bestscore_r1 = $p_log_evalue;
                        }
                    }
                    elsif ($read == 2) {
                        if (($p_log_evalue < $bestscore_r2) || (!$bestscore_r2)) {
                            $bestline_r2 = $passedline;
                            $bestscore_r2 = $p_log_evalue;
                        }
                    }
                }
                
                # Now implement the screening rule mentioned earlier.
                # This need not be complex - we simply require a sufficiently well aligned result for read 1 AND read 2. 
                if (($bestline_r1) && ($bestline_r2)) {
                    $count ++;
                }
                # That said, it would also be useful to count the number of instances where only read 1 or 2 is present separately.
                # Note that this is FizzBuzz. Funny stuff.
                elsif ($bestline_r1) { $r1only ++; }
                elsif ($bestline_r2) { $r2only ++; }
            }
            # CLEAR THE BUFFER!
            @buffer = ();
        }
        
        # Put the current query ID to $prev_query so we know if it's time to take a look in the buffer next time round.
        $prev_query = $query;
        
        # Filter out insufficiently long/good results
        if ($sp[2] ne '100.00') { next LINE; }
        if ($sp[3] < 50) { next LINE; }
        
        # If we've got this far, this alignment has passed the filters and we can add it to the buffer.
        push @buffer, \@sp;
    }
}
close INFILE;

print "Both reads:\t$count\nR1 only:\t$r1only\nR2 only:\t$r2only\n";

