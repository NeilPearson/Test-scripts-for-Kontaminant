#!/usr/bin/perl -w
use strict;
my $file = $ARGV[0];

# This follows a similar premise to blast_outfmt6_count_paired.pl, in that it is meant to count up paired reads in a SAM file. The process should be
# basically the same in principle, with some differences in practice; specifically, I'll be using SAMtools first, and the read details are recorded
# via SAM binary flags rather than as a tag in the query name. 
# ( samtools view -q 60 $results )

#my $sam = `source samtools-1.2.0; samtools view -q 60 $file`;
my $sam = `source samtools-1.2.0; samtools view $file`;
my @sam = split /\n/, $sam;

# Now I'll need to do the same sort of buffered loop as the one in the other script.
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
foreach my $line (@sam) {
    $i ++;
    chomp $line;
    if (($line) || ($i == $lc)) {
        my @sp = split /\t/, $line;
        $query = $sp[0];
        if ($prev_query eq 'Nothing yet!') { $prev_query = $query; }
        
        # This is where we check the query buffer.
        # If the previous query is not the same as the current one, it means we've moved to another one, which means we've
        # successfully picked out all the records (if any) that match the previous one.
        # Take a look, pick the best one, and add it to the results thus far.
        # 'Best' here is indicated by the flag marking the alignment as a primary alignment.
        if (($query ne $prev_query) || ($i == $lc)) {
            if (@buffer) {
                # Modified to have 2 concurrent counts, as described above!
                my $bestline_r1 = ();   my $bestline_r2 = ();
                ALIGNMENT: foreach my $passedline (@buffer) {
                    my $p_log_evalue = $passedline->[10];
                    # We get the read this query is in from the numeric flag field.
                    my $flag = $passedline->[1];
                    my $bin = sprintf("%012b", $flag);
                    my @bits = split //, $bin;
                    # If read isn't mapped, we're not interested.
                    if ($bits[2]) { next ALIGNMENT; }
                    
                    # Get only primary alignments - check flag
                    if($bits[8]) { next ALIGNMENT; }
                    
                    # The read number can be found from the 7th/8th bit; 7th being read 1, 8th being read 2.
                    if ($bits[6]) {
                        $bestline_r1 = $passedline;
                    }
                    elsif ($bits[7]) {
                        $bestline_r2 = $passedline;
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
        
        # If we've got this far, this alignment has passed the filters and we can add it to the buffer.
        push @buffer, \@sp;
    }
}

print "Both reads:\t$count\nR1 only:\t$r1only\nR2 only:\t$r2only\n";



