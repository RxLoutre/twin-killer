## removesmalls.pl
#!/usr/bin/perl
use strict;
use warnings;

#USAGE EXAMPLE
#perl removesmalls.pl 200 contigs.fasta > contigs-l200.fasta

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}
