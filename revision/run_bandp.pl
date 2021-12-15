#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @fq=<*_1.fastq.gz>;
my @cmds;
foreach my $fq (@fq) {
        (my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
        (my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
        my $bait=$spe."_100_single_copy";
#        my $cmd1="./rmrep.pl -taxalist=$ind";
        my $ind_results=$ind."_results";
        my $cmd2="./bandp.pl -query=$bait -subject=$ind";
        if (-d "$ind_results") {
                next;
        } else {
#                push @cmds, $cmd1;
                system($cmd2);
        }
}
