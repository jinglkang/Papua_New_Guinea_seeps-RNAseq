#!/usr/bin/perl
use strict;
use warnings;

my @fq=<*_1.fastq.gz>;
foreach my $fq (@fq) {
	(my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
	(my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
	my $bait=$spe."_100_single_copy";
	my $cmd1="./rmrep.pl -taxalist=$ind";
	my $cmd2="./bandp.pl -query=$bait -subject=$ind";
	if (-e "$spe.fas") {
		system($cmd2);
	} else {
		system($cmd1);
		system($cmd2);
	}
}

my $cmd3="Trinity.pl -input=. -output=trinity";
system($cmd3);
foreach my $fq (@fq) {
	(my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
	(my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
	my $bait=$spe."_100_single_copy";
	my $cmd4="./getbest.pl -query=$bait -subject=$ind";
	system($cmd4);
	chdir "result";
	my $resultnf=$bait.".resultnf";
	my $cmd5="./reblast.pl -query $resultnf -database $spe";
	system($cmd5);
}
