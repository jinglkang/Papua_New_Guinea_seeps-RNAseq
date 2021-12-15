#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

my $opt = GetOptions( 'query:s', \my $query,
	'subject:s', \my $subject); #set command line options

# getbest.pl -query="Acura.id" -subject="Acura10"
if (!$query && ! $subject) {#check for the required inputs
   die "pls give the query and subject\n";
}

my $resultdirnf = $query . ".resultnf_new";
my @sub = split (/\s/, $subject);
my @cmds;
foreach my $ind (@sub) {
	my $cmd="./getbest_new.pl -query=$query -subject=$ind";
	push @cmds, $cmd;
}

my $manager = new Parallel::ForkManager(12);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;

open QUERY, "< $query.fa" or die ("Cannot open $query for reading ($!)");
while (<QUERY>) {
	chomp;
	if (/>/) {
		s/>//;
		my $gene=$_;
		my $cmd="cat result/$resultdirnf/$gene*.fas >$gene.fas";
		system($cmd);
		system("rm result/$resultdirnf/$gene*.*.fas");
		system("mv $gene.fas result/$resultdirnf/");
	}
}
