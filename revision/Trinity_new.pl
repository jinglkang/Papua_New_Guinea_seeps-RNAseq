#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $usage=<<_EOH_;;
------------------------------------------------------------------------------------------
run trinity
Trinity.pl -species=Acura -output=trinity
																			Kang 2021-8-25
------------------------------------------------------------------------------------------
_EOH_
;

# assemble acoording gene

GetOptions('species:s', \my $species,
    'output:s', \ my $output);
if (!($output && $species)) {
    die $usage;
}

my @results=<$species*_results>;
foreach my $result (@results) {
        (my $spec)=basename($result)=~/(.*)_results/;
    (my $out)="$output/$spec";
	next if -d $out;
    system("mkdir -p $out");
    chdir $result;
    my $fa="trinity_out_dir.Trinity.fasta";
    my $map="trinity_out_dir.Trinity.fasta.gene_trans_map";
    my @fq=<*.fq>;
    foreach my $fastq (@fq) {
        my $cmd="docker run --rm -e LOCAL_USER_ID=`id -u \$USER` ";
                $cmd.="-u 1003:1003 -v \$SHARED_DIR:\$SHARED_DIR -w `pwd` sigenae\/drap ";
        $cmd.="Trinity --seqType fq --max_memory 300G --jaccard_clip --single $fastq --full_cleanup --CPU 24";
        system($cmd);
        (my $gene)=basename($fastq)=~/(.*)\.fq/;
        (my $new_fa)=$gene.".fasta";
        system("mv $fa ../$out/$new_fa");
        system("rm $map");
    }
        chdir("../");
}
