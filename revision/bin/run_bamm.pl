#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long 'HelpMessage';
use Parallel::ForkManager;

my $usage=<<_EOH_;
------------------------------------------------------------------------------------------
This script is used to run bamm

Usage:
perl run_bamm.pl --genes OG0015521 OG0020241 OG0011939 OG0002874 OG0007437 OG0004968 OG0020625 OG0022213 OG0004329 OG0011340 OG0033097\
--tree r8s_ultrametric_4.txt \
                                                Kang 2021-8-24
------------------------------------------------------------------------------------------
_EOH_
;


GetOptions(
        'genes:s{1,}', \ my @genes, # the genes id
        'tree:s', \my $tre, # the tree file
#        'cont:s', \my $cont, # the cont file
        'help', \ my $help
        );

if ($help || (! @genes) && (! $tre)) {
        die $usage; # all of these options are mandatory requirements
}

if ( (! -e "temp1.pl") || (! -e "PNG.TPM.TMM.sqrt.matrix") || (! -e "traitcontrol.txt")) {
        my $dir="~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm";
        system("cp $dir/temp1.pl ./");
        system("cp $dir/PNG.TPM.TMM.sqrt.matrix ./");
        system("cp $dir/traitcontrol.txt ./");
}

# to prepare files and mkdir dir for bamm

my @cmds;
foreach my $gene (@genes) {
        die "There is no temp1.pl\n" if ! -e "temp1.pl";
        my $cmd1.="perl temp1.pl $gene > $gene.txt;"; # create the expression data for the gene
        my $dir=$gene."_trait";
        $cmd1.="mkdir $dir;";
        $cmd1.="mv $gene.txt $dir;"; # move the expression data to the created dir
        $cmd1.="cp $tre $dir"; # copy the tree to the dir
        system($cmd1);

        # create the control file in the dir
        my $cont_new="$dir/traitcontrol.txt";
        open FIL1, "traitcontrol.txt" or die "can not open traitcontrol.txt\n";
        open FIL2, ">$cont_new" or die "can not create $cont_new\n";
        while (<FIL1>) {
                chomp;
                if (/traitfile\s=\sOG.*/) {
                        print FIL2 "traitfile = $gene.txt\n";
                } elsif (/treefile\s=/) {
                        if (-e "$dir/$tre") {
                                print FIL2 "treefile = $tre\n";
                        } else {
                                die "There is no $tre in $dir/\n";
                        }
                } else {
                        print FIL2 "$_\n";
                }
        }
        my $cmd2="cd $dir;bamm -c traitcontrol.txt";
        push @cmds, $cmd2;
}

my $manager = new Parallel::ForkManager(20);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
