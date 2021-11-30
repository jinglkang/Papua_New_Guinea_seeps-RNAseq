# identify the single copy genes to construct a phylogentic tree for BAMM
## annotate the six species genes: align to uniprot and selected the longest when sequences were aligned to uniprot id with the same genes
working dir (SNORLAX): ~/CO2-seeps/annotation    
```bash
for fa in *.fasta; do ./annotate --fasta ${fa}; done
```
\*.ano were the annotations of sequences with blast hits    
\*.ano.longest were the annotations of sequences with blast hits, and select the longest when sequences were aligned to uniprot id with the same genes     
\*.ano.longest.fasta is the fasta file of \*.ano.longest      
## Orthofinder: my own work station (~/Desktop/PapueNewGuinea-new/longest_pep)
copy all species longest pep into this directory
```bash
#Kang@fishlab3 Tue Nov 30 09:00:32 /media/HDD/cleaner_fish/genome/gene_family_2
cp Longest_Fugu_pep.fasta Longest_Japanese_Medaka_pep.fasta Longest_Platyfish_pep.fasta Longest_Spotted_gar_pep.fasta Longest_Stickleback_pep.fasta Longest_Zebrafish_pep.fasta ~/Desktop/PapueNewGuinea-new/longest_pep/
#Kang@fishlab3 Tue Nov 30 09:08:11 ~/Desktop/PapueNewGuinea-new/longest_pep
scp kang1234@147.8.76.155:~/CO2-seeps/annotation/*.ano.longest.fasta ./
```
## Change the name of each fasta file, such as Acura.ano.longest.fasta --> Acura.fast, Longest_Japanese_Medaka_pep.fasta --> Medaka.fasta
```bash
mkdir input_pep
cp *.fasta input_pep/
# run orthofinder
orthofinder -f input_pep -a 32
```
the results: ~/Desktop/PapueNewGuinea-new/longest_pep/input_pep/OrthoFinder/Results_Nov30/Single_Copy_Orthologue_Sequences    
4098 single copy genes     
## concatenated the fasta sequence to construct a phylogenetic tree by raxml
copy the single copy fasta file to a new dir "single_copy" and change the name   
```perl
#!/usr/bin/perl
use strict;
use warnings;

system("mkdir single_copy") unless -e "single_copy";
my @fas=<Single_Copy_Orthologue_Sequences/*.fa>;
foreach my $fa (@fas) {
        (my $name)=$fa=~/Single_Copy_Orthologue_Sequences\/(.*\.fa)/;
        open FIL1, "$fa" or die "can not open $fa";
        open FIL2, ">single_copy/$name" or die "can not create single_copy/$name";
        while (<FIL1>) {
                chomp;
                if (/>/) {
                        s/>//;
                        if (/ENSTRUG/) {
                                my $header="Fugu_".$_;
                                print FIL2 ">$header\n";
                        } elsif (/ENSORLG/) {
                                my $header="Medaka_".$_;
                                print FIL2 ">$header\n";
                        } elsif (/ENSXMAG/) {
                                my $header="Platyfish_".$_;
                                print FIL2 ">$header\n";
                        } elsif (/ENSLOCG/) {
                                my $header="Spottedgar_".$_;
                                print FIL2 ">$header\n";
                        } elsif (/ENSGACG/) {
                                my $header="Stickleback_".$_;
                                print FIL2 ">$header\n";
                        } elsif (/ENSDARG/) {
                                my $header="Zebrafish_".$_;
                                print FIL2 ">$header\n";
                        } else {
                                my $header=$_;
                                print FIL2 ">$header\n";
                        }
                } else {
                        print FIL2 "$_\n";
                }
        }
        close FIL1;
        close FIL2;
}
```
```bash
perl temp1.pl
```
single copy fasta file will be in single_copy/
## align and trimal to select the conserved regions (only keep the pep sequences longer than 50)
```perl
#!/usr/bin/perl
use strict;
use warnings;

my (%seq, %hash);
my $spe; my @spes;
my @fasta=<single_copy/*.fa>;
foreach my $fasta (@fasta) {
        open FAS, "$fasta" or die "can not open $fasta\n";
        my $flag;
        while (<FAS>) {
                chomp;
                unless (/>/) {
                        my $len=length($_);
                        if ($len<50) {
                                $flag=1;
                        }
                }
        }
        close FAS;
        if ($flag==1) {
                system("rm $fasta");
        } else {
                my $align=$fasta.".align";
                system("muscle -in $fasta -out $align");
                system("mv $align $fasta");
                my $new=$fasta."trimal";
                system("trimal -in $fasta -out $new -gt 0.8 -st 0.001 -cons 60");
                open FIL, "$new" or die "can not open $new\n";
                while (<FIL>) {
                        chomp;
                        if (/>(.*?)(_|-)/) {
                                $spe=$1;
                                $hash{$spe}++;
                                push @spes, $spe if $hash{$spe}==1;
                        } else {
                                $seq{$spe}.=$_;
                        }
                }
                system("rm $new");
        }
}

foreach my $spe (@spes) {
        print ">$spe\n$seq{$spe}\n";
}
```
```bash
perl temp2.pl >single_copy.concatenated.fasta
fasta2phy.pl single_copy.concatenated.fasta >single_copy.concatenated.phy
```
4091 single copy    
contruct phylogenetic tree     
```bash
#kang1234@celia-PowerEdge-T640 Tue Nov 30 11:46:22 ~/CO2-seeps/annotation
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s single_copy.concatenated.phy -o Spottedgar -n single_copy.concatenated -T 24 > raxml.process 2>&1 &
```
***
## second (annotate to same uniprot_id and select the longest)
working dir (SNORLAX): ~/CO2-seeps/annotation/second   

```bash
cp Acura.fasta Apoly.fasta Daru.fasta Ocomp.fasta Padel.fasta Pmol.fasta second/
for fa in *.fasta; do ./annotate1 --fasta ${fa}; done
```
## Orthofinder: my own work station (~/Desktop/PapueNewGuinea-new/longest_pep_sec)
```bash
orthofinder -f input_pep -a 32
cd ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30
cp ~/Desktop/PapueNewGuinea-new/longest_pep/input_pep/OrthoFinder/Results_Nov30/temp1.pl ./
cp ~/Desktop/PapueNewGuinea-new/longest_pep/input_pep/OrthoFinder/Results_Nov30/temp2.pl ./
perl temp1.pl
perl temp2.pl >single_copy.concatenated.fasta
```
2686 single-copy genes, 2683 single-copy genes in the final    

```bash
fasta2phy.pl single_copy.concatenated.fasta >single_copy.concatenated.phy
scp single_copy.concatenated.phy kang1234@147.8.76.155:~/CO2-seeps/annotation/second
# SNORLAX: ~/CO2-seeps/annotation/second
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s single_copy.concatenated.phy -o Spottedgar -n single_copy.concatenated -T 24 > raxml.process 2>&1 &
```
***
## MCMCtree
## Extract the nucleotide sequence
first to Get the id   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my @fasta=<single_copy/*.fa>;
foreach my $fa (@fasta) {
	open FIL, "$fa" or die "can not open $fa\n";
	while (<FIL>) {
		chomp;
		if (/>/) {
			s/>//;
			if (! /\_ENS/) {
				(my $spe)=$_=~/(.*?)\_/;
				my $output=$spe.".id.txt";
				open FILE, ">>$output" or die "can not create $output\n";
				print FILE "$_\n";
				close FILE;
			} else {
				(my $spe)=$_=~/(.*?)\_/;
				s/(.*)\_//;
				my $output=$spe.".id.txt";
				open FILE, ">>$output" or die "can not create $output\n";
				print FILE "$_\n";
				close FILE;
			}
		}
	}
}
```
the id were save in \*.id.txt   
```bash
perl temp3.pl
```
## prepare the cds of Platyfish
working dir: ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30   
```bash
perl prep_Ensembl_cds.pl --fasta Platyfish.cds.all.fa
```
the result: Platyfish_cds.fasta   
```bash
perl Ensemble_longest_pep.pl --fasta Fugu.id.fa Medaka.id.fa Platyfish_cds.fasta Spottedgar.id.fa Zebrafish.id.fa Stickleback.id.fa Spottedgar.id.fa
#Kang@fishlab3 Tue Nov 30 19:41:23 ~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_nuc
cp *.fa ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30
cat Longest_*id.fa Longest_Platyfish_cds.fasta *_nuc.fa >nuc_all_spe.fa
```
## put all cds nucleotide sequences in single_copy_nuc/
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $gene;
open FILE, "nuc_all_spe.fa" or die "can not open nuc_all_spe.fa\n";
while (<FILE>) {
        chomp;
        if (/>/) {
                s/>//;
                $gene=$_;
        } else {
                $hash{$gene}.=$_;
        }
}

my @fasta=<single_copy/*.fa>;
system("mkdir single_copy_nuc");
foreach my $fa (@fasta) {
        open FIL, "$fa" or die "can not open $fa\n";
        (my $name)=$fa=~/.*\/(.*)/;
        open FIL1, ">single_copy_nuc/$name" or die "can not create single_copy_nuc/$name\n";
        while (<FIL>) {
                chomp;
                if (/>/) {
                        s/>//;
                        if (/\_(ENS.*)/) {
                                (my $gene)=$1;
                                print FIL1 ">$_\n$hash{$gene}\n";
                        } else {
                                (my $gene)=$_;
                                print FIL1 ">$_\n$hash{$gene}\n";
                        }
                }
        }
        close FIL;
        close FIL1;
}
```
```bash
perl temp4.pl
```
all nucleotide sequences were in single_copy_nuc/    
```bash
cp temp2.pl temp5.pl
```
## according the steps of preparing paml input
working dir: ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30   
```bash
mkdir paml_input
less Longest_Fugu.id.fa|perl -alne 'if (/>/){s/>//;$na="Fugu_".$_;print ">$na"}else{print}' >paml_input/Longest_Fugu.id.fa
less Longest_Zebrafish.id.fa|perl -alne 'if (/>/){s/>//;$na="Zebrafish_".$_;print ">$na"}else{print}' >paml_input/Longest_Zebrafish.id.fa
less Longest_Spottedgar.id.fa|perl -alne 'if (/>/){s/>//;$na="Spottedgar_".$_;print ">$na"}else{print}' >paml_input/Longest_Spottedgar.id.fa
less Longest_Stickleback.id.fa|perl -alne 'if (/>/){s/>//;$na="Stickleback_".$_;print ">$na"}else{print}' >paml_input/Longest_Stickleback.id.fa
less Longest_Medaka.id.fa|perl -alne 'if (/>/){s/>//;$na="Medaka_".$_;print ">$na"}else{print}' >paml_input/Longest_Medaka.id.fa
less Longest_Platyfish_cds.fasta|perl -alne 'if (/>/){s/>//;$na="Platyfish_".$_;print ">$na"}else{print}' >paml_input/Longest_Platyfish_cds.fasta
cp *_nuc.fa paml_input/
less /media/HDD/cleaner_fish/genome/gene_family_2/Longest_Zebrafish_pep.fasta|perl -alne 'if (/>/){s/>//;$na="Zebrafish_".$_;print ">$na"}else{print}' >paml_input/Longest_Zebrafish_pep.fasta
```
prepare the orth list    
```bash
perl temp1.pl >ortho_list.txt
vi correlation.txt
perl prepare_input_paml.pl --input ortho_list.txt --seq_dir . --cor_list correlation.txt --output .
```
correlation.txt:    
```bash
refer   Longest_Zebrafish_pep.fasta
Zebrafish       Longest_Zebrafish.id.fa
Medaka  Longest_Medaka.id.fa
Spottedgar      Longest_Spottedgar.id.fa
Fugu    Longest_Fugu.id.fa
Platyfish       Longest_Platyfish_cds.fasta
Ocomp   Ocomp_nuc.fa
Stickleback     Longest_Stickleback.id.fa
Daru    Daru_nuc.fa
Apoly   Apoly_nuc.fa
Acura   Acura_nuc.fa
Pmol    Pmol_nuc.fa
Padel   Padel_nuc.fa
```
