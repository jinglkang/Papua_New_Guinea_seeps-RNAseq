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
```
4091 single copy    
contruct phylogenetic tree     
```bash
#kang1234@celia-PowerEdge-T640 Tue Nov 30 11:46:22 ~/CO2-seeps/annotation
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s single_copy.concatenated.phy -o Spottedgar -n single_copy.concatenated > raxml.process 2>&1 &
```
