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
### Estimate a phylogenetic tree based on the pep sequences of single copy genes between 6 PNG fish species
```bash
# Kang@fishlab3 Wed Dec 01 21:06:13 ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30
less single_copy.concatenated.phy|perl -lane '$i++;print "6\t1437114" if $i==1;if(/Apoly|Daru|Ocomp|Acura|Padel|Pmol/){print}' >single_copy.PNG.concatenated.phy
scp single_copy.PNG.concatenated.phy kang1234@147.8.76.155:~/CO2-seeps/annotation/second
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s  single_copy.PNG.concatenated.phy -o Ocomp -n single_copy.PNG.concatenated -T 24 > raxml.process 2>&1 &
```
***
### concatenate the aligned genes according to final_orth_input_paml.txt
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %seq; my @spe; my $i;
open LIST, "final_orth_input_paml.txt" or die "$!\n";
while (<LIST>) {
	chomp;
	$i++;
	my $orth=$_;
	open FASTA, "$orth/final_alignment.fa" or die "can not open $orth/final_alignment.fa\n";
	my $spe;
	while (<FASTA>) {
		chomp;
		if (/>/) {
			s/>//;
			$spe=$_;
			push @spe, $spe if $i==1;
		} else {
			$seq{$spe}.=$_;
		}
	}
}

foreach my $spe (@spe) {
	my $seqs=$seq{$spe};
	print ">$spe\n$seqs\n";
}
```
save the concatenated sequences in single_copy.cds.concatenated.fasta    
```bash
perl temp2.pl >single_copy.cds.concatenated.fasta
fasta2phy.pl single_copy.cds.concatenated.fasta >single_copy.cds.concatenated.phy
```
single_copy.cds.concatenated.phy (920,751) will be used as mcmctree input sequence file   

#### Rough estimation of the substitution rate
estimate the best nucleotide module (working dir: ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30/paml_input)    
```bash
java -jar ~/software/jmodeltest-2.1.10/jModelTest.jar -d single_copy.cds.concatenated.fasta -g 4 -i -f -AIC -BIC -a
```
Best Models:    
||Model|f(a)|f(c)|f(g)|f(t)|kappa|titv|Ra|Rb|Rc|Rd|Re|Rf|pInv|gamma|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|AIC|GTR+I+G|0.24|0.28|0.27|0.21|0.00|0.00|1.370|3.660|1.186|1.241|5.741|1.000|0.30|0.80|
|BIC|GTR+I+G|0.24|0.28|0.27|0.21|0.00|0.00|1.370|3.660|1.186|1.241|5.741|1.000|0.30|0.80|
##### put all files related to mcmctree to a same folder
the input tree file: PNG2.tree    
control file: baseml.ctl   
```
12 1
((Zebrafish,(Ocomp,((Fugu,Stickleback),((Platyfish,Medaka),(((Acura,Apoly),(Padel,Pmol)),Daru)))))'@(1.4885, 1.652)',Spottedgar);
```
run baseml   
```bash
# Kang@fishlab3 Wed Dec 01 10:16:19 ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep
mkdir mcmctree
cd mcmctree
mkdir baseml; cd baseml; baseml;
```
the control file of baseml: baseml.ctl in ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/mcmctree/baseml    
Substitution rate is per time unit (in file "mlb"): 0.017619 +- 3501.972727 (use this value to set the prior for the mean substitution rate in the Bayesian analysis)     

#### Estimation of the gradient and Hessian
the input tree file: mcmc.tree    
```
12 2
((Zebrafish,(Ocomp,((Fugu,Stickleback),((Platyfish,Medaka),(((Acura,Apoly),(Padel,Pmol)),Daru)))'B(.969, 1.509)'))'B(1.4885, 1.652)',Spottedgar);
```
in the control file (mcmctree.ctl): set "usedata = 3", "ndata = 1" (no partition)    
```bash
~/software/paml4.9j/src/mcmctree
```
The results were written to a file called out.BV   
Rename the out.BV file as in.BV   
```bash
mv out.BV in.BV
```
in the control file (mcmctree.ctl): set "usedata = 3" (perform Bayesian estimation of divergence times using in.BV)
```bash
~/software/paml4.9j/src/mcmctree
cat in.BV rst2 > 1.txt
mv 1.txt in.BV
```

#### Estimation of divergence times with the approximate likelihood method
in the control file (mcmctree.ctl): set "usedata = 2" (perform Bayesian estimation of divergence times using the approximate likelihood method), "ndata = 1" (no partition), "rgene_gamma = 1 55.6" (see the manual to calculate),     
```bash
~/software/paml4.9j/src/mcmctree
```
#### run again to check the convergence
```bash
#Kang@fishlab3 Wed Dec 01 22:44:10 ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/mcmctree/baseml
cp mcmctree2.ctl in.BV ../sec
#Kang@fishlab3 Wed Dec 01 22:45:21 ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/mcmctree/sec
nohup ~/software/paml4.9j/src/mcmctree mcmctree2.ctl >mcmctree.process 2>&1 &
```
#### based on RAxML_bestTree.single_copy.PNG.concatenated for a time-calibrated ultrametric phylogenetic tree
```bash
# Kang@fishlab3 Thu Dec 02 15:22:56 ~/Desktop/PapueNewGuinea-new/bamm
cp /media/HDD/cleaner_fish/genome/gene_family/cafetutorial_prep_r8s.py ./
scp kang1234@147.8.76.155:~/CO2-seeps/annotation/second/RAxML_bestTree.single_copy.PNG.concatenated ./
python cafetutorial_prep_r8s.py -i RAxML_bestTree.single_copy.PNG.concatenated -o r8s_ctl_file.txt -s 638853 -p 'Ocomp,Daru' -c '109.38'
r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > r8s_ultrametric.txt
cp ~/software/bamm/examples/traits/fishsize/traitcontrol.txt ./
```
##### based on PNG.TPM.TMM.sqrt.matrix to calculate the average reads nb per gene per species
less average_normal_expression.pl   
```perl
# ~/Documents/2019/香港大学/co2_seeps/EVE_release/run_again
#!/usr/bin/perl
use strict;
use warnings;

my %hash; my @spes=qw(Ocomp Daru Pmol Padel Acura Apoly);
print "\tOcomp\tDaru\tPmol\tPadel\tAcura\tApoly\n";
open FIL1, "PNG.TPM.TMM.sqrt.matrix" or die "can not open PNG.TPM.TMM.sqrt.matrix\n";
while (<FIL1>) {
	chomp;
	my @a=split /\t/;
	if (/^\s+/) {
		for (my $i = 1; $i < @a; $i++) {
			(my $spe)=$a[$i]=~/(\D+)/;
			$hash{$spe} = [] unless exists $hash{$spe};
			push @{$hash{$spe}}, $i;
		}
	} else {
		my $gene=$a[0]; my $info;
		foreach my $spe (@spes) {
			my $total;
			my $num=@{$hash{$spe}};
			foreach my $i (@{$hash{$spe}}) {
				$total+=$a[$i];
			}
			my $mean=$total/$num;
			$mean=sprintf("%.2f",$mean);
			$info.=$mean."\t";
		}
		$info=~s/\s+$//;
		print "$gene\t$info\n";
	}
}
```
```bash
# ~/Documents/2019/香港大学/co2_seeps/EVE_release/run_again
perl average_normal_expression.pl >PNG_average_normal_expression.txt
```
***
#### assemble transcriptome for each individual
kang1234@celia-PowerEdge-T640 Fri Dec 03 09:09:15 \~/CO2-seeps/high_index/paired/kraken/merge   
vi run_trinity.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my $spe=$ARGV[0];
my @fastq=<*_1.fastq.gz>;
my $ori1="trinity_out_dir.Trinity.fasta"; # Trinity original result fasta file
my $ori2="trinity_out_dir.Trinity.fasta.gene_trans_map"; # Trinity original gene transcript map file
foreach my $fastq (@fastq) {
	if (-e "trinity_out_dir") {
		system("rm -rf trinity_out_dir");
	}
	next if ! $fastq=~/${spe}/;
	(my $ind)=$fastq=~/(.*)_/;
	my $left=$fastq;
	my $right=$ind."_2.fastq.gz";
	system("export SHARED_DIR=\$PWD");
	my $result1=$ind.".Trinity.fasta";
	my $result2=$ind.".Trinity.fasta.gene_trans_map";
	my $cmd1="docker run --rm -e LOCAL_USER_ID=`id -u \$USER` ";
	$cmd1.="-u 1003:1003 -v \$SHARED_DIR:\$SHARED_DIR -w `pwd` sigenae\/drap ";
	$cmd1.="Trinity --full_cleanup --seqType fq --max_memory 400G --bflyHeapSpaceMax 4G --CPU 22 ";
	$cmd1.="--no_normalize_reads --left $left --right $right --SS_lib_type FR"; # exe Trinity
	my $cmd2="mv $ori1 $result1"; # change Trinity fasta file name
	my $cmd3="mv $ori2 $result2"; # change Trinity gene transcript map file name
	unless (-e $result1) {
		system($cmd1);
		if (-e $ori1) {
			system($cmd2);
		} else {
			die "there is no $ori1\n";
		}
		system($cmd3);
	}
}
```

```bash
nohup perl run_trinity.pl Acura >Acura_trinity.process 2>&1 &
```
\[1\] 14420   
***
also do the assembly in my own workstation
vi run_trinity.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my @inds=qw(Pmol7 Pmol8 Pmol9 Pmol22 Pmol23 Pmol24);
my $ori1="trinity_out_dir.Trinity.fasta"; # Trinity original result fasta file
my $ori2="trinity_out_dir.Trinity.fasta.gene_trans_map"; # Trinity original gene transcript map file
foreach my $ind (@inds) {
        if (-e "trinity_out_dir") {
                system("rm -rf trinity_out_dir");
        }
        my $left=$ind."_1.fastq.gz";;
        my $right=$ind."_2.fastq.gz";
        system("export SHARED_DIR=\$PWD");
        my $result1=$ind.".Trinity.fasta";
        my $result2=$ind.".Trinity.fasta.gene_trans_map";
        my $cmd1="docker run --rm -e LOCAL_USER_ID=`id -u \$USER` ";
		$cmd1.="-u 1002:1002 -v \$SHARED_DIR:\$SHARED_DIR -w `pwd` sigenae\/drap ";
        $cmd1.="Trinity --full_cleanup --seqType fq --max_memory 400G --bflyHeapSpaceMax 4G --CPU 28 ";
        $cmd1.="--no_normalize_reads --left $left --right $right --SS_lib_type FR"; # exe Trinity
        my $cmd2="mv $ori1 $result1"; # change Trinity fasta file name
        my $cmd3="mv $ori2 $result2"; # change Trinity gene transcript map file name
        unless (-e $result1) {
                system($cmd1);
                if (-e $ori1) {
                        system($cmd2);
                } else {
                        die "there is no $ori1\n";
                }
                system($cmd3);
        }
}
```
```bash
nohup perl run_trinity.pl >Pmol_trinity.process 2>&1 &
```
\[1\] 17065   

## Maybe we need not to assembly, build a phylogenetic tree based on SNPs dataset
use Apoly geneome as reference    
```bash
# Kang@fishlab3 Fri Dec 03 16:59:48 ~/Desktop/PapueNewGuinea-new/merge_clean
hisat2-build -p 30 -f apoly_primary_v1.fasta Apoly
mkdir Apoly-align
nohup perl hisat_sorted_bam.pl >map_sorted.process 2>&1 &
# [1] 6284
```
***
## Gene capture
select 100 single copy genes (~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30)   
```perl
#!/usr/bin/perl
use strict;
use warnings;

# select 100 single coply genes as the bait gene to capture
my %hash; my $gene;
open FIL1, "nuc_all_spe.fa" or die "can not open nuc_all_spe.fa\n";
while (<FIL1>) {
	chomp;
	if (/>/) {
		s/>//;
		$gene=$_;
	} else {
		$hash{$gene}.=$_;
	}
}

open FIL2, "paml_input/ortho_list.txt" or die "can not open paml_input/ortho_list.txt\n";
my $k;
while (<FIL2>) {
	chomp;
	$k++;
	last if $k==101;
	my @a=split;
	my $orth=$a[0];
	for (my $i = 1; $i < @a; $i++) {
		unless ($a[$i]=~/ENS/) {
			my $gene=$a[$i];
			(my $spe)=$a[$i]=~/(.*)\_/;
			my $file=$spe."_100_single_copy.fa";
			open SPE, ">>$file" or die "There is no $file\n";
			my $seq=$hash{$gene};
			my $header=$orth."-".$gene;
			print SPE ">$header\n$seq\n";
		}
	}
}
```

the sequences will be in *_100_single_copy.fa   
```bash
perl temp7.pl
```
### run in SNORLAX
working dir: ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture    
```bash
mkdir result
cp reblast.pl result/
cp *_100_single_copy.fa result/
cd result/
makeblastdb -in Acura_100_single_copy.fa -out Acura -dbtype nucl
makeblastdb -in Apoly_100_single_copy.fa -out Apoly -dbtype nucl
makeblastdb -in Daru_100_single_copy.fa -out Daru -dbtype nucl
```

```perl
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
```
**Start**    
```bash
nohup perl gene_capture.pl >gene_capture.process 2>&1 &
```
***
**parallel running**   
#### 1. rmrep.pl   
vi gene_capture_parallel.pl   
```perl
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
        my $cmd1="./rmrep.pl -taxalist=$ind";
        my $ind_results=$ind."_results";
#        my $cmd2="./bandp.pl -query=$bait -subject=$ind";
        if (-d "$ind_results") {
                next;
        } else {
                push @cmds, $cmd1;
#                system($cmd2);
        }
}

my $manager = new Parallel::ForkManager(5);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
```

```bash
nohup perl gene_capture_parallel.pl >gene_capture.process 2>&1 &
# [1] 21756
```
***
#### 2. bandp.pl   
vi run_bandp.pl   
```perl
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
```

```bash
nohup perl run_bandp.pl >run_bandp.process 2>&1 &
# [2] 25984
```
##### 只选择50个single copy genes
temp2.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash; my $i;
open FIL, "$ARGV[0]" or die "can not open $ARGV[0]\n";
while (<FIL>) {
	chomp;
	if (/>/) {
		$i++;
		print "$_\n" if $i<=50;
	} else {
		print "$_\n" if $i<=50;
	}
}
```

```bash
perl temp2.pl Daru_100_single_copy.fa >Daru_50_single_copy.fa
perl temp2.pl Ocomp_100_single_copy.fa >Ocomp_50_single_copy.fa
perl temp2.pl Padel_100_single_copy.fa >Padel_50_single_copy.fa
perl temp2.pl Pmol_100_single_copy.fa >Pmol_50_single_copy.fa
```

vi run_bandp1.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @fq=<*_1.fastq.gz>;
my @cmds;
foreach my $fq (@fq) {
        (my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
        (my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
        my $bait=$spe."_50_single_copy";
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
```

```bash
nohup perl run_bandp1.pl >run_bandp1.process 2>&1 &
# [1] 11347
```

***
#### 3. Trinity.pl (my own workstation)
Acura && Apoly   
```bash
nohup Trinity.pl -input=. -output=trinity >trinity.process 2>&1 &
# [1] 17900
# Daru
nohup Trinity.pl -input=. -output=trinity >trinity.process 2>&1 &
# [1] 5726
# Padel
nohup Trinity.pl -input=. -output=trinity >trinity.process 2>&1 &
# [1] 12675
# Pmol
nohup Trinity.pl -input=. -output=trinity >trinity.process 2>&1 &
# [1] 28220
```
***
### 4. getbest.pl
```bash
mkdir result
ll -d Acura*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
# Acura
nohup getbest.pl -query="Acura_100_single_copy" -subject="Acura10 Acura11 Acura1 Acura21 Acura22 Acura24 Acura25 Acura26 Acura27 Acura28 Acura29 Acura30 Acura3 Acura4 Acura5 Acura6 Acura7 Acura8 Acura9" >Acura_getbest.process 2>&1 &
# [2] 28612

# Apoly
ll -d Apoly*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
nohup getbest.pl -query="Apoly_100_single_copy" -subject="Apoly10 Apoly21 Apoly23 Apoly24 Apoly26 Apoly27 Apoly28 Apoly29 Apoly2 Apoly30 Apoly31 Apoly34 Apoly35 Apoly3 Apoly41 Apoly43 Apoly44 Apoly45 Apoly46 Apoly47 Apoly48 Apoly4 Apoly5 Apoly6 Apoly7 Apoly8 Apoly9" >Apoly_getbest.process 2>&1 &
# [3] 30501
```
#### parallel run getbest_new.pl
vi run_getbest.pl   
```perl
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

my $manager = new Parallel::ForkManager(5);
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
		my $cmd="cat $resultdirnf/$gene*.fas >$gene.fas";
		system($cmd);
		system("rm $resultdirnf/$gene*.*.fas");
	}
}
```

```bash
nohup perl run_getbest.pl -query="Apoly_100_single_copy" -subject="Apoly10 Apoly21 Apoly23 Apoly24 Apoly26 Apoly27 Apoly28 Apoly29 Apoly2 Apoly30 Apoly31 Apoly34 Apoly35 Apoly3 Apoly41 Apoly43 Apoly44 Apoly45 Apoly46 Apoly47 Apoly48 Apoly4 Apoly5 Apoly6 Apoly7 Apoly8 Apoly9" >Apoly_run_getbest.process 2>&1 & 
# [1] 4356
```

### 5. reblast.pl
```bash
cd result
makeblastdb -in Acura_100_single_copy.fa -out Acura -dbtype nucl
makeblastdb -in Apoly_100_single_copy.fa -out Apoly -dbtype nucl

```
***
### for the rest species (Ocomp, Padel, Pmol)
#### 1. rmrep.pl && bandp.pl
```perl
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
        my $cmd1="./rmrep.pl -taxalist=$ind";
        my $ind_results=$ind."_results";
#        my $cmd2="./bandp.pl -query=$bait -subject=$ind";
        if (-d "$ind_results") {
                next;
        } else {
                push @cmds, $cmd1;
#                system($cmd2);
        }
}

my $manager = new Parallel::ForkManager(7);
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;

foreach my $fq (@fq) {
        (my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
        (my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
        my $bait=$spe."_50_single_copy";
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
```

```bash
nohup perl run_rmrep_bandp.pl > run_rmrep_bandp.process 2>&1 &
# [1] 17581
```

run bandp.pl parallel   
vi run_bandp_parallel.pl (change the blastn parameter "-num_threads 24" to "-num_threads 8 (20)")   
```perl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @fq=<*_1.fastq.gz>;
my @cmds;
foreach my $fq (@fq) {
        (my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
        (my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
        my $bait=$spe."_50_single_copy";
#        my $cmd1="./rmrep.pl -taxalist=$ind";
        my $ind_results=$ind."_results";
        my $cmd2="./bandp.pl -query=$bait -subject=$ind";
        if (-d "$ind_results") {
                next;
        } else {
                push @cmds, $cmd2;
#                system($cmd2);
        }
}

my $manager = new Parallel::ForkManager(3); # 3 -> 2
foreach my $cmd (@cmds) {
        $manager->start and next;
        system($cmd);
        $manager->finish;
}
$manager -> wait_all_children;
```

```bash
nohup perl run_bandp_parallel.pl > un_bandp_parallel.process 2>&1 &
# [1] 585
```

#### 2. Trinity.pl
```bash
nohup ./Trinity_new.pl -species=Acura -output=trinity >Acura_trinity.process 2>&1 &
# [1] 18075
nohup ./Trinity_new.pl -species=Apoly -output=trinity >Apoly_trinity.process 2>&1 &
# [2] 22463
nohup ./Trinity_new.pl -species=Daru -output=trinity >Daru_trinity.process 2>&1 &
# [3] 23105
nohup ./Trinity_new.pl -species=Ocomp -output=trinity >Ocomp_trinity.process 2>&1 &
# [1] 16431
nohup ./Trinity_new.pl -species=Padel -output=trinity >Padel_trinity.process 2>&1 &
# [5] 26013
nohup ./Trinity_new.pl -species=Pmol -output=trinity >Pmol_trinity.process 2>&1 &
# [6] 28478
```

### 3. run_getbest.pl
```bash
# Apoly
nohup perl run_getbest.pl -query="Apoly_100_single_copy" -subject="Apoly10 Apoly21 Apoly23 Apoly24 Apoly26 Apoly27 Apoly28 Apoly29 Apoly2 Apoly30 Apoly31 Apoly34 Apoly35 Apoly3 Apoly41 Apoly43 Apoly44 Apoly45 Apoly46 Apoly47 Apoly48 Apoly4 Apoly5 Apoly6 Apoly7 Apoly8 Apoly9" >Apoly_run_getbest.process 2>&1 &
# [1] 27287

# Acura
nohup getbest.pl -query="Acura_50_single_copy" -subject="Acura10 Acura11 Acura1 Acura21 Acura22 Acura24 Acura25 Acura26 Acura27 Acura28 Acura29 Acura30 Acura3 Acura4 Acura5 Acura6 Acura7 Acura8 Acura9" >Acura_getbest.process 2>&1 &

# Daru
perl run_getbest.pl -query="Daru_50_single_copy" -subject="Daru10 Daru11 Daru12 Daru1 Daru21 Daru23 Daru24 Daru25 Daru26 Daru27 Daru29 Daru2 Daru30 Daru32 Daru33 Daru3 Daru4 Daru6 Daru7 Daru9"

# Ocomp
perl run_getbest.pl -query="Ocomp_50_single_copy" -subject="Ocomp1 Ocomp21 Ocomp22 Ocomp23 Ocomp24 Ocomp25 Ocomp26 Ocomp27 Ocomp28 Ocomp29 Ocomp2 Ocomp30 Ocomp3 Ocomp4 Ocomp5 Ocomp6 Ocomp7 Ocomp8 Ocomp9"

# Padel
perl run_getbest.pl -query="Padel_50_single_copy" -subject="Padel10 Padel11 Padel12 Padel15 Padel16 Padel17 Padel18 Padel21 Padel22 Padel23 Padel24 Padel25 Padel26 Padel28 Padel29 Padel30 Padel31 Padel32 Padel6 Padel7 Padel8 Padel9"

# Pmol
perl run_getbest.pl -query="Pmol_50_single_copy" -subject="Pmol10 Pmol11 Pmol13 Pmol21 Pmol22 Pmol23 Pmol24 Pmol25 Pmol26 Pmol27 Pmol28 Pmol29 Pmol31 Pmol32 Pmol41 Pmol42 Pmol43 Pmol44 Pmol45 Pmol46 Pmol7 Pmol8 Pmol9"
```

vi run_getbest.sh   
```bash
perl run_getbest.pl -query="Daru_50_single_copy" -subject="Daru10 Daru11 Daru12 Daru1 Daru21 Daru23 Daru24 Daru25 Daru26 Daru27 Daru29 Daru2 Daru30 Daru32 Daru33 Daru3 Daru4 Daru6 Daru7 Daru9"
perl run_getbest.pl -query="Ocomp_50_single_copy" -subject="Ocomp1 Ocomp21 Ocomp22 Ocomp23 Ocomp24 Ocomp25 Ocomp26 Ocomp27 Ocomp28 Ocomp29 Ocomp2 Ocomp30 Ocomp3 Ocomp4 Ocomp5 Ocomp6 Ocomp7 Ocomp8 Ocomp9"
perl run_getbest.pl -query="Padel_50_single_copy" -subject="Padel10 Padel11 Padel12 Padel15 Padel16 Padel17 Padel18 Padel21 Padel22 Padel23 Padel24 Padel25 Padel26 Padel28 Padel29 Padel30 Padel31 Padel32 Padel6 Padel7 Padel8 Padel9"
perl run_getbest.pl -query="Pmol_50_single_copy" -subject="Pmol10 Pmol11 Pmol13 Pmol21 Pmol22 Pmol23 Pmol24 Pmol25 Pmol26 Pmol27 Pmol28 Pmol29 Pmol31 Pmol32 Pmol41 Pmol42 Pmol43 Pmol44 Pmol45 Pmol46 Pmol7 Pmol8 Pmol9"
```

```bash
nohup sh run_getbest.sh > getbest.process 2>&1 &
# [1] 542
```
***
#### 4. reblast.pl
```bash
cd result
```

vi run_reblast.sh   
```bash
perl reblast.pl -query Apoly_100_single_copy.resultnf_new -database Apoly
perl reblast.pl -query Acura_100_single_copy.resultnf -database Acura
perl reblast.pl -query Daru_50_single_copy.resultnf_new -database Daru
perl reblast.pl -query Ocomp_50_single_copy.resultnf_new -database Ocomp
perl reblast.pl -query Pmol_50_single_copy.resultnf_new -database Pmol
perl reblast.pl -query Padel_50_single_copy.resultnf_new -database Padel
```

```bash
nohup sh run_reblast.sh > run_reblast.process 2>&1 &
# [1] 6803
```

the same orthogroup gene of all individuals per species copy to one fasta file   
```bash
ll *.fas|perl -alne 'if ($F[4]==0){`rm $F[-1]`}' # delete the empty file
ll *.fas|perl -alne '(my $na)=$F[-1]=~/(.*)-/;$hash{$na}++;push @na, $na if $hash{$na}==1;END{foreach my $na (@na){print"$na\t$hash{$na}"}}' # just select the genes were captured by all individuals
# OG0008841 was not captured in Apoly individuals, then now have 49 single-copy genes
# cat all sequences with the same orthogroup id together 
ll *.fas|perl -alne '(my $na)=$F[-1]=~/(.*)-/;$hash{$na}++;push @na, $na if $hash{$na}==1;END{foreach my $na (@na){print"$na\t$hash{$na}"}}'|perl -alne '`cat $F[0]*.fas > $F[0].fas`'
```

align, trim, concatenate   
vi temp1.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

# muscle -in OG0008797.fas -out OG0008797-align.fas
# trimal -in OG0008797-align.fas -out OG0008797-align-trim.fas -gt 0.9 -st 0.001 -cons 60

my @fasta=<*.fas>;
foreach my $fasta (@fasta) {
	(my $name)=$fasta=~/(.*)\.fas/;
	my $align=$name."-align.fas";
	my $trim=$name."-trim.fas";
	my $cmd1="muscle -in $fasta -out $align";
	my $cmd2="trimal -in $align -out $trim -gt 0.9 -st 0.001 -cons 60";
	system($cmd1);
	system($cmd2);
}

my %hash; 
my $gene; my @genes;
my $i;
my @trims=<*-trim.fas>;
foreach my $trim (@trims) {
	$i++;
	open TRIM, "$trim";
	while (<TRIM>) {
		chomp;
		if (/>/) {
			s/>//;
			$gene=$_;
			push @genes, $gene if $i==1;
		} else {
			$hash{$gene}.=$_;
		}
	}
}

open FASTA, ">All_gene_concatenated.fasta";
foreach my $gen (@genes) {
	my $seq=$hash{$gen};
	print FASTA ">$gen\n$seq\n";
}
```

extract the longest transcript of OG0008841 in Apoly, and then align   
vi temp1.pl (~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/trinity)    
```perl
#!/usr/bin/perl
use strict;
use warnings;

# extract the longest transcript of OG0008841 in Apoly
open FASTA, ">OG0008841-Apoly_10124.fasta" or die "can not create OG0008841-Apoly_10124.fasta\n";
my @dirs=<Apoly*>;
foreach my $dir (@dirs) {
	my $ind=$dir;
	my $fas="$dir/OG0008841-Apoly_10124.fasta";
	my %hash; my $len;
	open FAS, "$fas" or die "can not open $fas\n";
	while (<FAS>) {
		chomp;
		my @a=split;
		if (/>/) {
			s/>//;
			($len)=$a[1]=~/len=(\d+)/;
		} elsif ($hash{$ind}) {
			my $len_old=$hash{$ind}->{'len'};
			if ($len>$len_old) {
				$hash{$ind}={
					'len' => $len,
					'seq' => $_
				};
			}
		} else {
			$hash{$ind}={
				'len' => $len,
				'seq' => $_
			};
		}
	}
	my $seq=$hash{$ind}->{'seq'};
	print FASTA ">$ind\n$seq\n";
}
```

```bash
perl temp1.pl
cp OG0008841-Apoly_10124.fasta ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result
cd ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result
cat OG0008841-Apoly_10124.fasta Acura_100_single_copy.resultnf/OG0008841-Acura_177266.fas Daru_50_single_copy.resultnf_new/OG0008841-Daru_163989.fas Ocomp_50_single_copy.resultnf_new/OG0008841-Ocomp_153633.fas Padel_50_single_copy.resultnf_new/OG0008841-Padel_212527.fas Pmol_50_single_copy.resultnf_new/OG0008841-Pmol_87384.fas >OG0008841.fasta
muscle -in OG0008841.fasta -out OG0008841-align.fas
trimal -in OG0008841-align.fas -out OG0008841-trim.fas -gt 0.9 -st 0.001 -cons 60
cd genebin
mv ../OG0008841*.fas ./
perl temp1.pl
fasta2phy.pl All_gene_concatenated.fasta>All_gene_concatenated.phy

# remove OG0008841 OG0008839
mv OG0008820* OG0008839* ../
perl temp1.pl
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated.phy
# raxml
nohup raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m GTRCAT -T 24 -s All_gene_concatenated.phy -n All_gene_concatenated > Raxml.process 2>&1 &
[1] 9700
```

***
## run in my own workstation
working dir: ~/Desktop/PapueNewGuinea-new/merge_clean   
```bash
mkdir gene_capture
mv Ocomp*.fastq.gz Padel*.fastq.gz Pmol*.fastq.gz gene_capture/
cd gene_capture
scp kang1234@147.8.76.155:~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/*.pl ./
scp kang1234@147.8.76.155:~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/*100_single_copy.fa ./
scp kang1234@147.8.76.155:~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/*.pm ./
mkdir result
cp reblast.pl result/
cp *_100_single_copy.fa result/
cd result/
makeblastdb -in Ocomp_100_single_copy.fa -out Ocomp -dbtype nucl
makeblastdb -in Padel_100_single_copy.fa -out Padel -dbtype nucl
makeblastdb -in Pmol_100_single_copy.fa -out Pmol -dbtype nucl
cd ../
nohup perl gene_capture.pl >gene_capture.process 2>&1 &
# [1] 30780
```
RAM is not enough, just do the small size fastq files at present   
vi gene_capture_small.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open FIL, "large_ind.txt" or die "can not open large_ind.txt\n";
while (<FIL>) {
	chomp;
	$hash{$_}++;
}
my @fq=<*_1.fastq.gz>;
foreach my $fq (@fq) {
        (my $ind)=$fq=~/(.*)_1\.fastq\.gz/;
        (my $spe)=$fq=~/(\D+)\d+_1\.fastq\.gz/;
        my $bait=$spe."_100_single_copy";
        my $cmd1="./rmrep.pl -taxalist=$ind";
        my $cmd2="./bandp.pl -query=$bait -subject=$ind";
        if (-e "$ind.fas") {
                next;
        } elsif ($hash{$ind}) {
        	next;
        } else {
                system($cmd1);
                system($cmd2);
        }
}
```

```bash
ll -h *.fastq.gz|perl -alne '(my $ind)=$F[-1]=~/(.*)_\d\.fastq\.gz/;(my $size)=$F[4]=~/(.*)G/;print $ind if $size>=1.8'|sort -u >large_ind.txt
nohup perl gene_capture_small.pl >gene_capture.process 2>&1 &
# [1] 8763
```

get the longest transcript of each gene of all individuals per species   
vi ORFextract   
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long 'HelpMessage';

my $usage=<<_EOH_;;
------------------------------------------------------------------------------------------
run "TransDecoder.LongOrfs" to each fasta file in directory

Usage: 
ORFextract --dir Acura_100_single_copy.resultnf --out Acura_orf

													Kang 2021-8-24
------------------------------------------------------------------------------------------
_EOH_
;

GetOptions(
	'dir:s', \my $dir,	# the directory including fasta files to detect ORF
	'out:s', \ my $out,	# the output directory
	'help', \ my $help
	);

if ($help || (! $dir) || (! $out) ) {
	die $usage; # all of these options are mandatory requirements
}

unless (-d $out) {
	system("mkdir $out");
}
# TransDecoder.LongOrfs -t OG0008791-Acura_1.fas -O OG0008791-Acura_1
my @fasta=<$dir/*.fas>;
foreach my $fas (@fasta) {
#	open FAS, $fas or die "can not open $fas\n";
	(my $gene)=basename($fas)=~/(.*)\.fas/;
	my $outdir=$gene;
	my $cmd1="TransDecoder.LongOrfs -t $fas -O $gene";
	system($cmd1);
	my $flag=&extract_longest_transcript($gene);
	if ($flag) {
		system("rm -rf $gene");
		system("rm -rf $gene.__checkpoints_longorfs");
		system("rm pipeliner.*.cmds");
	}
}

################################
sub extract_longest_transcript {
	my ($gene)=@_;
	my $flag;
	my $orf=$gene."/longest_orfs.cds";
	my $longorf=$gene.".fas";
	open LONGORF, ">$longorf" or die "can not create $longorf\n";
	my ($gene1, $len); my @genes;
	my %hash;
	if (-e $orf) {
		open ORF, $orf or die "can not open $orf\n"; # extract the longest ORF per individual
		while (<ORF>) {
			chomp;
			my @a=split;
			if (/>/) {
				($gene1)=$a[-1]=~/(.*)\:/;
				($len)=$a[-2]=~/len\:(.*)/;
			} elsif ($hash{$gene1}) {
				my $len_old=$hash{$gene1}->{'len'};
				my $seq=$_;
				if ($len > $len_old) {
					$hash{$gene1}={
						'len' => $len,
						'seq' => $seq
					};
				}
			} else {
				my $seq=$_;
				$hash{$gene1}={
					'len' => $len,
					'seq' => $seq
				};
			}
		}
		foreach my $gene2 (sort keys %hash) {
			my $seq=$hash{$gene2}->{'seq'};
			print LONGORF ">$gene2\n$seq\n";
		}
		close LONGORF;
		system("mv $longorf $out/");
		$flag=1;
		return($flag);
	} else {
		$flag=0;
		return($flag);
	}
}
```

```bash
perl ORFextract --dir Acura_100_single_copy.resultnf --out Acura_orf
perl ORFextract --dir Apoly_100_single_copy.resultnf_new --out Apoly_orf
perl ORFextract --dir Daru_50_single_copy.resultnf_new --out Daru_orf
perl ORFextract --dir Ocomp_50_single_copy.resultnf_new --out Ocomp_orf
perl ORFextract --dir Padel_50_single_copy.resultnf_new --out Padel_orf
perl ORFextract --dir Pmol_50_single_copy.resultnf_new --out Pmol_orf
```

cat sequences of all species according orthogroup id   
```bash
# Kang@fishlab3 Fri Dec 10 21:54:56 ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result
ll Ocomp_orf/*.fas|perl -alne '(my $n)=$F[-1]=~/\/(.*)-/;print $n' >1.tx
```

vi temp1.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my @spes=qw(Daru Ocomp Pmol Padel Acura Apoly);
open TXT, "1.txt" or die "can not open 1.txt";
while (<TXT>) {
	chomp;
	my $cmd2="cat ";
	my $gene=$_;
	foreach my $spe (@spes) {
		my $dir=$spe."_orf";
		my $gene1=$dir."/".$gene."\*".".fas";
		$cmd2.=$gene1." ";
	}
	$cmd2.=">$gene.fas";
	`mkdir Total_orf` unless -d "Total_orf";
	system($cmd2);
	system("mv $gene.fas Total_orf/");
}
```

```bash
perl temp1.pl
```
the results in (~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf)    
only select the individuals that were detected ORF in all orthogroup id   
```bash
for n in *.fas;do grep '>' ${n};done|perl -alne 's/>//;$hash{$_}++;END{foreach my $key (keys %hash) {print $key if $hash{$key}==50}}' >ind_id.txt
```

vi temp1.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my @fas=<*.fas>;
my %hash1;
open TXT, "ind_id.txt" or die "can not open 1.txt";
while (<TXT>) {
	chomp;
	$hash1{$_}++;
}
close TXT;

foreach my $fas (@fas) {
	my $ind; my %hash2;
	open FAS, $fas or die "can not open $fas\n";
	while (<FAS>) {
		chomp;
		if (/>/) {
			s/>//;
			$ind=$_;
		} else {
			$hash2{$ind}.=$_;
		}
	}
	system("mkdir Final_ind") unless -d "Final_ind";
	open FASTA, ">Final_ind/$fas";
	foreach my $key (sort keys %hash2) {
		if ($hash1{$key}) {
			my $seq=$hash2{$key};
			print FASTA ">$key\n$seq\n";
		}
	}
	close FASTA;
}
```

```bash
perl temp1.pl
```

the final sequences were in ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind    
align and then trim   
```bash
perl temp2.pl
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated.phy
# run raxml in SNORLAX
# Kang@fishlab3 Fri Dec 10 23:15:38 ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind
scp All_gene_concatenated.fasta All_gene_concatenated.phy kang1234@147.8.76.155:~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result/Final_ind
```
No Apoly individuals, select the genes were captured by all individuals   
```bash
# Kang@fishlab3 Fri Dec 10 23:35:49 ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf
mkdir Final_ind_2
# the gene should be captured by all individuals
for n in *.fas;do echo -n "${n}  ";grep '>' ${n}|wc -l;done|perl -alne '`cp $F[0] Final_ind_2` if $F[1]==130'
cp ../Final_ind/temp2.pl ./
# check if any individuals were removed in the trimmed file
# Kang@fishlab3 Sat Dec 11 00:19:36 ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2
for n in *-trim.fas;do echo -n "${n} "; grep '>' ${n}|wc -l;done|perl -alne 'print if $F[1] != 130'
# OG0008820-trim.fas 129
# OG0008843-trim.fas 129
# still have 31 genes to be concatenated
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated.phy
scp All_gene_concatenated.fasta All_gene_concatenated.phy kang1234@147.8.76.155:~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result/Final_ind
# (base) kang1234@celia-PowerEdge-T640 Sat Dec 11 00:26:17 ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result/Final_ind
# 45343 bp
nohup raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m GTRCAT -T 24 -s All_gene_concatenated.phy -n All_gene_concatenated > Raxml.process 2>&1 &
# [1] 16227

# set all Ocomp individuals as root
nohup raxmlHPC -f a -x 12345 -p 12345 -# 1000 -m GTRCAT -T 24 -s All_gene_concatenated.phy -o Ocomp21,Ocomp1,Ocomp5,Ocomp23,Ocomp28,Ocomp25,Ocomp9,Ocomp6,Ocomp2,Ocomp27,Ocomp30,Ocomp8,Ocomp22,Ocomp26,Ocomp24,Ocomp29,Ocomp3,Ocomp7,Ocomp4 -n All_gene_concatenated > Raxml.process 2>&1 &
# [1] 25108
```
***
based on RAxML_bestTree.All_gene_concatenated for a time-calibrated ultrametric phylogenetic tree
```bash
# Kang@fishlab3 Sat Dec 11 10:22:56 ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm
cp /media/HDD/cleaner_fish/genome/gene_family/cafetutorial_prep_r8s.py ./
cp ../RAxML_bestTree.All_gene_concatenated ./
python cafetutorial_prep_r8s.py -i RAxML_bestTree.All_gene_concatenated -o r8s_ctl_file.txt -s 638853 -p 'Ocomp1,Daru1' -c '109.38'
r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt

less r8s_ultrametric.txt|perl -alne 's/omparu//g;print'
```

```bash
python cafetutorial_prep_r8s.py -i RAxML_bestTree.All_gene_concatenated -o r8s_ctl_file_1.txt -s 638853 -p 'Ocomp21,Ocomp1,Ocomp5,Ocomp23,Ocomp28,Ocomp25,Ocomp9,Ocomp6,Ocomp2,Ocomp27,Ocomp30,Ocomp8,Ocomp22,Ocomp26,Ocomp24,Ocomp29,Ocomp3,Ocomp7,Ocomp4,Daru10,Daru11,Daru12,Daru1,Daru21,Daru23,Daru24,Daru25,Daru26,Daru27,Daru29,Daru2,Daru30,Daru32,Daru33,Daru3,Daru4,Daru6,Daru7,Daru9' -c '109.38'
r8s -b -f r8s_ctl_file_1.txt > r8s_tmp_1.txt
tail -n 1 r8s_tmp_1.txt | cut -c 16- > r8s_ultrametric.txt
less r8s_ultrametric.txt|perl -alne 's/0\.000000/0\.000001/g;s/p21mp1//ig;print' >r8s_ultrametric_1.txt
```

### bamm
working dir: "~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm"   
```bash
perl temp1.pl OG0015521 >OG0015521.txt # S4A4
perl temp1.pl OG0033097 >OG0033097.txt # CAH10
perl temp1.pl OG0020241 >OG0020241.txt # CREB1
mkdir OG0020241_trait
mv OG0020241.txt OG0020241_trait/
cd OG0020241_trait/
cp ../OG0033097_trait/r8s_ultrametric_1.txt ./
cp ../OG0033097_trait/traitcontrol.txt ./

perl temp1.pl OG0011939 > OG0011939.txt # S4A8
mkdir OG0011939_trait
mv OG0011939.txt OG0011939_trait/
cd OG0011939_trait/
cp ../OG0033097_trait/r8s_ultrametric_1.txt ./
cp ../OG0033097_trait/traitcontrol.txt ./
bamm -c traitcontrol.txt

perl temp1.pl OG0002874 > OG0002874.txt # CLOCK
mkdir OG0002874_trait
mv OG0002874.txt OG0002874_trait/
cd G0002874_trait/
cp ../OG0033097_trait/r8s_ultrametric_1.txt ./
cp ../OG0033097_trait/traitcontrol.txt ./
bamm -c traitcontrol.txt
```

Ion transport   
OG0007437: CAC1H   
OG0004968: CACB1   
OG0020625: KCNH4   
OG0022213: KCNK9   

vi temp2.pl   
```perl
#!/usr/bin/perl -w
use strict;
use warnings;

# to prepare files and mkdir dir for bamm
my $gene=$ARGV[0];

my $cmd1="perl temp1.pl $gene > $gene.txt"; # create the expression data for the gene

my $dir=$gene."_trait";
my $cmd2="mkdir $dir";
my $cmd3="mv $gene.txt $dir"; # move the expression data to the created dir
my $cmd4="cp OG0033097_trait/r8s_ultrametric_1.txt $dir"; # copy the tree to the dir

system($cmd1);
system($cmd2) unless -d $dir;
system($cmd3);
system($cmd4);

# create the control file in the dir
my $cont_f="OG0033097_trait/traitcontrol.txt";
my $cont_n="$dir/traitcontrol.txt";
open FIL1, $cont_f or die "can not open $cont_f\n";
open FIL2, ">$cont_n" or die "can not create $cont_n\n";

while (<FIL1>) {
	chomp;
	if (/traitfile\s=\sOG.*/) {
		print FIL2 "traitfile = $gene.txt\n";
	} else {
		print FIL2 "$_\n";
	}
}
```

```bash
perl temp2.pl OG0007437 # CAC1H
cd OG0007437_trait/
bamm -c traitcontrol.txt

perl temp2.pl OG0004968 # CACB1
cd OG0004968_trait/
bamm -c traitcontrol.txt

perl temp2.pl OG0020625 # KCNH4
cd OG0020625_trait/
bamm -c traitcontrol.txt

perl temp2.pl OG0022213 # KCNK9
cd OG0022213_trait/
bamm -c traitcontrol.txt

# pick a gene that is not positively selected
# OG0005247: BMAL2
perl temp2.pl OG0005247
cd OG0005247_trait/
bamm -c traitcontrol.txt
```
***
### prepare the table of final 32 single copy genes
working dir: ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/OrthoFinder/Results_Nov30   
vi temp8.pl   
```perl
#!/usr/bin/perl
use strict;
use warnings;

my %hash1;
open FIL1, "Final_select_32_orth.txt" or die "can not open Final_select_32_orth.txt\n";
while (<FIL1>) {
	chomp;
	$hash1{$_}++;
}
open FIL2, "paml_input/ortho_list.txt" or die "can not open paml_input/ortho_list.txt\n";
while (<FIL2>) {
	chomp;
	my @a=split /\t/;
	if ($hash1{$a[0]}) {
		(my $gene)=$a[1]=~/Zebrafish_(.*)/;
		print "$gene\n";
	}
}
```

```bash
perl temp8.pl >Final_select_32_orth_zebrafish.txt
```
And submit to ensembl the get the gene description information    
