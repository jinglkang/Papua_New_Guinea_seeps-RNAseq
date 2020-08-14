#!/usr/bin/perl -w
my (%hash1, %hash2, %hash1_1, %hash2_1);
my (@gene1, @gene2);

# enrichment file (co2 seep vs. control), not using "reduced to most specific", not using "filter < 0.05", include the test DEGs.
open fil1, "all_enrichment_apoly_co2_control.txt"; 
while (<fil1>) {
        chomp;
        next if /^Tags/;
        my @a=split /\t/;
        $hash1_1{$a[2]}=$a[4];
        if ($a[-1]=~/OG/) {
                my @gene=split /\;/, $a[-1];
                my $gene=\@gene;
                $hash1{$a[2]}=$gene;
        } else {
                $hash1{$a[2]}=0;
        }
}

# enrichment file (co2 seep vs. mid), not using "reduced to most specific", not using "filter < 0.05", include the test DEGs.
open fil2, "all_enrichment_apoly_co2_mid.txt";
while (<fil2>) {
        chomp;
        next if /^Tags/;
        my @a=split /\t/;
        $hash2_1{$a[2]}=$a[4];
        if ($a[-1]=~/OG/) {
                my @gene=split /\;/, $a[-1];
                my $gene=\@gene;
                $hash2{$a[2]}=$gene;
        } else {
                $hash2{$a[2]}=0;
        }
}

print "GO_name\tsig_in_high_vs.low\tgene_num_in_high_vs.low\tsig_num_in_high_vs.mid\tgene_num_in_high_vs.mid\tCommon\tGap\tType\tRatio\n";

# total_enrichment.txt: include all significant functions (name) in co2 seep vs. control and co2 seep vs. mid.
open fil3, "total_enrichment.txt";
while (<fil3>) {
        s/\r$//g;
        chomp;
        if (($hash1{$_} !=0) && ($hash2{$_} !=0)) {
                my @gene1_1=@{$hash1{$_}};
                my $num1=@gene1_1;
                my @gene2_1=@{$hash2{$_}};
                my $num2=@gene2_1;
                my %hash3=map {$_=>1} @gene1_1;
                my @common=grep {$hash3{$_}} @gene2_1;
                my $num3=@common;
                my $gap=$num2-$num1;
                my $ratio=$num3/$num1;
                $ratio=sprintf "%.2f",$ratio;
                if ($hash1_1{$_}<=0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tsig\t$num1\tsig\t$num2\t$num3\t$gap\tCommon\t$ratio\n";
                } elsif ($hash1_1{$_}>0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tunsig\t$num1\tsig\t$num2\t$num3\t$gap\tCO2_seep_vs._Mid\t$ratio\n";
                } elsif ($hash1_1{$_}<=0.05 && $hash2_1{$_}>0.05) {
                        print "$_\tsig\t$num1\tunsig\t$num2\t$num3\t$gap\tCO2_seep_vs_Control\t$ratio\n";
                }
        } elsif ($hash1{$_}==0 && ($hash2{$_} !=0)) {
                my $num1=0;
                my @gene2_1=@{$hash2{$_}};
                my $num2=@gene2_1;
                my $num3=0;
                my $gap=$num2-$num1;
                my $ratio=0;
                if ($hash1_1{$_}<=0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tsig\t$num1\tsig\t$num2\t$num3\t$gap\tCommon\t$ratio\n";
                } elsif ($hash1_1{$_}>0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tunsig\t$num1\tsig\t$num2\t$num3\t$gap\tCO2_seep_vs._Mid\t$ratio\n";
                } elsif ($hash1_1{$_}<=0.05 && $hash2_1{$_}>0.05) {
                        print "$_\tsig\t$num1\tunsig\t$num2\t$num3\t$gap\tCO2_seep_vs_Control\t$ratio\n";
                }
        } elsif (($hash1{$_} != 0) && ($hash2{$_}==0)) {
                my @gene1_1=@{$hash1{$_}};
                my $num1=@gene1_1;
                my $num2=0;
                my $num3=0;
                my $gap=$num2-$num1;
                my $ratio=0;
                if ($hash1_1{$_}<=0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tsig\t$num1\tsig\t$num2\t$num3\t$gap\tCommon\t$ratio\n";
                } elsif ($hash1_1{$_}>0.05 && $hash2_1{$_}<=0.05) {
                        print "$_\tunsig\t$num1\tsig\t$num2\t$num3\t$gap\tCO2_seep_vs._Mid\t$ratio\n";
                } elsif ($hash1_1{$_}<=0.05 && $hash2_1{$_}>0.05) {
                        print "$_\tsig\t$num1\tunsig\t$num2\t$num3\t$gap\tCO2_seep_vs_Control\t$ratio\n";
                }
        }
}
