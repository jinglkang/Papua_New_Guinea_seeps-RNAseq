@file=<*_no_reduce_enrichment.txt>;
foreach $file (@file) {
        $i++;
        open FI, $file;
        ($spe)=$file=~/(.*)\_no\_reduce\_enrichment\.txt/ if $file!=~/_/;
        $spe="Apoly" if $file=~/_/;
        @header=split /\t/, <FI>;
        @header=@header[1..9];
        print join ("\t", @header) if $i==1;
        print "\tRatio\tSpecies\n" if $i==1;
        while (<FI>) {
                chomp;
                next if /^Tags/;
                s/OG\d+.*;//g;
                s/OG\d+//g;
                s/\s+$//g;
                s/^\[.*\]\s+//;
                @a=split /\t/;
                $info=join "\t", @a;
                $hash1{$a[1]}++;
                $hash2{$a[1]}=$a[0];
                push @term, $a[1] if $hash1{$a[1]}==1;
                push @{$a[1]}, $a[3];
                $hash4{$spe}->{$a[1]}=$info;
        }
}

@txt=<*_reduce_enrichment.txt>;
foreach $txt (@txt) {
        open TXT, $txt;
        while (<TXT>) {
                chomp;
                next if /^Tags/;
                @a=split /\t/;
                $hash3{$a[2]}++;
        }
}

@spe=("Apoly", "Daru", "Ocomp", "Padel", "Pmol");
foreach $term (@term) {
        @fdr=@{$term};
        @fdr=sort {$a<=>$b} @fdr;
        if (@fdr==5 && $fdr[0]<=0.05 && $fdr[1]<=0.05) {
                $num=@fdr;
                foreach $spe (@spe) {
                        $info=$hash4{$spe}->{$term};
                        @a=split /\t/, $info;
                        $ratio=$a[-4]/($a[-2]+$a[-4]);
                        print "$info\t$ratio\t$spe\n" if $hash3{$term};
                }
        }
}
