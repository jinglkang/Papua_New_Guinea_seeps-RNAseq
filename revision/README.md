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
