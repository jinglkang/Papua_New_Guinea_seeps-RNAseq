The data files informations  
===========================  
#### \*\_no_reduce_enrichment.txt  
enrichment results that didn't remove the categories (p>0.05) and didn't use "reduced to most specific terms".  
#### \*\_reduce_enrichment.txt
enrichment results that remove the categories (p>0.05) and use "reduced to most specific terms".

## 1. extract_sig_more_than_two_species.pl: extract the information (from \*\_no_reduce_enrichment.txt) of the significant functions (in at least two speices) (from \*\_reduce_enrichment.txt). 
Usage: perl extract_sig_more_than_two_species.pl    
  
## 2. extract_species-specific-sig.pl: extract the sigficant functions that were only in one species but have no any DEGs from other species.  
Usage: perl extract_species-specific-sig.pl
