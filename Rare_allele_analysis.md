# Rare alleles and Fine Scale Population Genomics
This is a place for me to leave notes about how useful rare alleles can be for fine scale population genomics. 

**Background**: Previous studies highlight that rare alleles can identify hidden population structure in human genomics, while other studies reveal that standard filtering of minor alleles at low frequency (often 5%), can bias downstream analyses since a large proportion of variation is held within each population's minor allele frequencies. In terms of **coral genomics**, many studies that aim to identify thermally associated genes find minimal signals for adaptation, possibly due to high connectivity within the marine environment, which can elevate selected alleles beyond standard filtering thresholds. However, we find that rare alleles can reveal hidden population structure within a singular population of *Porites rus* at a fine spatial scale of only 40km. 

We used plink/2.0a7 to filter and evaluate the effects of different MAF filtering thresholds and find 
  1) Removal of standard MAF thresholds *and* singletons **reveals hidden population structure associated with warm thermal environments**, and
  2) Standard filtering PCA results are **susceptible** to Linkage Disquilibrium (LD), while **Rare Allele PCA results** (with removed or a limited MAF threshold) are **robust** to LD. 
