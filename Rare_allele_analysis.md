# Rare alleles and Fine Scale Population Genomics
This is a place for me to leave notes about how useful rare alleles can be for fine scale population genomics. 

**Background**: Previous studies highlight that rare alleles can identify hidden population structure in human genomics, while other studies reveal that standard filtering of minor alleles at low frequency (often 5%), can bias downstream analyses since a large proportion of variation is held within each population's minor allele frequencies. In terms of **coral genomics**, many studies that aim to identify thermally associated genes find minimal signals for adaptation, possibly due to high connectivity within the marine environment, which can elevate selected alleles beyond standard filtering thresholds. However, we find that rare alleles can reveal hidden population structure within a singular population of *Porites rus* at a fine spatial scale of only 40km. 

We used plink/2.0a7 to filter and evaluate the effects of different MAF filtering thresholds. Here we find 
  1) Removal of standard MAF thresholds *and* singletons **reveals hidden population structure** associated with warm thermal environments, and
  2) Standard filtering PCA results are **susceptible** to Linkage Disquilibrium (LD), while **Rare Allele PCA results** (with removed or a limited MAF threshold) are **robust** to LD. 

## Methods ##

We generated different genomic datasets across different MAF thresholds including: No Minor Allele Filtering (no MAF), Minor Allele Counts 2-5 (MAC 2-5; doubletons - quintupletons), Rare Alleles (MAF < 0.05), Standard Filtering at 5% (MAF > 0.05), and Standard Filtering at 10% (MAF > 0.10). For each dataset we also created thier LD pruned counterpart, to evaluate differences in structure between pruned (LD) and unpruned datasets. Singletons (MAC = 1) were excluded from all analyses. 

A PCA for each dataset was created, and we evaluated the presence of population structure via the overlap of 95% confident ellipses. 
