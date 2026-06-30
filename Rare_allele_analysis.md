# Rare alleles and Fine Scale Population Genomics
This is a place for me to leave notes about how useful rare alleles can be for fine scale population genomics.  

**Background**: Previous studies highlight that rare alleles can identify hidden population structure in human genomics, while other studies reveal that standard filtering of minor alleles at low frequency (often 5%), can bias downstream analyses since a large proportion of variation is held within each population's minor allele frequencies. In terms of **coral genomics**, many studies that aim to identify thermally associated genes find minimal signals for adaptation, possibly due to high connectivity within the marine environment, which can elevate selected alleles beyond standard filtering thresholds. However, we find that rare alleles can reveal hidden population structure within a singular population of *Porites rus* at a fine spatial scale of only 40km. 

We used plink/2.0a7 to filter and evaluate the effects of different MAF filtering thresholds. Here we find 
  1) Removal of standard MAF thresholds *and* singletons **reveals hidden population structure** associated with warm thermal environments, and
  2) Standard filtering PCA results are **susceptible** to Linkage Disquilibrium (LD), while **Rare Allele PCA results** (with removed or a limited MAF threshold) are **robust** to LD. 

## Separation of datasets by MAF threshold and LD-pruning  ##

We generated different genomic datasets across different MAF thresholds including: No Minor Allele Filtering (no MAF; all rare alleles are kept), Rare Alleles (MAF < 0.05; only rare alleles at frequency less than 5% are kept), Minor Allele Counts 2-5 (MAC 2-5 [doubletons - quintupletons]; very rare alleles at frequency less than 2%), Standard Filtering at 5% (MAF > 0.05; some rare alleles filtered), and Standard Filtering at 10% (MAF > 0.10; all rare alleles filtered). For each dataset we also created thier LD pruned counterpart, to evaluate differences in structure between pruned (LD) and unpruned datasets. Singletons (MAC = 1) were excluded from all analyses. 

## Results ##

A PCA for each dataset was created, and we evaluated the presence of population structure via the overlap of 95% confident ellipses. Here we find evidence of population structure with limited MAF filtering that decreases as filtering increases (top left -> top right), and the absence of structure with standard MAF filtering thresholds (bottom).
<br><br>

<img width="6000" height="3600" alt="combined_ld_plots" src="https://github.com/user-attachments/assets/3ab180d0-bf68-475c-bf51-bc522894a312" />

#######################################################################################################
<br><br>

Additonally, side-by-side comparisons between unpruned (left) and pruned (right) datasets shows that **No MAF Filtering is robust to LD-pruning** while Standard Filtering is less so. 
<br><br>
<img width="3000" height="3000" alt="standard_combined_minor_allele_pca_4panel" src="https://github.com/user-attachments/assets/da5ff5df-231a-4af1-bedd-802e37008c89" />

#######################################################################################################
<br><br>

This holds true with the Rare Allele dataset (where only alleles at less than 5% frequency were kept), however we do report that allele counts alone (MAC 2-5) is **susceptible** to LD pruning differences, possibly reflecting highly related individuals that share many quintupletons specific to their familial group (unvalidated).
<br><br>
<img width="2400" height="1800" alt="combined_minor_allele_pca_4panel" src="https://github.com/user-attachments/assets/9d2d187e-543c-45ef-adc8-38490e0e48f1" />
