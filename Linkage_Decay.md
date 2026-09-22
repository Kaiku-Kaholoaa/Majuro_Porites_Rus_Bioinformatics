# Linkage Decay Analysis

The purpose of this analysis is to identify the optimal distance for LD pruning. 

Conceptually, the closer snps are to each other in the genome, the higher the chance they are inherited together. This analysis helps us identify at what distance do these average corrleations drop, and more importantly plateau. The plateau will indicate that at x distance or greater, the average correlation among linked genes is even (and thus negliable).

In this doc, we will use plink to gather the distances and correlations between snps, and then we will use our own python script to bin distance sizes and calculate average correlations (r^2) for each bin range. Finally, we will plot our results to identify the optimal distance treshold for ld pruning. Knowing this distance will be essential for the doubleton analysis, for which we'll include minor alleles (excluded in this analysis), but will still need to properly filter for linkage. In that analysis, we'll call doubletons and ensure that no doubletons are within x distance of each other.

## Step 1: Use plink to obtain our .vcor file:
Starting with our fully qa/qc'ed dataset (prus_qc_noclones), lets use plink to identify correlated snps. 

`cat plink_linkage_decay.sbatch `

```bash
#!/bin/bash
#SBATCH --job-name=prus_pgen_ld
#SBATCH -p serc,hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=600G
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

ml plink/2.0a7

plink2 --pfile prus_qc_noclones --recode vcf --out prus_qc_noclones

plink2 \
  --vcf prus_qc_noclones.vcf \
  --r2-unphased \
  --ld-window-kb 1000 \
  --out test_25_pct_may26_linkage_decay_1000kb \
  --ld-window-r2 0 \
  --threads 64 \
  --thin 0.25
```

The output will look something like this (below is an old result, not the actually-used output):

`head linkage_decay2.vcor`

#CHROM_A	POS_A	ID_A	CHROM_B	POS_B	ID_B	PHASED_R2
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44011	OZ037992.1:44011:C:G	0.209434
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44570	OZ037992.1:44570:C:A	0.220901
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7736	OZ037992.1:7736:T:G	1
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7778	OZ037992.1:7778:T:A	0.340884

Great, now we can write our own python script that will calculate the average r^2 value based on the distances between SNPs. 

## Step 2: Calculate average r^2 values across distances (stepsize = 500)
`cat calculate_decay.py`

```python
import sys
import math

bin_dict={}

for i in range(1,50000,500):
	bin_dict[i, i+499] = [0,0]
	
#print(bin_dict)
#So far we properly set up our distance bins, so that we can later evaluate the avg r^2 for each bin range. 

with open(sys.argv[1]) as file: 
	next(file)
	for line in file:
		l=line.strip().split()	
		posa=int(l[1])
		posb=int(l[4])
		r2=float(l[6])
		dist=abs(posa-posb)

        #for the fencepost issue (disances of exactly base 500) 

        bin_mult = math.ceil(dist/500)
		bin_end = 500*bin_mult
		bin_start = bin_end - 499
		bin_dict[bin_start, bin_end][0]+=r2
		bin_dict[bin_start, bin_end][1]+=1

        #calcuated the distances, and added the current r^2 to the total r^2 and 1 to the counter (+=1) 

for i in bin_dict:
	value=bin_dict[i]
	r2=value[0]
	n=value[1]
	avg=r2/n
	print(i[0],i[1],avg)

#calculated average for each bin :)
```

## Step 3: Plot!
```cat plot.py```
```
 cat plot.py 
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

g = defaultdict(list)
#OZ037994.1	500	0.232311	0.313767
with open(sys.argv[1],'r') as f:
	next(f)
	for line in f:
		l=line.strip().split()
		dist = int(l[1])
		r2 = float(l[2])
		g[dist].append(r2)


x,y=[],[]

for k in g:
	x.append(k)
	y.append(np.average(g[k]))


plt.scatter(x,y,s=1,color='black') 
plt.savefig('prus_ld',dpi=720)
```

Great work! :) Here in our dataset, we see that decay drops to 0.05 at around 100kb. 

<img width="2166" height="1833" alt="decay1000kb_25pct" src="https://github.com/user-attachments/assets/5b33ca2b-0172-40b8-996a-5d71c9626625" />

Note: pruning in the no_maf_no_singletons (all rare alleles except singletons) dataset is tricky, because --indep-pairwise will only take the first 100 snps if you do not provide the 'kb' modifier. Essentially, this will make the window lengths between the pruned and unpruned datasets vastly different and thus uncomparable. To resolve this, ensure you do --indep-pairwise 100kb 1 0.1. This standardizes the window to 100kb, with a stepsize of 1, and an r^2 threshold of 0.1. Use this for comparing between pruned and unpruned datasets in future analyses (like the rare allele pca). Great work!

# Appendix

To make accurate comparisons between the standard-ld-pruned and rare allele-ld-pruned datasets, we specifically use plink2's --indep-pairwise function to set a distance window of 100kb with a r^2 threshold of 0.1 (which are both based on our decay plots). Thus, our command for pruning both standard and rare allele datasets is: 

```bash
--indep-pairwise 100kb 1 0.1.
```

Although it is unlikely to support our data, we also evaluated the standard r^2 threshold of 0.2, which is normally used by the scientific field and literature. Here are our discoveries of SNPs retained between datasets: 

```bash
#		standard-maf_unpruned		standard-maf_ld_r2_0.1		standard-maf_ld_r2_0.2	
		
#snps		13,841,415						533,739						1,021,612

#		prus_rare_unpruned			prus_rare_ld_r2_0.1			prus_rare_ld_r2_0.2	

#snps	 	22,953,951						776,257						19,45,499

```


	
