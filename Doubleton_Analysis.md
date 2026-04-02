# Doubleton Analysis of Majuro P.rus

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
  --r2-phased \
  --ld-window-kb 50 \
  --out linkage_decay2
```

The output will look something like this:

`head linkage_decay2.vcor`

#CHROM_A	POS_A	ID_A	CHROM_B	POS_B	ID_B	PHASED_R2
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44011	OZ037992.1:44011:C:G	0.209434
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44570	OZ037992.1:44570:C:A	0.220901
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7736	OZ037992.1:7736:T:G	1
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7778	OZ037992.1:7778:T:A	0.340884
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7827	OZ037992.1:7827:T:C	0.305002
OZ037992.1	7736	OZ037992.1:7736:T:G	OZ037992.1	7778	OZ037992.1:7778:T:A	0.340959
OZ037992.1	7736	OZ037992.1:7736:T:G	OZ037992.1	7827	OZ037992.1:7827:T:C	0.305099
OZ037992.1	7739	OZ037992.1:7739:G:A	OZ037992.1	7749	OZ037992.1:7749:A:G	0.748072
OZ037992.1	7739	OZ037992.1:7739:G:A	OZ037992.1	7806	OZ037992.1:7806:T:G	0.555487

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
Great work! :) All that's left to do is plot!

