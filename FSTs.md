# Fst analysis of *Porites rus* 

## Prep the data

1. For FST analysis we'll need to take our non-clonal files and assign their populations to them. 

Here we have our population tsv, followed by our fam file. 

` head Population_IDs.tsv `

```
#IID	SITE
P-rus_022_S25_	SR4
P-rus_023_S26_	SR4
P-rus_077_S28_	SR2
P-rus_078_S29_	SR2
P-rus_079_S30_	SR2
P-rus_080_S31_	SR2
P-rus_081_S32_	SR2
P-rus_082_S33_	SR2
P-rus_083_S34_	SR2

```
#183 samples (includes clones i guess)

` head prus_qc_noclones_ld.fam `

```
0 P-rus_022_S25 0 0 0 -9
0 P-rus_023_S26 0 0 0 -9
0 P-rus_077_S28 0 0 0 -9
0 P-rus_078_S29 0 0 0 -9
0 P-rus_079_S30 0 0 0 -9
0 P-rus_080_S31 0 0 0 -9
0 P-rus_081_S32 0 0 0 -9
0 P-rus_082_S33 0 0 0 -9
0 P-rus_083_S34 0 0 0 -9
```
#159 samples (clones removed) 

Notice there are some discrepancies between IDs, so lets do some qc to ensure:
1) all samples are named the same, and 
2) sames in fam file are in population file (but other way around not necessary)

We can validate those using this:

``` 
awk 'BEGIN{OFS="\t"} NR==1{print; next} {sub(/_$/, "", $1); print}' Population_IDs.tsv > population_IDs.fixed.tsv

cut -f1 population_IDs.fixed.tsv | sort > pop_ids.txt
awk '{print $2}' prus_qc_noclones_ld.fam | sort > fam_ids.txt
comm -3 pop_ids.txt fam_ids.txt
```
For which we find a missing 24 samples, but that is accounted for given 183 in pop file - 159 in the fam file = 24 samples that were removed due to quality or clonality. Great!

Lets run plink to get pairwise fsts:

`plink_fsts.sbatch`

```
#!/bin/bash
#SBATCH --job-name=fsts
#SBATCH -p serc,spalumbi,hns
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --mem=96G
#SBATCH --output=fst.out
#SBATCH --error=fst.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -eo pipefail

module load plink/2.0a7

plink2 \
  --bfile prus_qc_noclones_ld \
  --fst SITE \
  --pheno population_IDs.fixed.tsv \
  --pheno-name SITE \
  --out prus_site_fst
```

And let's check our result: 

`head sorted_site_fsts.summary `

```
BIKIRIN	CMI	0.00814853
CMI	RANDOM	0.00805306
CMI	SITE	0.00599581
AJELTAKE	CMI	0.00589973
CMI	SHERWOODS	0.00401378
AIRPORT	CMI	0.00382265
CMI	SR3	0.00120066
CMI	TOBLAR	0.0011752
CMI	SR2	0.000931007
CMI	SR6	-0.00157603
```
