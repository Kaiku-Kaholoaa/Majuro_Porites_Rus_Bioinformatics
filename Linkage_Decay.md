# Documentation for Linkage Decay and Doubleton Analysis

Linkage Decay

## Proper Filtering and Clone Removal (Business as usual) 

`cat pipeline.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=prus_final_qc_pca
#SBATCH -p serc,hns,spalumbi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

# ===============================
# INPUT
# ===============================
VCF_IN="bcf_total_compressed_prus_dec25.vcf.gz"

ml biology
ml bcftools
ml vcftools
ml plink/2.0a7

echo "========================================"
echo " STEP 1: BCFTOOLS VARIANT FILTERING"
echo "========================================"

# --- 1A: Biallelic SNPs ---
bcftools view -m2 -M2 -v snps \
    ${VCF_IN} \
    -Oz -o tmp_biallelic.vcf.gz
bcftools index tmp_biallelic.vcf.gz

# --- 1B: QUAL + DP filters ---
bcftools filter -i 'QUAL>30 && DP>5' \
    tmp_biallelic.vcf.gz \
    -Oz -o tmp_filtered.vcf.gz
bcftools index tmp_filtered.vcf.gz

# --- 1C: Site-level missingness ---
vcftools \
    --gzvcf tmp_filtered.vcf.gz \
    --max-missing 0.80 \
    --recode --recode-INFO-all \
    --out tmp_mm80

bgzip -f tmp_mm80.recode.vcf
tabix -p vcf tmp_mm80.recode.vcf.gz

FILTERED_VCF="tmp_mm80.recode.vcf.gz"

echo "Done filtering. Output = $FILTERED_VCF"
echo ""

echo "========================================"
echo " STEP 2: CONVERT TO PGEN"
echo "========================================"

plink2 \
    --vcf ${FILTERED_VCF} \
    --make-pgen \
    --out tmp_cleaned

echo ""

echo "========================================"
echo " STEP 3: SET UNIQUE VARIANT IDS"
echo "========================================"

plink2 \
    --pfile tmp_cleaned \
    --set-all-var-ids @:#:\$r:\$a \
    --make-pgen \
    --out prus_cleaned_unique

echo ""

echo "========================================"
echo " STEP 4: PLINK2 SAMPLE + SNP FILTERS"
echo "========================================"
# These match best-practice thresholds

plink2 \
    --pfile prus_cleaned_unique \
    --geno 0.05 \
    --mind 0.20 \
    --maf 0.01 \
    --max-alleles 2 \
    --make-pgen \
    --out prus_qc_final

echo ""

echo "========================================"
echo " STEP 5: KING RELATEDNESS (cutoff = 0.40)"
echo "========================================"

# Identify clones (retain <= 0.4)
plink2 \
    --pfile prus_qc_final \
    --king-cutoff 0.40 \
    --out prus_qc_final_king

# Create clone/duplicate-removed dataset
plink2 \
    --pfile prus_qc_final \
    --keep prus_qc_final_king.king.cutoff.in.id \
    --make-pgen \
    --out prus_qc_noclones

echo ""

echo "========================================"
echo " CLEANING TEMP FILES"
echo "========================================"

rm -f tmp_biallelic.vcf.gz* \
      tmp_filtered.vcf.gz* \
      tmp_mm80.recode.vcf.gz* \
      tmp_cleaned.*

echo "Temp files removed."
echo ""

echo "🎉 DONE! Final outputs:"
echo "  prus_qc_final.*           = fully QC'd dataset"
echo "  prus_qc_noclones.*        = clone/relative filtered dataset"
echo ""
```

Great, now with our fully qa/qc'ed dataset, lets use plink to identify correlated snps. 

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
`#CHROM_A	POS_A	ID_A	CHROM_B	POS_B	ID_B	PHASED_R2
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44011	OZ037992.1:44011:C:G	0.209434
OZ037992.1	7731	OZ037992.1:7731:T:C	OZ037992.1	44570	OZ037992.1:44570:C:A	0.220901
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7736	OZ037992.1:7736:T:G	1
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7778	OZ037992.1:7778:T:A	0.340884
OZ037992.1	7734	OZ037992.1:7734:T:C	OZ037992.1	7827	OZ037992.1:7827:T:C	0.305002
OZ037992.1	7736	OZ037992.1:7736:T:G	OZ037992.1	7778	OZ037992.1:7778:T:A	0.340959
OZ037992.1	7736	OZ037992.1:7736:T:G	OZ037992.1	7827	OZ037992.1:7827:T:C	0.305099
OZ037992.1	7739	OZ037992.1:7739:G:A	OZ037992.1	7749	OZ037992.1:7749:A:G	0.748072
OZ037992.1	7739	OZ037992.1:7739:G:A	OZ037992.1	7806	OZ037992.1:7806:T:G	0.555487`

Great, now we can write our own python script that will calculate the average r^2 value based on the distances between SNPs. 
Conceptually, snps in the first bin (distance of 0-500) will have a very high avg correlation (r^2), and as you increase 
your distance your avg correlations will start to decrease, and eventually plateau. The point of the linkage decay analysis 
is to identify the minimum distance upon which your correlations stabilize, so that you can optimize your pruning. Knowing the
optimal distance to prune is essential for the doubleton analysis, since we'll be including our minor alleles (excluded in this 
analysis), but need to properly filter for linkage. 

`cat calculate_decay.py`

```bash
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

Great work! :) All that's left to do is plot!


