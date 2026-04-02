# Documentation for QA/QC

## Pipeline for quality filters:

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

Great work! :)

