Summary stats for our snpEff output: top_5_pct_output.snpEff.vcf

cd /scratch/users/kaiku/snpEff/prus_analysis/top5_pct_gemma

INPUT_VCF="top_5_pct_output.snpEff.vcf"
PREFIX="top5_pct_exp_gemma"

# ------------------------------------------------------------
# Extract VCF records that contain at least one HIGH annotation
# ------------------------------------------------------------

grep "^#" ${INPUT_VCF} > ${PREFIX}.HIGH.vcf

grep -v "^#" ${INPUT_VCF} \
  | awk '$8 ~ /HIGH/' \
  >> ${PREFIX}.HIGH.vcf


# ---------------------------------------------------------------------
# Extract VCF records that contain at least one HIGH or MODERATE annotation
# ---------------------------------------------------------------------

grep "^#" ${INPUT_VCF} > ${PREFIX}.HIGH_MODERATE.vcf

grep -v "^#" ${INPUT_VCF} \
  | awk '$8 ~ /HIGH|MODERATE/' \
  >> ${PREFIX}.HIGH_MODERATE.vcf


# ---------------------------------------------------------------------
# Create HIGH annotation-only table
# This prints only ANN records where IMPACT == HIGH
# ---------------------------------------------------------------------

grep -v "^#" ${INPUT_VCF} \
  | awk 'BEGIN {
      FS=OFS="\t"
      print "CHROM","POS","ID","REF","ALT","QUAL","EFFECT","IMPACT",
            "GENE_NAME","GENE_ID","TRANSCRIPT_ID","HGVS_C","HGVS_P"
    }
    {
      ann=$8
      sub(/^.*ANN=/,"",ann)
      sub(/;.*/,"",ann)
      n=split(ann,a,",")

      for (i=1; i<=n; i++) {
        split(a[i],f,"|")

        if (f[3] == "HIGH") {
          print $1,$2,$3,$4,$5,$6,
                f[2],f[3],f[4],f[5],f[7],f[10],f[11]
        }
      }
    }' \
  > ${PREFIX}.HIGH.annotations_only.tsv


# ---------------------------------------------------------------------
# Create HIGH + MODERATE annotation-only table
# This prints only ANN records where IMPACT == HIGH or MODERATE
# ---------------------------------------------------------------------

grep -v "^#" ${INPUT_VCF} \
  | awk 'BEGIN {
      FS=OFS="\t"
      print "CHROM","POS","ID","REF","ALT","QUAL","EFFECT","IMPACT",
            "GENE_NAME","GENE_ID","TRANSCRIPT_ID","HGVS_C","HGVS_P"
    }
    {
      ann=$8
      sub(/^.*ANN=/,"",ann)
      sub(/;.*/,"",ann)
      n=split(ann,a,",")

      for (i=1; i<=n; i++) {
        split(a[i],f,"|")

        if (f[3] == "HIGH" || f[3] == "MODERATE") {
          print $1,$2,$3,$4,$5,$6,
                f[2],f[3],f[4],f[5],f[7],f[10],f[11]
        }
      }
    }' \
  > ${PREFIX}.HIGH_MODERATE.annotations_only.tsv


# ---------------------------------------------------------------------
# Create expanded table containing every annotation for HIGH/MODERATE sites
# This is useful for context, but not for candidate gene filtering
# ---------------------------------------------------------------------

grep -v "^#" ${PREFIX}.HIGH_MODERATE.vcf \
  | awk 'BEGIN {
      FS=OFS="\t"
      print "CHROM","POS","ID","REF","ALT","QUAL","EFFECT","IMPACT",
            "GENE_NAME","GENE_ID","TRANSCRIPT_ID","HGVS_C","HGVS_P"
    }
    {
      ann=$8
      sub(/^.*ANN=/,"",ann)
      sub(/;.*/,"",ann)
      n=split(ann,a,",")

      for (i=1; i<=n; i++) {
        split(a[i],f,"|")

        print $1,$2,$3,$4,$5,$6,
              f[2],f[3],f[4],f[5],f[7],f[10],f[11]
      }
    }' \
  > ${PREFIX}.HIGH_MODERATE.all_annotations_for_sites.tsv


# ---------------------------------------------------------------------
# Create unique gene-ID and transcript-ID lists from corrected tables
# ---------------------------------------------------------------------

cut -f10 ${PREFIX}.HIGH.annotations_only.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > ${PREFIX}.HIGH.gene_ids.txt

cut -f11 ${PREFIX}.HIGH.annotations_only.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > ${PREFIX}.HIGH.transcript_ids.txt

cut -f10 ${PREFIX}.HIGH_MODERATE.annotations_only.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > ${PREFIX}.HIGH_MODERATE.gene_ids.txt

cut -f11 ${PREFIX}.HIGH_MODERATE.annotations_only.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > ${PREFIX}.HIGH_MODERATE.transcript_ids.txt


# ---------------------------------------------------------------------
# Write readable summary report
# ---------------------------------------------------------------------

{
  echo "SnpEff GEMMA TOP 5% EXPERIENCED RESULTS"
  echo "======================================="
  echo

  echo "INPUT VCF"
  echo "---------"
  echo "${INPUT_VCF}"
  echo

  echo "TOTAL VARIANTS"
  echo "--------------"
  grep -vc "^#" ${INPUT_VCF}
  echo

  echo "VARIANTS WITH ANN ANNOTATIONS"
  echo "-----------------------------"
  grep -v "^#" ${INPUT_VCF} | grep -c "ANN="
  echo

  echo "PRIMARY EFFECT PER SNP"
  echo "----------------------"
  grep -v "^#" ${INPUT_VCF} \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f2 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "PRIMARY IMPACT PER SNP"
  echo "----------------------"
  grep -v "^#" ${INPUT_VCF} \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f3 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "ALL ANNOTATION EFFECTS"
  echo "----------------------"
  grep -v "^#" ${INPUT_VCF} \
    | grep -o "ANN=[^;]*" \
    | sed 's/^ANN=//' \
    | tr ',' '\n' \
    | cut -d"|" -f2 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "ALL ANNOTATION IMPACTS"
  echo "----------------------"
  grep -v "^#" ${INPUT_VCF} \
    | grep -o "ANN=[^;]*" \
    | sed 's/^ANN=//' \
    | tr ',' '\n' \
    | cut -d"|" -f3 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "NUMBER OF VARIANT SITES WITH AT LEAST ONE HIGH ANNOTATION"
  echo "---------------------------------------------------------"
  grep -vc "^#" ${PREFIX}.HIGH.vcf
  echo

  echo "NUMBER OF VARIANT SITES WITH AT LEAST ONE HIGH OR MODERATE ANNOTATION"
  echo "---------------------------------------------------------------------"
  grep -vc "^#" ${PREFIX}.HIGH_MODERATE.vcf
  echo

  echo "NUMBER OF HIGH ANNOTATION RECORDS"
  echo "---------------------------------"
  tail -n +2 ${PREFIX}.HIGH.annotations_only.tsv | wc -l
  echo

  echo "NUMBER OF HIGH + MODERATE ANNOTATION RECORDS"
  echo "--------------------------------------------"
  tail -n +2 ${PREFIX}.HIGH_MODERATE.annotations_only.tsv | wc -l
  echo

  echo "HIGH/MODERATE IMPACT COUNTS FROM CORRECTED ANNOTATION TABLE"
  echo "-----------------------------------------------------------"
  cut -f8 ${PREFIX}.HIGH_MODERATE.annotations_only.tsv \
    | sort \
    | uniq -c
  echo

  echo "HIGH/MODERATE EFFECT COUNTS FROM CORRECTED ANNOTATION TABLE"
  echo "-----------------------------------------------------------"
  cut -f7 ${PREFIX}.HIGH_MODERATE.annotations_only.tsv \
    | tail -n +2 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "NUMBER OF UNIQUE HIGH-IMPACT GENES"
  echo "----------------------------------"
  wc -l < ${PREFIX}.HIGH.gene_ids.txt
  echo

  echo "NUMBER OF UNIQUE HIGH-IMPACT TRANSCRIPTS"
  echo "----------------------------------------"
  wc -l < ${PREFIX}.HIGH.transcript_ids.txt
  echo

  echo "NUMBER OF UNIQUE HIGH + MODERATE GENES"
  echo "--------------------------------------"
  wc -l < ${PREFIX}.HIGH_MODERATE.gene_ids.txt
  echo

  echo "NUMBER OF UNIQUE HIGH + MODERATE TRANSCRIPTS"
  echo "--------------------------------------------"
  wc -l < ${PREFIX}.HIGH_MODERATE.transcript_ids.txt
  echo

  echo "HIGH-IMPACT ANNOTATION TABLE"
  echo "----------------------------"
  column -t ${PREFIX}.HIGH.annotations_only.tsv
  echo

  echo "OUTPUT FILES CREATED"
  echo "--------------------"
  echo "${PREFIX}.HIGH.vcf"
  echo "${PREFIX}.HIGH_MODERATE.vcf"
  echo "${PREFIX}.HIGH.annotations_only.tsv"
  echo "${PREFIX}.HIGH_MODERATE.annotations_only.tsv"
  echo "${PREFIX}.HIGH_MODERATE.all_annotations_for_sites.tsv"
  echo "${PREFIX}.HIGH.gene_ids.txt"
  echo "${PREFIX}.HIGH.transcript_ids.txt"
  echo "${PREFIX}.HIGH_MODERATE.gene_ids.txt"
  echo "${PREFIX}.HIGH_MODERATE.transcript_ids.txt"

} > ${PREFIX}.snpeff.summary.txt

echo "Finished. Read the report with:"
echo "cat ${PREFIX}.snpeff.summary.txt"
