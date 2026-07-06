# Evaluating SNPs of Interest from Unpruned, No MAF-Filtering
Using the unpruned no MAF filtering dataset, we created three different datasets for evaluation, upon which we find:

  1. (Priority) ~1,000 alleles at low frequency (less than 10%) in the main population, that are at high frequency (greater than 80%) in the thermally-associated subpopulation, 
  2. (Also very interesting) ~4,000 private alleles (not present in the main population) that are at a frequency of greater than 20%,
  3. (Non-priority, but worth exploring) ~400 alleles that are fixed in the subpopulation, and have a 20% increase in AF relative to the main population.

PCA: 
<img width="3000" height="2100" alt="mac2_no_maf_unpruned" src="https://github.com/user-attachments/assets/0e0859fc-4b25-498e-b8f3-e512dae2e910" />

## Annotations and Summary Statistics

### Alleles at high frequency in subpopulation, but low frequency in major population

For this analysis, our vcf contains 1,039 alleles

```bash
cd /scratch/users/kaiku/snpEff/prus_analysis

ml java/21.0.4  

java -Xmx8g -jar ../snpEff.jar ann \
  -noStats \ #Tells SnpEff not to generate the HTML summary/statistics report. This avoids extra files and issues.
  -v \ #verbose
  prus \ #database
  major_allele_differences_no_singletons_no_maf_unpruned.vcf.gz \ #input vcf
  > major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf #output vcf
```

And we can count the number of successful annotations:

```bash
echo "Total variants:"
grep -vc "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf 
#counts all lines that are not do not start (^) with '#'.
```

All 1039 snps were successfully annotated! Nice. 

As a reminder, the output is a vcf so that means that column 8 contains the annotations. 

To view only the annotations, we can isolate the variant lines, then remove everything before the annotations via 
`sed s/.*ANN=//`
and anything after the annotations by removing anything following the semicolon:
`sed s/;.*//`

Great, now these are all the annotations. But most snps will have multiple annotations, so we can cut just the first field which is the most important, remembering that snpEff separates annotations via comma:

`cut -f1 -d ","`

This should leave each snpEff field separated by "|". we can replace this with a new line for ease:

`tr "|" "\n" `

For which we can number our new lines using:

`nl -ba`

Awesome, this is complete! As a reminder, snpEff usually follows the following format:

```bash
1   Allele #alt allele being annotated
2   Effect #biological consequence, like missense or stop_gain_variant
3   Impact #predicted impact based the biological consequence (high impact, moderate, low, etc.)
4   Gene name #name of the associated gene from the database
5   Gene ID #unique gene ID that will be used for the functional annotation
6   Feature type #transcript, gene, cds, intron, etc. Since maker2 predicted based on transcripts, most would be transcripts. 
7   Transcript ID #id for connecting transcipt id to blast results
8   Transcript biotype #class of transcript. ex: trna, mirna, rrna or protein_coding. 
9   Exon/intron rank #exon1 = first exon; intron 1 = first intron; intron5 = 5th intron and is between exon 4 and 5. 
10  HGVS.c #human genome variation soc nomenclature, coding region pos. ex: c.1964C>G = at pos 1964 C turns to G
11  HGVS.p #HGVS nomenclature, protein pos. ex: p.Ser1002* = AA change of Serine to stop codon at AA 1002. 
12  cDNA position/length #position along cDNA transcript. ex: 1803 (pos) /2877 (total cDNA length)
13  CDS position/length #position along CDS. ex: 1694 (pos) / 2631 (total cds length) 
14  Protein position/length #position along AA string (protein). ex: 565 (AA) / 876 (total AA length)
15  Distance #from gene feature (gene, cds, etc.), usually used for up/downstream variants. ex: 1808 = bp up/downstream. 
16  Warnings/errors # transcript specific warning. ex: WARNING_NO_START_CODON, WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS
```
Usually, we care about fields 2, 3, 5, 7, 10, 11, 14, 16

Cool! 

## This script will help with our annotation summary: 

```bash

# Extract HIGH-impact variants
grep "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  > major_allele_differences_NSNMU.snpeff.HIGH.vcf

grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  \
  | awk '$8 ~ /HIGH/' \
  >> major_allele_differences_NSNMU.HIGH.vcf

# Extract HIGH + MODERATE variants
grep "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf > major_allele_differences_NSNMU.HIGH_MODERATE.vcf

grep -v "^#" private_snps.snpeff.vcf \
  | awk '$8 ~ /HIGH|MODERATE/' \
  >> major_allele_differences_NSNMU.HIGH_MODERATE.vcf

# Create HIGH main-annotation table
grep -v "^#" major_allele_differences_NSNMU.HIGH.vcf \
  | awk 'BEGIN {
      FS=OFS="\t"
      print "CHROM","POS","ID","REF","ALT","QUAL","EFFECT","IMPACT",
            "GENE_NAME","GENE_ID","TRANSCRIPT_ID","HGVS_C","HGVS_P"
    }
    {
      ann=$8
      sub(/^.*ANN=/,"",ann)
      sub(/;.*/,"",ann)
      split(ann,a,",")
      split(a[1],f,"|")

      print $1,$2,$3,$4,$5,$6,
            f[2],f[3],f[4],f[5],f[7],f[10],f[11]
    }' \
  > major_allele_differences_NSNMU.HIGH.main_annotation.tsv

# Create HIGH + MODERATE main-annotation table
grep -v "^#" major_allele_differences_NSNMU.HIGH_MODERATE.vcf \
  | awk 'BEGIN {
      FS=OFS="\t"
      print "CHROM","POS","ID","REF","ALT","QUAL","EFFECT","IMPACT",
            "GENE_NAME","GENE_ID","TRANSCRIPT_ID","HGVS_C","HGVS_P"
    }
    {
      ann=$8
      sub(/^.*ANN=/,"",ann)
      sub(/;.*/,"",ann)
      split(ann,a,",")
      split(a[1],f,"|")

      print $1,$2,$3,$4,$5,$6,
            f[2],f[3],f[4],f[5],f[7],f[10],f[11]
    }' \
  > major_allele_differences_NSNMU.HIGH_MODERATE.main_annotation.tsv

# Create expanded table containing every annotation
grep -v "^#" major_allele_differences_NSNMU.HIGH_MODERATE.vcf \
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
  > major_allele_differences_NSNMU.HIGH_MODERATE.all_annotations.tsv

# Create unique gene-ID lists
cut -f10 major_allele_differences_NSNMU.HIGH.main_annotation.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > major_allele_differences_NSNMU.HIGH.gene_ids.txt

cut -f10 major_allele_differences_NSNMU.HIGH_MODERATE.all_annotations.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > major_allele_differences_NSNMU.HIGH_MODERATE.gene_ids.txt

# Write all results into one readable text file
{
  echo "SnpEff PRIVATE-SNP RESULTS"
  echo "=========================="
  echo

  echo "TOTAL VARIANTS"
  echo "--------------"
  grep -vc "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf 
  echo

  echo "VARIANTS WITH ANN ANNOTATIONS"
  echo "-----------------------------"
  grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  | grep -c "ANN="
  echo

  echo "PRIMARY EFFECT PER SNP"
  echo "----------------------"
  grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f2 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "PRIMARY IMPACT PER SNP"
  echo "----------------------"
  grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f3 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "ALL ANNOTATION EFFECTS"
  echo "----------------------"
  grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  \
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
  grep -v "^#" major_allele_differences_no_singletons_no_maf_unpruned.snpeff.vcf  \
    | grep -o "ANN=[^;]*" \
    | sed 's/^ANN=//' \
    | tr ',' '\n' \
    | cut -d"|" -f3 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "NUMBER OF HIGH-IMPACT VARIANTS"
  echo "------------------------------"
  grep -vc "^#" major_allele_differences_NSNMU.HIGH.vcf
  echo

  echo "NUMBER OF HIGH + MODERATE VARIANTS"
  echo "----------------------------------"
  grep -vc "^#" major_allele_differences_NSNMU.HIGH_MODERATE.vcf
  echo

  echo "NUMBER OF UNIQUE HIGH-IMPACT GENES"
  echo "----------------------------------"
  wc -l < major_allele_differences_NSNMU.HIGH.gene_ids.txt
  echo

  echo "NUMBER OF UNIQUE HIGH + MODERATE GENES"
  echo "--------------------------------------"
  wc -l < pmajor_allele_differences_NSNMU.HIGH_MODERATE.gene_ids.txt
  echo

  echo "HIGH-IMPACT VARIANT TABLE"
  echo "-------------------------"
  column -t major_allele_differences_NSNMU.HIGH.main_annotation.tsv
  echo

  echo "OUTPUT FILES CREATED"

} > major_allele_differences_NSNMU.snpeff.summary.txt

echo "Finished. Read the report with:"
echo "cat major_allele_differences_NSNMU.snpeff.summary.txt"
```
## Summary Results 

```bash
cat private_snps.snpeff.summary.txt

SnpEff PRIVATE-SNP RESULTS
==========================

Total variants:
799

Variants with ANN annotations:
799

PRIMARY EFFECT PER SNP:
    300 upstream_gene_variant
    219 intron_variant
    135 downstream_gene_variant
     62 missense_variant
     30 intergenic_region
     29 synonymous_variant
      9 3_prime_UTR_variant
      4 stop_gained
      4 splice_region_variant&intron_variant
      2 splice_donor_variant&intron_variant
      2 missense_variant&splice_region_variant
      2 5_prime_UTR_variant
      1 5_prime_UTR_premature_start_codon_gain_variant

PRIMARY IMPACT PER SNP:
    695 MODIFIER
     64 MODERATE
     34 LOW
      6 HIGH

ALL ANNOTATION EFFECTS:
    466 downstream_gene_variant
    455 upstream_gene_variant
    428 intron_variant
    283 intergenic_region
     62 missense_variant
     30 synonymous_variant
      9 3_prime_UTR_variant
      4 stop_gained
      4 splice_region_variant&intron_variant
      3 5_prime_UTR_variant
      2 splice_donor_variant&intron_variant
      2 missense_variant&splice_region_variant
      1 5_prime_UTR_premature_start_codon_gain_variant

ALL ANNOTATION IMPACTS:
   1644 MODIFIER
     64 MODERATE
     35 LOW
      6 HIGH

HIGH-impact variants:
6

HIGH + MODERATE variants:
70

Unique HIGH-impact genes:
6

Unique HIGH + MODERATE genes:
69
```
