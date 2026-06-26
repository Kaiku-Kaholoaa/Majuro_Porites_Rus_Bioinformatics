# Using snpEff For Variant Subset Annotation

Earlier we conducted a GWAS via gemma to identify peaks associated with the top5% heat experienced by corals. We also used pcadapt to
evaluate outlier snps associated with our first two principle components from our PCA, which showed a slight grouping with a thermal outgroup. 

We also learned that minor allele frequencies can identify hidden trends in human population structure, and can be extremely useful for
our analysis given our fine spatial scale (< 40km) and our highly-related, mostly panmictic population. The minor allele frequency pca 
identified that there is a clear subgroup among our thermal individuals, to which we grouped together based on a PC1 threshold. 
The alleles of this subgroup (n_indiv = 19; total_possible_alleles = 38), were then filtered for private alleles at high frequency. This
represents alelles that are not present in the greater population, display a frequency of > 20%, and have a total mimimum allele count of 8.
In the end we idenitfied ~800 alleles at high frequency not present in the main population. 

The snps from each analysis (GEMMA-GWAS, PCAdapt, and Private Minor Allele Analysis) were then converted into two different vcf files by filtering the comprehensive population vcf by the positions of SNPs of interest. 

snpEff can then annotate these vcf files by mapping snp positions to an annotated database which can idenitfy snp locations (intergenic, intron, upstream/downstream variants, etc.), and mutational changes of interest (including synonymous, missense, stop_gained, splice_acceptor/donor variants, etc.). 

However, the database must first be created using the gff that we made of the Porites rus genome in our Gene_Annotation_Pipeline.md.

## Create the Database

`cd /scratch/users/kaiku/snpEff `
`mkdir data/prus`

In the prus directory, add both the reference.fna and our GFF file. 

Note: GFF file must not include ##fasta lines at the end, and must only contain the final main features.
* Final main features include: exon, CDS, mRNA, gene, three_prime_UTR, five_prime_UTR, etc. 
* Non-gene model features are the lines of evidence from maker2, including: match_part, protein_match, expressed_sequence_match, match, translated_nucleotide_match, etc.

* From our maker2.gff file, we had:
  
```
328774 exon
324579 CDS
92944 mRNA
92944 gene
13429 three_prime_UTR
13388 five_prime_UTR
```

Note: Esure that the IDs are consistent between the .vcf, .gff and .fna files

```bash
# VCF contigs
bcftools view -H private_alleles_for_snpEff.vcf.gz \
  | cut -f1 \
  | sort -u \
  > vcf.contigs.txt

# GFF contigs
grep -v "^#" /scratch/users/kaiku/snpEff/data/prus/genes.gff \
  | cut -f1 \
  | sort -u \
  > gff.contigs.txt

# FASTA contigs
grep "^>" /scratch/users/kaiku/snpEff/data/prus/sequences.fa \
  | sed 's/^>//' \
  | cut -d" " -f1 \
  | sort -u \
  > fasta.contigs.txt
```

Note: Finally, snpEff expects specifc names for the reference and gff file, which we can symlink using

```bash
ln -s genes.clean.gff genes.gff
ln -s porites_rus_reference_genome.fna sequences.fa 
```

And now we can build our database of annotations! Do this using the .jar file with the build command. 

Here the -Xmx16g flag requests 16G of memory, with the -gff3 flag since that is what we're building our database on, and -v for verbose messaging. Additionally, snpEff normally checks created databases against a known CDS or Protein file, but since we do not have that, we just omit and continue. Finally we call it prus for consistency. 

```bash
cd /scratch/users/kaiku/snpEff

java -Xmx16g -jar snpEff.jar build \
  -gff3 \
  -v \
  -noCheckCds \
  -noCheckProtein \
  prus
```

And let's check the output exists: 

`ls -lh data/prus/snpEffectPredictor.bin`

Great! Now we can start annotating:

## Annotations and Summary Statistics

To annotate using snpEff, we once again call the .jar file but this time use the ann command. Here we'll just focus on the private snp analysis for simplicity. Our vcf contains 799 snps. 

```bash
cd /scratch/users/kaiku/snpEff/prus_analysis

java -Xmx8g -jar ../snpEff.jar ann \
  -noStats \ #Tells SnpEff not to generate the HTML summary/statistics report. This avoids extra files and issues.
  -v \ #verbose
  prus \ #database
  private_alleles_for_snpEff.vcf.gz \ #input vcf
  > private_snps.snpeff.vcf #output vcf
```

And we can count the number of successful annotations:

```bash
echo "Total variants:"
grep -vc "^#" private_snps.snpeff.vcf
#counts all lines that are not do not start (^) with '#'.
```

All 799 snps were successfully annotated! Nice. 

This script will help with our annotation summary: 

```bash
cd /scratch/users/kaiku/snpEff/prus_analysis

# Extract HIGH-impact variants
grep "^#" private_snps.snpeff.vcf > private_snps.HIGH.vcf

grep -v "^#" private_snps.snpeff.vcf \
  | awk '$8 ~ /HIGH/' \
  >> private_snps.HIGH.vcf

# Extract HIGH + MODERATE variants
grep "^#" private_snps.snpeff.vcf > private_snps.HIGH_MODERATE.vcf

grep -v "^#" private_snps.snpeff.vcf \
  | awk '$8 ~ /HIGH|MODERATE/' \
  >> private_snps.HIGH_MODERATE.vcf

# Create HIGH main-annotation table
grep -v "^#" private_snps.HIGH.vcf \
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
  > private_snps.HIGH.main_annotation.tsv

# Create HIGH + MODERATE main-annotation table
grep -v "^#" private_snps.HIGH_MODERATE.vcf \
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
  > private_snps.HIGH_MODERATE.main_annotation.tsv

# Create expanded table containing every annotation
grep -v "^#" private_snps.HIGH_MODERATE.vcf \
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
  > private_snps.HIGH_MODERATE.all_annotations.tsv

# Create unique gene-ID lists
cut -f10 private_snps.HIGH.main_annotation.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > private_snps.HIGH.gene_ids.txt

cut -f10 private_snps.HIGH_MODERATE.main_annotation.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > private_snps.HIGH_MODERATE.gene_ids.txt

# Write all results into one readable text file
{
  echo "SnpEff PRIVATE-SNP RESULTS"
  echo "=========================="
  echo

  echo "TOTAL VARIANTS"
  echo "--------------"
  grep -vc "^#" private_snps.snpeff.vcf
  echo

  echo "VARIANTS WITH ANN ANNOTATIONS"
  echo "-----------------------------"
  grep -v "^#" private_snps.snpeff.vcf | grep -c "ANN="
  echo

  echo "PRIMARY EFFECT PER SNP"
  echo "----------------------"
  grep -v "^#" private_snps.snpeff.vcf \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f2 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "PRIMARY IMPACT PER SNP"
  echo "----------------------"
  grep -v "^#" private_snps.snpeff.vcf \
    | grep -o "ANN=[^;]*" \
    | cut -d"|" -f3 \
    | sort \
    | uniq -c \
    | sort -nr
  echo

  echo "ALL ANNOTATION EFFECTS"
  echo "----------------------"
  grep -v "^#" private_snps.snpeff.vcf \
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
  grep -v "^#" private_snps.snpeff.vcf \
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
  grep -vc "^#" private_snps.HIGH.vcf
  echo

  echo "NUMBER OF HIGH + MODERATE VARIANTS"
  echo "----------------------------------"
  grep -vc "^#" private_snps.HIGH_MODERATE.vcf
  echo

  echo "NUMBER OF UNIQUE HIGH-IMPACT GENES"
  echo "----------------------------------"
  wc -l < private_snps.HIGH.gene_ids.txt
  echo

  echo "NUMBER OF UNIQUE HIGH + MODERATE GENES"
  echo "--------------------------------------"
  wc -l < private_snps.HIGH_MODERATE.gene_ids.txt
  echo

  echo "HIGH-IMPACT VARIANT TABLE"
  echo "-------------------------"
  column -t private_snps.HIGH.main_annotation.tsv
  echo

  echo "OUTPUT FILES CREATED"
  echo "--------------------"
  echo "private_snps.HIGH.vcf"
  echo "private_snps.HIGH_MODERATE.vcf"
  echo "private_snps.HIGH.main_annotation.tsv"
  echo "private_snps.HIGH_MODERATE.main_annotation.tsv"
  echo "private_snps.HIGH_MODERATE.all_annotations.tsv"
  echo "private_snps.HIGH.gene_ids.txt"
  echo "private_snps.HIGH_MODERATE.gene_ids.txt"

} > private_snps.snpeff.summary.txt

echo "Finished. Read the report with:"
echo "cat private_snps.snpeff.summary.txt"
```

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
