To evaluate gene functions, we will first gather the transcript_ids from snpEff.

```bash
cut -f11 major_allele_differences_NSNMU.HIGH_MODERATE.all_annotations.tsv \
  | tail -n +2 \
  | grep -v "^$" \
  | sort -u \
  > major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt

wc -l major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt
head major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt
```

We'll also just double check that the ids match the interpro.tsv file, using both the gene and transcript ids because sometimes one produces more than the other:

```wc -l *MODERATE*_ids.txt 
 # 26 major_allele_differences_NSNMU.HIGH_MODERATE.gene_ids.txt
 # 26 major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt
#Seems great, both show our 26 genes.

grep -F -f major_allele_differences_NSNMU.HIGH_MODERATE.gene_ids.txt ../../Past.interpro.all.tsv | wc -l #147
grep -F -f major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt ../../Past.interpro.all.tsv | wc -l #147
#and both id files produce the same number of interpro IDS which is awesome. This is likely due to multiple databases having a result for each gene. It #could also be because multiple genes could be producing different transcripts, although we dont really see that here.
```

Awesome! Looks like the IDs match the interpro.tsv, and we can now pull interpro ids using the snpEff transcript ids.

```bash
grep -F -f major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt Past.interpro.all.tsv \
  > major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv
```

and we can check our results
```bash
wc -l major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv #147, nice!
cut -f1 major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv | sort -u | wc -l
#19 unique genes. Interesting, 7 genes are not well known apparently. 
cut -f1 major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv | sort -u | grep -c "38003"
#10 unique genes. Wow, 10 of our genes are within this interesting chr 12 region. 

```
Looks like of our 26 genes, only 19 were identified via interproscan. 

Then we'll pull interproscan ids using the transcript ids from snpEFF
```
PROTEINS="/scratch/users/kaiku/may18_26_pipeline/transcriptomics/dec25_maker2_mpi/functional_annotation/no_mpi_round3.2.all.maker.proteins.fasta"

awk -v IDS="major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt" '
BEGIN {
  while ((getline id < IDS) > 0) keep[id]=1
}
^>/ {
  id=$1
  sub(/^>/, "", id)
  print_record = (id in keep)
}
print_record {
  print
}
' "$PROTEINS" > major_allele_differences_NSNMU.HIGH_MODERATE.proteins.faa

grep -c "^>" major_allele_differences_NSNMU.HIGH_MODERATE.proteins.faa
```

check missing IDs

```bash
grep "^>" major_allele_differences_NSNMU.HIGH_MODERATE.proteins.faa \
  | sed 's/^>//; s/[[:space:]].*//' \
  | sort -u \
  > recovered_HIGH_MODERATE.protein_ids.txt

comm -23 \
  <(sort major_allele_differences_NSNMU.HIGH_MODERATE.transcript_ids.txt) \
  <(sort recovered_HIGH_MODERATE.protein_ids.txt) \
  > missing_HIGH_MODERATE.protein_ids.txt

wc -l missing_HIGH_MODERATE.protein_ids.txt
head missing_HIGH_MODERATE.protein_ids.txt
```
