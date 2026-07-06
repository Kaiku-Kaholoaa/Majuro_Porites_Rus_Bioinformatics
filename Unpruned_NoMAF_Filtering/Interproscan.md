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
