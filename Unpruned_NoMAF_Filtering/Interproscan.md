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
#10 unique genes.

#Wow, 10 of our genes are within this haplotype block(?) on that chr 12 region, which is especially interesting since the 8 high-variant alleles lie within only 5 genes. 
```
Summarize what databases/domains are represented:

```bash
cut -f4 major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv \
  | sort \
  | uniq -c \
  | sort -nr
```

Summarize the most common descriptions:

```bash
cut -f4,5,6 major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv \
  | sort \
  | uniq -c \
  | sort -nr
```

Finally, summarize candidate transcripts with InterPro descriptions (important!)

```bash
cut -f1,4,5,6,12,13 major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv \
  | sort -u \
  > summary_major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv
```


### Taking a closer look at the 5 genes from the high-impact variants (no moderate-variants)

```bash
grep -F -f major_allele_differences_NSNMU.HIGH.gene_ids.txt major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv > interpro_high_impact_variants_only.tsv
```


