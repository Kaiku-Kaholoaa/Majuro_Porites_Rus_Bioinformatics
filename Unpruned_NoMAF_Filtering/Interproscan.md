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

### Results: 
```bash
cat summary_major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv
```

augustus_masked-OZ037999.1-processed-gene-295.10-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037993.1-exonerate_protein2genome-gene-300.91-mRNA-1	PANTHER	PTHR21301	REVERSE TRANSCRIPTASE	-	-
maker-OZ037993.1-snap-gene-142.47-mRNA-1	Coils	Coil	Coil	-	-
maker-OZ037993.1-snap-gene-142.47-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037995.1-exonerate_protein2genome-gene-334.16-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037995.1-exonerate_protein2genome-gene-334.56-mRNA-1	PANTHER	PTHR21301	REVERSE TRANSCRIPTASE	-	-
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	Coils	Coil	Coil	-	-
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	FunFam	G3DSA:2.60.120.260:FF:000009	SUN domain-containing protein 1 isoform X1	-	-
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	Gene3D	G3DSA:2.60.120.260	-	-	-
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	PANTHER	PTHR12911	SAD1/UNC-84-LIKE PROTEIN-RELATED	IPR045119	SUN domain-containing protein 1-5
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	Pfam	PF07738	Sad1 / UNC-like C-terminal	IPR012919	SUN domain
maker-OZ037997.1-exonerate_protein2genome-gene-246.32-mRNA-1	SUPERFAMILY	SSF49785	Galactose-binding domain-like	IPR008979	Galactose-binding-like domain superfamily
maker-OZ037997.1-snap-gene-246.4-mRNA-1	CDD	cd03263	ABC_subfamily_A	-	-
maker-OZ037997.1-snap-gene-246.4-mRNA-1	FunFam	G3DSA:3.40.50.300:FF:000335	ATP binding cassette subfamily A member 5	-	-
maker-OZ037997.1-snap-gene-246.4-mRNA-1	FunFam	G3DSA:3.40.50.300:FF:000436	ATP binding cassette subfamily A member 9	-	-
maker-OZ037997.1-snap-gene-246.4-mRNA-1	Gene3D	G3DSA:3.40.50.300	-	IPR027417	P-loop containing nucleoside triphosphate hydrolase
maker-OZ037997.1-snap-gene-246.4-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037997.1-snap-gene-246.4-mRNA-1	PANTHER	PTHR19229	ATP-BINDING CASSETTE TRANSPORTER SUBFAMILY A  ABCA	IPR026082	ABC transporter A
maker-OZ037997.1-snap-gene-246.4-mRNA-1	Pfam	PF00005	ABC transporter	IPR003439	ABC transporter-like, ATP-binding domain
maker-OZ037997.1-snap-gene-246.4-mRNA-1	Pfam	PF12698	ABC-2 family transporter protein	-	-
maker-OZ037997.1-snap-gene-246.4-mRNA-1	SMART	SM00382	AAA_5	IPR003593	AAA+ ATPase domain
maker-OZ037997.1-snap-gene-246.4-mRNA-1	SUPERFAMILY	SSF52540	P-loop containing nucleoside triphosphate hydrolases	IPR027417	P-loop containing nucleoside triphosphate hydrolase
maker-OZ037999.1-augustus-gene-295.23-mRNA-1	Coils	Coil	Coil	-	-
maker-OZ037999.1-augustus-gene-295.23-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ037999.1-augustus-gene-295.23-mRNA-1	PANTHER	PTHR47331	PHD-TYPE DOMAIN-CONTAINING PROTEIN	-	-
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	CDD	cd01644	RT_pepA17	-	-
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Gene3D	G3DSA:1.10.340.70	-	-	-
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Gene3D	G3DSA:3.10.10.10	HIV Type 1 Reverse Transcriptase, subunit A, domain 1	-	-
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Gene3D	G3DSA:3.30.420.10	-	IPR036397	Ribonuclease H superfamily
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Gene3D	G3DSA:3.30.70.270	-	IPR043128	Reverse transcriptase/Diguanylate cyclase domain
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	PANTHER	PTHR22955	RETROTRANSPOSON	-	-
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Pfam	PF05380	Pao retrotransposon peptidase	IPR008042	Retrotransposon, Pao
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Pfam	PF17921	Integrase zinc binding domain	IPR041588	Integrase zinc-binding domain
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	Pfam	PF18701	Family of unknown function (DUF5641)	IPR040676	Domain of unknown function DUF5641
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	SUPERFAMILY	SSF53098	Ribonuclease H-like	IPR012337	Ribonuclease H-like superfamily
maker-OZ037999.1-exonerate_protein2genome-gene-295.29-mRNA-1	SUPERFAMILY	SSF56672	DNA/RNA polymerases	IPR043502	DNA/RNA polymerase superfamily
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	CDD	cd00051	EFh	IPR002048	EF-hand domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Coils	Coil	Coil	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Gene3D	G3DSA:1.10.238.10	-	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	PANTHER	PTHR10869	PROLYL 4-HYDROXYLASE ALPHA SUBUNIT	IPR045054	Prolyl 4-hydroxylase
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Pfam	PF13499	EF-hand domain pair	IPR002048	EF-hand domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Pfam	PF13640	2OG-Fe(II) oxygenase superfamily	IPR044862	Prolyl 4-hydroxylase alpha subunit, Fe(2+) 2OG dioxygenase domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	SUPERFAMILY	SSF47473	EF-hand	IPR011992	EF-hand domain pair
maker-OZ038003.1-augustus-gene-230.42-mRNA-1	CDD	cd10975	CE4_CDA_like_2	-	-
maker-OZ038003.1-augustus-gene-230.42-mRNA-1	Gene3D	G3DSA:3.20.20.370	Glycoside hydrolase/deacetylase	-	-
maker-OZ038003.1-augustus-gene-230.42-mRNA-1	PANTHER	PTHR45985	-	-	-
maker-OZ038003.1-augustus-gene-230.42-mRNA-1	Pfam	PF01522	Polysaccharide deacetylase	IPR002509	NodB homology domain
maker-OZ038003.1-augustus-gene-230.42-mRNA-1	SUPERFAMILY	SSF88713	Glycoside hydrolase/deacetylase	IPR011330	Glycoside hydrolase/deacetylase, beta/alpha-barrel
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Gene3D	G3DSA:1.10.238.10	-	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	PANTHER	PTHR10869	PROLYL 4-HYDROXYLASE ALPHA SUBUNIT	IPR045054	Prolyl 4-hydroxylase
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Pfam	PF13640	2OG-Fe(II) oxygenase superfamily	IPR044862	Prolyl 4-hydroxylase alpha subunit, Fe(2+) 2OG dioxygenase domain
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	SMART	SM00702	p4hc	IPR006620	Prolyl 4-hydroxylase, alpha subunit
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	SUPERFAMILY	SSF47473	EF-hand	IPR011992	EF-hand domain pair
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	CDD	cd00051	EFh	IPR002048	EF-hand domain
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	Gene3D	G3DSA:1.10.238.10	-	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	PANTHER	PTHR10869	PROLYL 4-HYDROXYLASE ALPHA SUBUNIT	IPR045054	Prolyl 4-hydroxylase
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	Pfam	PF13640	2OG-Fe(II) oxygenase superfamily	IPR044862	Prolyl 4-hydroxylase alpha subunit, Fe(2+) 2OG dioxygenase domain
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	SMART	SM00702	p4hc	IPR006620	Prolyl 4-hydroxylase, alpha subunit
maker-OZ038003.1-exonerate_est2genome-gene-230.11-mRNA-1	SUPERFAMILY	SSF47473	EF-hand	IPR011992	EF-hand domain pair
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	CDD	cd01723	LSm4	IPR034101	Sm-like protein Lsm4
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	Gene3D	G3DSA:2.30.30.100	-	-	-
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	PANTHER	PTHR23338	SMALL NUCLEAR RIBONUCLEOPROTEIN SM	IPR027141	Like-Sm (LSM) domain containing protein, LSm4/SmD1/SmD3
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	Pfam	PF01423	LSM domain	IPR001163	LSM domain, eukaryotic/archaea-type
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	SMART	SM00651	Sm3	IPR001163	LSM domain, eukaryotic/archaea-type
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	SUPERFAMILY	SSF50182	Sm-like ribonucleoproteins	IPR010920	LSM domain superfamily
maker-OZ038003.1-snap-gene-230.15-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-snap-gene-230.15-mRNA-1	PANTHER	PTHR12117	HISTONE ACETYLTRANSFERASE COMPLEX	-	-
maker-OZ038003.1-snap-gene-230.15-mRNA-1	Pfam	PF10637	Oxoglutarate and iron-dependent oxygenase degradation C-term	IPR019601	Oxoglutarate/iron-dependent oxygenase, C-terminal degradation domain
maker-OZ038003.1-snap-gene-230.15-mRNA-1	Pfam	PF13661	2OG-Fe(II) oxygenase superfamily	IPR039558	Prolyl 3,4-dihydroxylase TPA1/OFD1, N-terminal domain
maker-OZ038003.1-snap-gene-230.15-mRNA-1	SMART	SM00702	p4hc	IPR006620	Prolyl 4-hydroxylase, alpha subunit
maker-OZ038003.1-snap-gene-230.52-mRNA-1	CDD	cd10975	CE4_CDA_like_2	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	Gene3D	G3DSA:3.20.20.370	Glycoside hydrolase/deacetylase	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	PANTHER	PTHR45985	-	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	Pfam	PF01522	Polysaccharide deacetylase	IPR002509	NodB homology domain
maker-OZ038003.1-snap-gene-230.52-mRNA-1	SUPERFAMILY	SSF88713	Glycoside hydrolase/deacetylase	IPR011330	Glycoside hydrolase/deacetylase, beta/alpha-barrel
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:2.60.120.260	-	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:3.30.2410.10	Hect, E3 ligase catalytic domain	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:3.90.1750.10	Hect, E3 ligase catalytic domains	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	PANTHER	PTHR46654	E3 UBIQUITIN-PROTEIN LIGASE HECTD3	IPR042469	E3 ubiquitin-protein ligase HECTD3
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Pfam	PF00632	HECT-domain (ubiquitin-transferase)	IPR000569	HECT domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Pfam	PF03256	Anaphase-promoting complex, subunit 10 (APC10)	IPR004939	APC10/DOC domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SMART	SM00119	hect_3	IPR000569	HECT domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SMART	SM01337	APC10_2	IPR004939	APC10/DOC domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SUPERFAMILY	SSF49785	Galactose-binding domain-like	IPR008979	Galactose-binding-like domain superfamily
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SUPERFAMILY	SSF56204	Hect, E3 ligase catalytic domain	IPR035983	HECT, E3 ligase catalytic domain
maker-OZ038003.1-snap-gene-230.57-mRNA-1	CDD	cd00637	7tm_classA_rhodopsin-like	-	-
maker-OZ038003.1-snap-gene-230.57-mRNA-1	Gene3D	G3DSA:1.20.1070.10	-	-	-
maker-OZ038003.1-snap-gene-230.57-mRNA-1	PANTHER	PTHR24240	OPSIN	-	-
maker-OZ038003.1-snap-gene-230.57-mRNA-1	Pfam	PF00001	7 transmembrane receptor (rhodopsin family)	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-230.57-mRNA-1	PRINTS	PR00237	Rhodopsin-like GPCR superfamily signature	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-230.57-mRNA-1	SMART	SM01381	7TM_GPCR_Srsx_2	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-230.57-mRNA-1	SUPERFAMILY	SSF81321	Family A G protein-coupled receptor-like	-	-
maker-OZ038003.1-snap-gene-231.7-mRNA-1	CDD	cd00637	7tm_classA_rhodopsin-like	-	-
maker-OZ038003.1-snap-gene-231.7-mRNA-1	Gene3D	G3DSA:1.20.1070.10	-	-	-
maker-OZ038003.1-snap-gene-231.7-mRNA-1	PANTHER	PTHR24240	OPSIN	-	-
maker-OZ038003.1-snap-gene-231.7-mRNA-1	Pfam	PF00001	7 transmembrane receptor (rhodopsin family)	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-231.7-mRNA-1	PRINTS	PR00237	Rhodopsin-like GPCR superfamily signature	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-231.7-mRNA-1	SMART	SM01381	7TM_GPCR_Srsx_2	IPR000276	G protein-coupled receptor, rhodopsin-like
maker-OZ038003.1-snap-gene-231.7-mRNA-1	SUPERFAMILY	SSF81321	Family A G protein-coupled receptor-like	-	-

### Taking a closer look at the 5 genes from the high-impact variants (no moderate-variants)

```bash
grep -F -f major_allele_differences_NSNMU.HIGH.gene_ids.txt summary_major_allele_differences_NSNMU.HIGH_MODERATE.interpro.tsv > high_only_summary_major_allele_differences_NSNMU.HIGH.interpro.tsv

wc -l high_only_summary_major_allele_differences_NSNMU.HIGH.interpro.tsv 
#38 interpro annotations
cut -f1 high_only_summary_major_allele_differences_NSNMU.HIGH.interpro.tsv | sort -u | wc -l
#5 unique genes
```

```bash
cat high_only_summary_major_allele_differences_NSNMU.HIGH.interpro.tsv
```
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	CDD	cd00051	EFh	IPR002048	EF-hand domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Coils	Coil	Coil	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Gene3D	G3DSA:1.10.238.10	-	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	PANTHER	PTHR10869	PROLYL 4-HYDROXYLASE ALPHA SUBUNIT	IPR045054	Prolyl 4-hydroxylase
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Pfam	PF13499	EF-hand domain pair	IPR002048	EF-hand domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	Pfam	PF13640	2OG-Fe(II) oxygenase superfamily	IPR044862	Prolyl 4-hydroxylase alpha subunit, Fe(2+) 2OG dioxygenase domain
maker-OZ038003.1-augustus-gene-230.19-mRNA-1	SUPERFAMILY	SSF47473	EF-hand	IPR011992	EF-hand domain pair
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Gene3D	G3DSA:1.10.238.10	-	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Gene3D	G3DSA:2.60.120.620	q2cbj1_9rhob like domain	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	PANTHER	PTHR10869	PROLYL 4-HYDROXYLASE ALPHA SUBUNIT	IPR045054	Prolyl 4-hydroxylase
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	Pfam	PF13640	2OG-Fe(II) oxygenase superfamily	IPR044862	Prolyl 4-hydroxylase alpha subunit, Fe(2+) 2OG dioxygenase domain
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	SMART	SM00702	p4hc	IPR006620	Prolyl 4-hydroxylase, alpha subunit
maker-OZ038003.1-exonerate_est2genome-gene-230.10-mRNA-1	SUPERFAMILY	SSF47473	EF-hand	IPR011992	EF-hand domain pair
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	CDD	cd01723	LSm4	IPR034101	Sm-like protein Lsm4
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	Gene3D	G3DSA:2.30.30.100	-	-	-
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	MobiDBLite	mobidb-lite	consensus disorder prediction	-	-
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	PANTHER	PTHR23338	SMALL NUCLEAR RIBONUCLEOPROTEIN SM	IPR027141	Like-Sm (LSM) domain containing protein, LSm4/SmD1/SmD3
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	Pfam	PF01423	LSM domain	IPR001163	LSM domain, eukaryotic/archaea-type
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	SMART	SM00651	Sm3	IPR001163	LSM domain, eukaryotic/archaea-type
maker-OZ038003.1-exonerate_protein2genome-gene-230.49-mRNA-1	SUPERFAMILY	SSF50182	Sm-like ribonucleoproteins	IPR010920	LSM domain superfamily
maker-OZ038003.1-snap-gene-230.52-mRNA-1	CDD	cd10975	CE4_CDA_like_2	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	Gene3D	G3DSA:3.20.20.370	Glycoside hydrolase/deacetylase	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	PANTHER	PTHR45985	-	-	-
maker-OZ038003.1-snap-gene-230.52-mRNA-1	Pfam	PF01522	Polysaccharide deacetylase	IPR002509	NodB homology domain
maker-OZ038003.1-snap-gene-230.52-mRNA-1	SUPERFAMILY	SSF88713	Glycoside hydrolase/deacetylase	IPR011330	Glycoside hydrolase/deacetylase, beta/alpha-barrel
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:2.60.120.260	-	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:3.30.2410.10	Hect, E3 ligase catalytic domain	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Gene3D	G3DSA:3.90.1750.10	Hect, E3 ligase catalytic domains	-	-
maker-OZ038003.1-snap-gene-230.53-mRNA-1	PANTHER	PTHR46654	E3 UBIQUITIN-PROTEIN LIGASE HECTD3	IPR042469	E3 ubiquitin-protein ligase HECTD3
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Pfam	PF00632	HECT-domain (ubiquitin-transferase)	IPR000569	HECT domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	Pfam	PF03256	Anaphase-promoting complex, subunit 10 (APC10)	IPR004939	APC10/DOC domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SMART	SM00119	hect_3	IPR000569	HECT domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SMART	SM01337	APC10_2	IPR004939	APC10/DOC domain
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SUPERFAMILY	SSF49785	Galactose-binding domain-like	IPR008979	Galactose-binding-like domain superfamily
maker-OZ038003.1-snap-gene-230.53-mRNA-1	SUPERFAMILY	SSF56204	Hect, E3 ligase catalytic domain	IPR035983	HECT, E3 ligase catalytic domain

