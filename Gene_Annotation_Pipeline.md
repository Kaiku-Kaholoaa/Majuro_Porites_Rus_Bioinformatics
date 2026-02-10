# Gene Annotation of *Porites rus* 
# With Maker2, SNAP, and Augustus

## Download required databases (Swissprot)

'cat maker_opts.ctl'

```bash
#-----Genome (these are always required)
genome=inputs_for_maker2/porites_rus_reference_genome.fna
organism_type=eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff=
est_pass=0
altest_pass=0
protein_pass=0
rm_pass=0
model_pass=0
pred_pass=0
other_pass=0

#-----EST Evidence (for best results provide a file for at least one)
est=inputs_for_maker2/porites_rus_trinity.fasta
altest=inputs_for_maker2/porites_transcripts/porites_combined_altest_transcripts.fasta
est_gff=
altest_gff=

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=inputs_for_maker2/all_porites_proteins.fasta
protein_gff=

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=
rmlib=
repeat_protein=
rm_gff=
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1

#-----Gene Prediction
snaphmm=
gmhmm=
augustus_species=
fgenesh_par_file= #FGENESH parameter file
pred_gff=
model_gff=
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1
protein2genome=1
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1

#-----MAKER Behavior Options
max_dna_len=100000
min_contig=1000

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1
min_protein=20
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=1

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files

```
