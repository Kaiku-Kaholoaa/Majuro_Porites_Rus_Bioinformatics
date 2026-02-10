# Functional Annotation of *Porites rus*

### Functional Annotation Pipeline Adopted From Researcher Kevin Wong: https://kevinhwong1.github.io/KevinHWong_Notebook/Porites-astreoides-genome-annotation/ 
## Download required databases (Swissprot)

`cat build_swissprot_db.sh`

```bash
#!/bin/bash
#SBATCH --job-name=build_swissprot_db
#SBATCH -p hns,spalumbi,serc
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --output=build_swissprot_db.out
#SBATCH --error=build_swissprot_db.err

set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# ============================
# CONFIGURATION
# ============================

DB_DIR=/scratch/users/kaiku/databases/swissprot
UNIPROT_URL=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
FASTA_GZ=uniprot_sprot.fasta.gz
FASTA=uniprot_sprot.fasta
DB_NAME=swissprot

# ============================
# LOAD MODULES
# ============================

ml ncbi-blast+/2.16.0

# ============================
# SETUP DIRECTORY
# ============================

mkdir -p ${DB_DIR}
cd ${DB_DIR}

echo "Working directory: $(pwd)"

# ============================
# DOWNLOAD SWISS-PROT
# ============================

if [ ! -f ${FASTA_GZ} ]; then
    echo "Downloading Swiss-Prot FASTA from UniProt..."
    wget ${UNIPROT_URL}
else
    echo "Swiss-Prot FASTA archive already exists, skipping download."
fi

# ============================
# UNZIP FASTA
# ============================

if [ ! -f ${FASTA} ]; then
    echo "Unzipping Swiss-Prot FASTA..."
    gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta
else
    echo "Uncompressed FASTA already exists, skipping unzip."
fi

# ============================
# BUILD BLAST DATABASE
# ============================

echo "Building BLAST database..."
makeblastdb \
  -in ${FASTA} \
  -dbtype prot \
  -parse_seqids \
  -out ${DB_NAME}

# ============================
# VERIFY DATABASE
# ============================

echo "Verifying BLAST database..."
blastdbcmd -db ${DB_NAME} -info

echo "Swiss-Prot BLAST database successfully built."
echo "Job finished at $(date)"
```

## Align query protein sequences against databases

#### 1. BLAST the protein sequences against Swiss-Prot (Highest Quality Protein Matches)

`nano swiss_prot.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=SwPrtBlast
#SBATCH -p hns,spalumbi,serc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -t 168:00:00
#SBATCH --mem=100GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu
#SBATCH --output=swiss_prot_out
#SBATCH --error=swiss_prot_err

# Load BUSCO

ml biology ncbi-blast+/2.16.0

echo "Blast against swissprot database" $(date)

blastp -max_target_seqs 5 \
-num_threads 20 \
-db /scratch/users/kaiku/databases/swissprot/swissprot \
-query no_mpi_round3.2.all.maker.proteins.busco.fasta \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_sprot_1e-5_max5.out


echo "STOP" $(date)
```

i) Get the best hit for each Gene Model (protein) Swiss-Prot

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PastGeneModels_vs_sprot_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_sprot_1e-5_besthit.out

wc -l PastGeneModels_vs_sprot_1e-5_max5.out
#202,765 homologous protein matches

wc -l PastGeneModels_vs_sprot_1e-5_besthit.out
#23,678 best hits
```

ii) Select the gene model proteins without hits in Swiss-Prot


first use awk to print a list of all the Gene Model names from besthits.out
```
awk '{print $1}' PastGeneModels_vs_sprot_1e-5_besthit.out > list_of_Pastgenemodelproteins_sprot.txt

wc -l list_of_Pastgenemodelproteins_sprot.txt #30,487
```

Then exclude these Gene Model names from your original fasta/.faa/protein file.
* needed to load the module that has the script with the -exclude command in it

#first, loaded the newest module for kentUtils/416-foss-2020b

`module load kentUtils/416-foss-2020b`


I selected to the prepend-path /opt/software/kentUtils/416-foss-2020b/bin to see if it took me to the 'faSomeRecords' script which it did

`/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords`

I then ran the -exclude command to exclude the blasted Gene Models from the .faa file

```
/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude {PATH}/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta list_of_Pastgenemodelproteins_sprot.txt Past_proteins_names_v1.0.faa.prot4trembl
```

Checking the number of Gene Models:

`grep -c ">" Past_proteins_names_v1.0.faa.prot4trembl #34,149`

* use this file to blast against trembl

Downloading the .xml file for Blast2Go

`nano swissprot_blast_xml.sh`

```bash
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein-xml"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_error"
#SBATCH --output="xml_blastp_out"
#SBATCH --exclusive
#SBATCH -D {PATH}/past_struc_annotations_v1/functional_anno_v1  

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database with xml format out" $(date)
blastp -max_target_seqs 5 \
-num_threads 20 \
-db {PATH}/databases/swiss_db/swissprot_20211022 \
-query {PATH}/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_maxhit.xml

echo "STOP" $(date)
```

#### 10. BLAST the remaining protein sequences against Trembl

`nano trembl_blastp.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="trembl_blastp_out_error"
#SBATCH --output="trembl_blastp_out"
#SBATCH --exclusive
#SBATCH -D {PATH}/past_struc_annotations_v1/functional_anno_v1

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database" $(date)
blastp -max_target_seqs 5 \
-num_threads 20 \
-db {PATH}/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_trembl_1e-5_max5.out

echo "STOP" $(date)
```

`nano trembl_blastp_hit1.sh`

```bash
#!/bin/bash
#SBATCH --job-name="trembl-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="trembl_hit1_blastp_out_error"
#SBATCH --output="trembl_hit1_blastp_out"
#SBATCH --exclusive
#SBATCH -D {PATH}/past_struc_annotations_v1/functional_anno_v1
#SBATCH -c 36

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database" $(date)
blastp -max_target_seqs 1 \
-num_threads 20 \
-db {PATH}/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_vs_trembl_1e-5_max1.out

echo "STOP" $(date)
```

Get the best hit for each Gene Model (protein) Trembl
```
cat PastGeneModels_vs_trembl_1e-5_max1.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_trembl_1e-5_besthit.out

wc -l PastGeneModels_vs_trembl_1e-5_besthit.out #24525
```

Select the gene model proteins without hits in Trembl
```bash
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PastGeneModels_vs_trembl_1e-5_besthit.out > list_of_Pastgenemodelproteins_trembl.txt

#load the newest module for kentUtils/416-foss-2020b

module load kentUtils/416-foss-2020b

#then exclude these Gene Model names from your original fasta/.faa/protein file

/opt/software/kentUtils/416-foss-2020b/bin/faSomeRecords -exclude Past_proteins_names_v1.0.faa.prot4trembl list_of_Pastgenemodelproteins_trembl.txt Past_proteins_names_v1.0.faa.prot4nr


#check the number of Gene Models

grep -c ">" Past_proteins_names_v1.0.faa.prot4nr #9624

#using this file to blast against nr database

```

`nano trembl_blastp_xml.sh`

```bash
#!/bin/bash
#SBATCH --job-name="xml-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="xml_blastp_out_terror"
#SBATCH --output="xml_blastp_tout"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against trembl database for xml" $(date)
blastp -max_target_seqs 1 \
-num_threads 20 \
-db {PATH}/databases/trembl_db/trembl_20211022 \
-query Past_proteins_names_v1.0.faa.prot4trembl \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out Past_protein_blastp_trembl.xml

echo "STOP" $(date)
```

#### 11. BLAST the remaining protein sequences against nr

`nano ncbi_blastp_out.sh`

```bash
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein-out"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/shared/ncbi-nr/nr \
-query Past_proteins_names_v1.0.faa.prot4nr \
-evalue 1e-5 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_ncbi_max5hits.out

echo "STOP" $(date)
```


`nano ncbi_blastp.sh`

```bash
#!/bin/bash
#SBATCH --job-name="ncbi-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="ncbi_blastp_out_error"
#SBATCH --output="ncbi_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against ncbi database" $(date)
blastp -max_target_seqs 5 \
-num_threads $SLURM_CPUS_ON_NODE \
-db /data/shared/ncbi-nr/nr \
-query Past_proteins_names_v1.0.faa.prot4nr \
-evalue 1e-5 \
-outfmt '5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
-out PastGeneModels_ncbi.xml

echo "STOP" $(date)
```


#### 12. Interproscan

`sbatch Past_InterProScan.sh`

```bash
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=100GB
#SBATCH --error="interproscan_out_error"
#SBATCH --output="interproscan_out"
#SBATCH --exclusive
#SBATCH -D {PATH}/past_struc_annotations_v1/functional_anno_v1/InterProScan

echo "START $(date)"

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...
interproscan.sh -version
interproscan.sh -f XML -i {PATH}/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta -b Past.interpro.20220113  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i Past.interpro.20220113.xml -b Past.interpro.20220113

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"
```


#### 13. Statistics with AGAT

#### Porites lutea

`module load AGAT/0.8.1-foss-2020b`

`agat_sp_statistics.pl --gff plut2v1.1.genes.gff3`


```bash
Reading file plut2v1.1.genes.gff3
Parsing: 100% [======================================================]D 0h01m12sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              31126
Number of mrnas                              31126
Number of mrnas with utr both sides          9305
Number of mrnas with at least one utr        15997
Number of cdss                               31126
Number of exons                              204364
Number of five_prime_utrs                    12044
Number of three_prime_utrs                   13258
Number of exon in cds                        197258
Number of exon in five_prime_utr             16114
Number of exon in three_prime_utr            16195
Number of intron in cds                      166132
Number of intron in exon                     173238
Number of intron in five_prime_utr           4070
Number of intron in three_prime_utr          2937
Number gene overlapping                      1069
Number of single exon gene                   5296
Number of single exon mrna                   5296
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          6.6
mean five_prime_utrs per mrna                0.4
mean three_prime_utrs per mrna               0.4
mean exons per cds                           6.3
mean exons per five_prime_utr                1.3
mean exons per three_prime_utr               1.2
mean introns in cdss per mrna                5.3
mean introns in exons per mrna               5.6
mean introns in five_prime_utrs per mrna     0.1
mean introns in three_prime_utrs per mrna    0.1
Total gene length                            255612213
Total mrna length                            255612213
Total cds length                             43135478
Total exon length                            56249572
Total five_prime_utr length                  2866380
Total three_prime_utr length                 10247714
Total intron length per cds                  190476728
Total intron length per exon                 199362641
Total intron length per five_prime_utr       6152481
Total intron length per three_prime_utr      2502850
mean gene length                             8212
mean mrna length                             8212
mean cds length                              1385
mean exon length                             275
mean five_prime_utr length                   237
mean three_prime_utr length                  772
mean cds piece length                        218
mean five_prime_utr piece length             177
mean three_prime_utr piece length            632
mean intron in cds length                    1146
mean intron in exon length                   1150
mean intron in five_prime_utr length         1511
mean intron in three_prime_utr length        852
Longest gene                                 143011
Longest mrna                                 143011
Longest cds                                  52230
Longest exon                                 12544
Longest five_prime_utr                       12253
Longest three_prime_utr                      11162
Longest cds piece                            12180
Longest five_prime_utr piece                 12253
Longest three_prime_utr piece                6747
Longest intron into cds part                 37855
Longest intron into exon part                37855
Longest intron into five_prime_utr part      28545
Longest intron into three_prime_utr part     13676
Shortest gene                                159
Shortest mrna                                159
Shortest cds                                 9
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     20
Shortest intron into three_prime_utr part    20

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              31126
Number of mrnas                              31126
Number of mrnas with utr both sides          9305
Number of mrnas with at least one utr        15997
Number of cdss                               31126
Number of exons                              204364
Number of five_prime_utrs                    12044
Number of three_prime_utrs                   13258
Number of exon in cds                        197258
Number of exon in five_prime_utr             16114
Number of exon in three_prime_utr            16195
Number of intron in cds                      166132
Number of intron in exon                     173238
Number of intron in five_prime_utr           4070
Number of intron in three_prime_utr          2937
Number gene overlapping                      1069
Number of single exon gene                   5296
Number of single exon mrna                   5296
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          6.6
mean five_prime_utrs per mrna                0.4
mean three_prime_utrs per mrna               0.4
mean exons per cds                           6.3
mean exons per five_prime_utr                1.3
mean exons per three_prime_utr               1.2
mean introns in cdss per mrna                5.3
mean introns in exons per mrna               5.6
mean introns in five_prime_utrs per mrna     0.1
mean introns in three_prime_utrs per mrna    0.1
Total gene length                            255612213
Total mrna length                            255612213
Total cds length                             43135478
Total exon length                            56249572
Total five_prime_utr length                  2866380
Total three_prime_utr length                 10247714
Total intron length per cds                  190476728
Total intron length per exon                 199362641
Total intron length per five_prime_utr       6152481
Total intron length per three_prime_utr      2502850
mean gene length                             8212
mean mrna length                             8212
mean cds length                              1385
mean exon length                             275
mean five_prime_utr length                   237
mean three_prime_utr length                  772
mean cds piece length                        218
mean five_prime_utr piece length             177
mean three_prime_utr piece length            632
mean intron in cds length                    1146
mean intron in exon length                   1150
mean intron in five_prime_utr length         1511
mean intron in three_prime_utr length        852
Longest gene                                 143011
Longest mrna                                 143011
Longest cds                                  52230
Longest exon                                 12544
Longest five_prime_utr                       12253
Longest three_prime_utr                      11162
Longest cds piece                            12180
Longest five_prime_utr piece                 12253
Longest three_prime_utr piece                6747
Longest intron into cds part                 37855
Longest intron into exon part                37855
Longest intron into five_prime_utr part      28545
Longest intron into three_prime_utr part     13676
Shortest gene                                159
Shortest mrna                                159
Shortest cds                                 9
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     20
Shortest intron into three_prime_utr part    20

--------------------------------------------------------------------------------

Bye Bye.
```


#### Porites astreoides

`agat_sp_statistics.pl --gff Pastreoides_all_v1.gff`

```bash
Reading file Pastreoides_all_v1.gff
Parsing: 100% [======================================================]D 0h10m33sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute match_part with isoforms if any

Number of matchs                             905337
Number of protein_matchs                     300777
Number of introns                            1062102
Number of match_parts                        2268216
Number gene overlapping                      0
Number of single exon match                  905337
Number of single exon protein_match          300777
mean introns per match                       1.2
mean introns per protein_match               3.5
mean match_parts per match                   2.5
mean match_parts per protein_match           7.5
Total match length                           939607611
Total protein_match length                   676247819
Total intron length                          1166327186
Total match_part length                      449528244
mean match length                            1037
mean protein_match length                    2248
mean intron length                           1098
mean match_part length                       198
Longest match                                104218
Longest protein_match                        140218
Longest intron                               46623
Longest match_part                           29944
Shortest match                               6
Shortest protein_match                       66
Shortest intron                              4
Shortest match_part                          2

Re-compute match_part without isoforms asked. We remove shortest isoforms if any

Number of matchs                             905337
Number of protein_matchs                     300777
Number of introns                            1062102
Number of match_parts                        2268216
Number gene overlapping                      0
Number of single exon match                  905337
Number of single exon protein_match          300777
mean introns per match                       1.2
mean introns per protein_match               3.5
mean match_parts per match                   2.5
mean match_parts per protein_match           7.5
Total match length                           939607611
Total protein_match length                   676247819
Total intron length                          1166327186
Total match_part length                      449528244
mean match length                            1037
mean protein_match length                    2248
mean intron length                           1098
mean match_part length                       198
Longest match                                104218
Longest protein_match                        140218
Longest intron                               46623
Longest match_part                           29944
Shortest match                               6
Shortest protein_match                       66
Shortest intron                              4
Shortest match_part                          2

--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              64636
Number of mrnas                              64636
Number of mrnas with utr both sides          1623
Number of mrnas with at least one utr        8185
Number of cdss                               64636
Number of exons                              317304
Number of five_prime_utrs                    5227
Number of three_prime_utrs                   4581
Number of exon in cds                        315399
Number of exon in five_prime_utr             6439
Number of exon in three_prime_utr            5240
Number of intron in cds                      250763
Number of intron in exon                     252668
Number of intron in five_prime_utr           1212
Number of intron in three_prime_utr          659
Number gene overlapping                      1149
Number of single exon gene                   6221
Number of single exon mrna                   6221
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          4.9
mean five_prime_utrs per mrna                0.1
mean three_prime_utrs per mrna               0.1
mean exons per cds                           4.9
mean exons per five_prime_utr                1.2
mean exons per three_prime_utr               1.1
mean introns in cdss per mrna                3.9
mean introns in exons per mrna               3.9
mean introns in five_prime_utrs per mrna     0.0
mean introns in three_prime_utrs per mrna    0.0
Total gene length                            279311982
Total mrna length                            279311982
Total cds length                             58377165
Total exon length                            60411427
Total five_prime_utr length                  564947
Total three_prime_utr length                 1469315
Total intron length per cds                  216586326
Total intron length per exon                 218900555
Total intron length per five_prime_utr       1735255
Total intron length per three_prime_utr      537735
mean gene length                             4321
mean mrna length                             4321
mean cds length                              903
mean exon length                             190
mean five_prime_utr length                   108
mean three_prime_utr length                  320
mean cds piece length                        185
mean five_prime_utr piece length             87
mean three_prime_utr piece length            280
mean intron in cds length                    863
mean intron in exon length                   866
mean intron in five_prime_utr length         1431
mean intron in three_prime_utr length        815
Longest gene                                 105022
Longest mrna                                 105022
Longest cds                                  26259
Longest exon                                 10349
Longest five_prime_utr                       1687
Longest three_prime_utr                      2206
Longest cds piece                            10349
Longest five_prime_utr piece                 1501
Longest three_prime_utr piece                2206
Longest intron into cds part                 32611
Longest intron into exon part                32611
Longest intron into five_prime_utr part      11040
Longest intron into three_prime_utr part     7173
Shortest gene                                39
Shortest mrna                                39
Shortest cds                                 9
Shortest exon                                1
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     5
Shortest intron into three_prime_utr part    5

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              64636
Number of mrnas                              64636
Number of mrnas with utr both sides          1623
Number of mrnas with at least one utr        8185
Number of cdss                               64636
Number of exons                              317304
Number of five_prime_utrs                    5227
Number of three_prime_utrs                   4581
Number of exon in cds                        315399
Number of exon in five_prime_utr             6439
Number of exon in three_prime_utr            5240
Number of intron in cds                      250763
Number of intron in exon                     252668
Number of intron in five_prime_utr           1212
Number of intron in three_prime_utr          659
Number gene overlapping                      1149
Number of single exon gene                   6221
Number of single exon mrna                   6221
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          4.9
mean five_prime_utrs per mrna                0.1
mean three_prime_utrs per mrna               0.1
mean exons per cds                           4.9
mean exons per five_prime_utr                1.2
mean exons per three_prime_utr               1.1
mean introns in cdss per mrna                3.9
mean introns in exons per mrna               3.9
mean introns in five_prime_utrs per mrna     0.0
mean introns in three_prime_utrs per mrna    0.0
Total gene length                            279311982
Total mrna length                            279311982
Total cds length                             58377165
Total exon length                            60411427
Total five_prime_utr length                  564947
Total three_prime_utr length                 1469315
Total intron length per cds                  216586326
Total intron length per exon                 218900555
Total intron length per five_prime_utr       1735255
Total intron length per three_prime_utr      537735
mean gene length                             4321
mean mrna length                             4321
mean cds length                              903
mean exon length                             190
mean five_prime_utr length                   108
mean three_prime_utr length                  320
mean cds piece length                        185
mean five_prime_utr piece length             87
mean three_prime_utr piece length            280
mean intron in cds length                    863
mean intron in exon length                   866
mean intron in five_prime_utr length         1431
mean intron in three_prime_utr length        815
Longest gene                                 105022
Longest mrna                                 105022
Longest cds                                  26259
Longest exon                                 10349
Longest five_prime_utr                       1687
Longest three_prime_utr                      2206
Longest cds piece                            10349
Longest five_prime_utr piece                 1501
Longest three_prime_utr piece                2206
Longest intron into cds part                 32611
Longest intron into exon part                32611
Longest intron into five_prime_utr part      11040
Longest intron into three_prime_utr part     7173
Shortest gene                                39
Shortest mrna                                39
Shortest cds                                 9
Shortest exon                                1
Shortest five_prime_utr                      1
Shortest three_prime_utr                     1
Shortest cds piece                           1
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                4
Shortest intron into exon part               4
Shortest intron into five_prime_utr part     5
Shortest intron into three_prime_utr part    5

--------------------------------------------------------------------------------

Bye Bye.
```


#### Porites rus

`agat_sp_statistics.pl --gff Prus_gene_prediction.gff`

```bash
Reading file Prus_gene_prediction.gff
Parsing: 100% [======================================================]D 0h01m13sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute mrna with isoforms if any

Number of genes                              39453
Number of mrnas                              39453
Number of cdss                               39453
Number of exons                              195533
Number of exon in cds                        195533
Number of intron in cds                      156080
Number of intron in exon                     156080
Number gene overlapping                      0
Number of single exon gene                   10066
Number of single exon mrna                   10066
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.0
mean exons per cds                           5.0
mean introns in cdss per mrna                4.0
mean introns in exons per mrna               4.0
Total gene length                            237912951
Total mrna length                            237912951
Total cds length                             46736451
Total exon length                            46736451
Total intron length per cds                  191176500
Total intron length per exon                 191176500
mean gene length                             6030
mean mrna length                             6030
mean cds length                              1184
mean exon length                             239
mean cds piece length                        239
mean intron in cds length                    1224
mean intron in exon length                   1224
Longest gene                                 149486
Longest mrna                                 149486
Longest cds                                  36834
Longest exon                                 9978
Longest cds piece                            9978
Longest intron into cds part                 10000
Longest intron into exon part                10000
Shortest gene                                300
Shortest mrna                                300
Shortest cds                                 300
Shortest exon                                3
Shortest cds piece                           3
Shortest intron into cds part                20
Shortest intron into exon part               20

Re-compute mrna without isoforms asked. We remove shortest isoforms if any

Number of genes                              39453
Number of mrnas                              39453
Number of cdss                               39453
Number of exons                              195533
Number of exon in cds                        195533
Number of intron in cds                      156080
Number of intron in exon                     156080
Number gene overlapping                      0
Number of single exon gene                   10066
Number of single exon mrna                   10066
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.0
mean exons per cds                           5.0
mean introns in cdss per mrna                4.0
mean introns in exons per mrna               4.0
Total gene length                            237912951
Total mrna length                            237912951
Total cds length                             46736451
Total exon length                            46736451
Total intron length per cds                  191176500
Total intron length per exon                 191176500
mean gene length                             6030
mean mrna length                             6030
mean cds length                              1184
mean exon length                             239
mean cds piece length                        239
mean intron in cds length                    1224
mean intron in exon length                   1224
Longest gene                                 149486
Longest mrna                                 149486
Longest cds                                  36834
Longest exon                                 9978
Longest cds piece                            9978
Longest intron into cds part                 10000
Longest intron into exon part                10000
Shortest gene                                300
Shortest mrna                                300
Shortest cds                                 300
Shortest exon                                3
Shortest cds piece                           3
Shortest intron into cds part                20
Shortest intron into exon part               20

--------------------------------------------------------------------------------
```


#### Porites australiensis

`agat_sp_statistics.pl --gff paus_mRNA.gff`

```bash
Reading file paus_mRNA.gff
Parsing: 100% [======================================================]D 0h03m03sParsing Finished
Compute statistics
--------------------------------------------------------------------------------

Compute transcript with isoforms if any

Number of genes                              30603
Number of transcripts                        35910
Number of mrnas with utr both sides          35517
Number of mrnas with at least one utr        35907
Number of cdss                               35910
Number of exons                              314893
Number of five_prime_utrs                    35681
Number of introns                            252017
Number of start_codons                       34552
Number of stop_codons                        34675
Number of three_prime_utrs                   35743
Number of transcription_end_sites            34542
Number of exon in cds                        285689
Number of exon in five_prime_utr             56405
Number of exon in three_prime_utr            43838
Number of intron in cds                      249779
Number of intron in exon                     278983
Number of intron in five_prime_utr           20724
Number of intron in intron                   222367
Number of intron in three_prime_utr          8095
Number gene overlapping                      593
Number of single exon gene                   3970
Number of single exon transcript             4021
mean transcripts per gene                    1.2
mean cdss per transcript                     1.0
mean exons per transcript                    8.8
mean five_prime_utrs per transcript          1.0
mean introns per transcript                  7.0
mean start_codons per transcript             1.0
mean stop_codons per transcript              1.0
mean three_prime_utrs per transcript         1.0
mean transcription_end_sites per transcript  1.0
mean exons per cds                           8.0
mean exons per five_prime_utr                1.6
mean exons per three_prime_utr               1.2
mean introns in cdss per transcript          7.0
mean introns in exons per transcript         7.8
mean introns in five_prime_utrs per transcript0.6
mean introns in introns per transcript       6.2
mean introns in three_prime_utrs per transcript0.2
Total gene length                            348974335
Total transcript length                      463448447
Total cds length                             60898458
Total exon length                            99450105
Total five_prime_utr length                  10676138
Total intron length                          301789183
Total start_codon length                     103656
Total stop_codon length                      104025
Total three_prime_utr length                 27875509
Total transcription_end_site length          34542
Total intron length per cds                  297498344
Total intron length per exon                 363998342
Total intron length per five_prime_utr       49561106
Total intron length per intron               37321291
Total intron length per three_prime_utr      15710850
mean gene length                             11403
mean transcript length                       12905
mean cds length                              1695
mean exon length                             315
mean five_prime_utr length                   299
mean intron length                           1197
mean start_codon length                      3
mean stop_codon length                       3
mean three_prime_utr length                  779
mean transcription_end_site length           1
mean cds piece length                        213
mean five_prime_utr piece length             189
mean three_prime_utr piece length            635
mean intron in cds length                    1191
mean intron in exon length                   1304
mean intron in five_prime_utr length         2391
mean intron in intron length                 167
mean intron in three_prime_utr length        1940
Longest gene                                 232931
Longest transcript                           232931
Longest cds                                  63957
Longest exon                                 18417
Longest five_prime_utr                       16693
Longest intron                               87093
Longest start_codon                          3
Longest stop_codon                           3
Longest three_prime_utr                      18129
Longest transcription_end_site               1
Longest cds piece                            14759
Longest five_prime_utr piece                 16693
Longest three_prime_utr piece                18129
Longest intron into cds part                 87093
Longest intron into exon part                87093
Longest intron into five_prime_utr part      57732
Longest intron into intron part              14759
Longest intron into three_prime_utr part     42549
Shortest gene                                174
Shortest transcript                          174
Shortest cds                                 5
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest intron                              6
Shortest start_codon                         3
Shortest stop_codon                          3
Shortest three_prime_utr                     1
Shortest transcription_end_site              1
Shortest cds piece                           3
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                50
Shortest intron into exon part               50
Shortest intron into five_prime_utr part     50
Shortest intron into intron part             4
Shortest intron into three_prime_utr part    50

Re-compute transcript without isoforms asked. We remove shortest isoforms if any

Number of genes                              30603
Number of transcripts                        30603
Number of mrnas with utr both sides          30242
Number of mrnas with at least one utr        30601
Number of cdss                               30603
Number of exons                              231831
Number of five_prime_utrs                    30396
Number of introns                            179972
Number of start_codons                       29381
Number of stop_codons                        29514
Number of three_prime_utrs                   30447
Number of transcription_end_sites            29393
Number of exon in cds                        208588
Number of exon in five_prime_utr             46979
Number of exon in three_prime_utr            36768
Number of intron in cds                      177985
Number of intron in exon                     201228
Number of intron in five_prime_utr           16583
Number of intron in intron                   155513
Number of intron in three_prime_utr          6321
Number gene overlapping                      489
Number of single exon gene                   3976
Number of single exon transcript             3976
mean transcripts per gene                    1.0
mean cdss per transcript                     1.0
mean exons per transcript                    7.6
mean five_prime_utrs per transcript          1.0
mean introns per transcript                  5.9
mean start_codons per transcript             1.0
mean stop_codons per transcript              1.0
mean three_prime_utrs per transcript         1.0
mean transcription_end_sites per transcript  1.0
mean exons per cds                           6.8
mean exons per five_prime_utr                1.5
mean exons per three_prime_utr               1.2
mean introns in cdss per transcript          5.8
mean introns in exons per transcript         6.6
mean introns in five_prime_utrs per transcript0.5
mean introns in introns per transcript       5.1
mean introns in three_prime_utrs per transcript0.2
Total gene length                            348974335
Total transcript length                      346221435
Total cds length                             48222964
Total exon length                            80303427
Total five_prime_utr length                  9038176
Total intron length                          215643389
Total start_codon length                     88143
Total stop_codon length                      88542
Total three_prime_utr length                 23042287
Total transcription_end_site length          29393
Total intron length per cds                  211904194
Total intron length per exon                 265918008
Total intron length per five_prime_utr       39962084
Total intron length per intron               27047457
Total intron length per three_prime_utr      13007009
mean gene length                             11403
mean transcript length                       11313
mean cds length                              1575
mean exon length                             346
mean five_prime_utr length                   297
mean intron length                           1198
mean start_codon length                      3
mean stop_codon length                       3
mean three_prime_utr length                  756
mean transcription_end_site length           1
mean cds piece length                        231
mean five_prime_utr piece length             192
mean three_prime_utr piece length            626
mean intron in cds length                    1190
mean intron in exon length                   1321
mean intron in five_prime_utr length         2409
mean intron in intron length                 173
mean intron in three_prime_utr length        2057
Longest gene                                 232931
Longest transcript                           232931
Longest cds                                  63957
Longest exon                                 18417
Longest five_prime_utr                       16693
Longest intron                               87093
Longest start_codon                          3
Longest stop_codon                           3
Longest three_prime_utr                      18129
Longest transcription_end_site               1
Longest cds piece                            14759
Longest five_prime_utr piece                 16693
Longest three_prime_utr piece                18129
Longest intron into cds part                 87093
Longest intron into exon part                87093
Longest intron into five_prime_utr part      57732
Longest intron into intron part              14759
Longest intron into three_prime_utr part     42549
Shortest gene                                174
Shortest transcript                          174
Shortest cds                                 5
Shortest exon                                3
Shortest five_prime_utr                      1
Shortest intron                              6
Shortest start_codon                         3
Shortest stop_codon                          3
Shortest three_prime_utr                     1
Shortest transcription_end_site              1
Shortest cds piece                           3
Shortest five_prime_utr piece                1
Shortest three_prime_utr piece               1
Shortest intron into cds part                50
Shortest intron into exon part               50
Shortest intron into five_prime_utr part     50
Shortest intron into intron part             4
Shortest intron into three_prime_utr part    50

--------------------------------------------------------------------------------

Bye Bye.
```
