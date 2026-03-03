# Functional Annotation of *Porites rus*

Functional Annotation Pipeline Adopted From Researcher Kevin Wong: 

https://kevinhwong1.github.io/KevinHWong_Notebook/Porites-astreoides-genome-annotation/ 

*Note: Please also be sure that you downloaded your databases!*

https://github.com/Kaiku-Kaholoaa/Majuro_Porites_Rus_Bioinformatics/blob/main/Downloading_Protein_Databases.md

## Align query protein sequences against databases

### 1. BLAST the protein sequences against Swiss-Prot (Highest Quality Protein Matches)

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

wc -l list_of_Pastgenemodelproteins_sprot.txt
#23,678
```

Then exclude these Gene Model names from your original fasta/.faa/protein file.
* needed to load the module that has the script with the -exclude command in it

#first, load the modules needed for exclude

```
ml system libpng 

ml biology ucsc-utils
```

Double check if the 'faSomeRecords' script is loaded as a command:

```faSomeRecords```

* faSomeRecords - Extract multiple fa records
* usage:
   * faSomeRecords in.fa listFile out.fa
* options:
   * -exclude - output sequences not in the list file.


Perfect! I then ran the -exclude command to exclude the blasted Gene Models from the .faa file

```faSomeRecords -exclude no_mpi_round3.2.all.maker.proteins.busco.fasta list_of_Pastgenemodelproteins_sprot.txt Past_proteins_names_v1.0.faa.prot4trembl ```


Checking to make sure it worked: 

```
grep "augustus_masked-CAXIVM010000010.1" list_of_Pastgenemodelproteins_sprot.txt | head
#augustus_masked-CAXIVM010000010.1-processed-gene-0.3-mRNA-1

grep ">augustus_masked-CAXIVM010000010.1" no_mpi_round3.2.all.maker.proteins.busco.fasta | head
#>augustus_masked-CAXIVM010000010.1-processed-gene-0.3-mRNA-1 protein AED:1.00 eAED:1.00 QI:0|-1|0|0|-1|1|1|0|1539
#>augustus_masked-CAXIVM010000010.1-processed-gene-0.2-mRNA-1 protein AED:1.00 eAED:1.00 QI:0|-1|0|0|-1|1|1|0|260
#>augustus_masked-CAXIVM010000010.1-processed-gene-0.5-mRNA-1 protein AED:1.00 eAED:1.00 QI:0|-1|0|0|-1|1|1|0|70

grep ">augustus_masked-CAXIVM010000010.1" list_of_Pastgenemodelproteins_sprot.txt | head
#file not present - good!
```

*Please also note that an AED of 1 means no transcriptional evidence for the gene model. This is ok because it could mean that you have the genetic evidence for it there, but the rna data could be missing from your transcriptome. As a result, this means that if your protein matches well to swissprot, then in reality you have a high confidence gene model since even though your transcripts are missing from your transcriptome, both gene and protein are present and are aligning to this high quality database.*

Checking the number of Gene Models:

`grep -c ">" no_mpi_round3.2.all.maker.proteins.busco.fasta 
#92,944 [total gene model count]

grep -c ">" Past_proteins_names_v1.0.faa.prot4trembl
#69,266 [unaligned gene models to swissprot]

wc -l list_of_Pastgenemodelproteins_sprot.txt
#23,678 [number of gene models aligned to swissprot]

#92,944 - 69,266 = 23,678
`
*Amazing fact check! Use this file to blast against trembl after we obtain the xml*

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

### 2. BLAST the remaining protein sequences against Trembl

#### First chunk the fasta file of protein models so that we can submit an array, and run all chuncks in parallel ###

```bash
mkdir trembl_chunks

awk -v n=500 '                   # Set variable n=500 (number of sequences per chunk)

/^>/ {                           # If the line starts with ">" (FASTA header line)

    if (++seq % n == 1) {        # Increment sequence counter.
                                  # Every time we hit a new header, seq increases by 1.
                                  # If seq modulo n equals 1, we start a new chunk file.
                                  # This happens for sequence 1, 501, 1001, etc.

        file=sprintf("trembl_chunks/chunk_%04d.faa", ++filecount)
                                  # Create new output filename.
                                  # %04d pads filecount with leading zeros (0001, 0002, ...)
                                  # filecount increments each time a new chunk starts.
    }
}

{ print >> file }                # Print the current line (header OR sequence line)
                                 # into the current chunk file.
                                 # >> means append to that file.

' Past_proteins_names_v1.0.faa.prot4trembl
```

Great, with your chunked files now we can blast them to tremble. 

`cat trembl_array_max5.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=M5array_trembl_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%30
#SBATCH --output=trembl_logs/%x_%A_%a.out
#SBATCH --error=trembl_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/trembl/trembl"
CHUNK_DIR="trembl_chunks"
OUT_DIR="trembl_out_max5"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Running chunk: $chunk_file"

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.out" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

For this dataset, we have 139 chunks (faa files) and thus 139 outputs in our Trembl_output dir. 
We can add them directly using cat, since there are no headers. 

`cat chunk_*out > PastGeneModels_vs_trembl_1e-5_max5.out`

Let's also do the same for the best hit proteins.

`cat trembl_array_max1.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=M1_array_trembl_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%40
#SBATCH --output=trembl_logs/%x_%A_%a.out
#SBATCH --error=trembl_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/trembl/trembl"
CHUNK_DIR="trembl_chunks"
OUT_DIR="trembl_out_max1"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Running chunk: $chunk_file"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.out" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

After navigating to the output dir, we can once again combine our chunks: 

`cat chunk_*out > PastGeneModels_vs_trembl_1e-5_max1.out`
`cat PastGeneModels_vs_trembl_1e-5_max1.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_trembl_1e-5_besthit.out`

We can then double check the amount of gene models that aligned to the database. 

```bash
wc -l PastGeneModels_vs_trembl_1e-5_besthit.out #42,918
```

Select the gene model proteins without hits in Trembl

```bash
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PastGeneModels_vs_trembl_1e-5_besthit.out > list_of_Pastgenemodelproteins_trembl.txt

* needed to load the module that has the script with the -exclude command in it

#first, load the modules needed for exclude

```
ml system libpng 

ml biology ucsc-utils
```

Double check if the 'faSomeRecords' script is loaded as a command:

```faSomeRecords```

Perfect! I then ran the -exclude command to exclude the blasted Gene Models from the .faa file

#then exclude these Gene Model names from your original fasta/.faa/protein file

faSomeRecords -exclude Past_proteins_names_v1.0.faa.prot4trembl list_of_Pastgenemodelproteins_trembl.txt Past_proteins_names_v1.0.faa.prot4nr

#check the number of Gene Models

grep -c ">" Past_proteins_names_v1.0.faa.prot4nr #26,348

#using this file to blast against nr database
```

Let's quickly double check:

```
grep -c ">" no_mpi_round3.2.all.maker.proteins.busco.fasta 
#92,944 [total gene model count]

wc -l list_of_Pastgenemodelproteins_sprot.txt
#23,678 [number of gene models aligned to swissprot]

grep -c ">" Past_proteins_names_v1.0.faa.prot4trembl
#69,266 [unaligned gene models to swissprot, for trembl]

cut -f1 PastGeneModels_vs_trembl_1e-5_max1.out |sort -u| wc -l
#42,918 [number of gene models aligned to trembl] 

#fact check: 69,266 - 42,918 = 26,348 unmapped protein models for NR

grep -c ">" Past_proteins_names_v1.0.faa.prot4nr #26,348 Perfect!
```

Awesome, let's create the xml file and continue with NR: 

`cat trembl_array_max5_xml.sbatch `

```bash
#!/bin/bash
#SBATCH --job-name=xml_M5array_trembl_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%30
#SBATCH --output=trembl_logs/%x_%A_%a.out
#SBATCH --error=trembl_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/trembl/trembl"
CHUNK_DIR="trembl_chunks"
OUT_DIR="trembl_out_max5_"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Running chunk: $chunk_file"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.xml" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

`cat chunk_*xml > PastGeneModels_vs_trembl_1e-5_max5.xml`


### 3. BLAST the remaining protein sequences against NR


#### First chunk the fasta file of protein models so that we can submit an array, and run all chuncks in parallel ###

```bash
mkdir nr_chunks

awk -v n=500 '                   # Set variable n=500 (number of sequences per chunk)

/^>/ {                           # If the line starts with ">" (FASTA header line)

    if (++seq % n == 1) {        # Increment sequence counter.
                                  # Every time we hit a new header, seq increases by 1.
                                  # If seq modulo n equals 1, we start a new chunk file.
                                  # This happens for sequence 1, 501, 1001, etc.

        file=sprintf("nr_chunks/chunk_%04d.faa", ++filecount)
                                  # Create new output filename.
                                  # %04d pads filecount with leading zeros (0001, 0002, ...)
                                  # filecount increments each time a new chunk starts.
    }
}

{ print >> file }                # Print the current line (header OR sequence line)
                                 # into the current chunk file.
                                 # >> means append to that file.

' Past_proteins_names_v1.0.faa.prot4nr
```

```bash
awk -v n=500 ' /^>/ {                          
    if (++seq % n == 1) { 
        file=sprintf("nr_chunks/chunk_%04d.faa", ++filecount)
    }
}
{ print >> file }
' Past_proteins_names_v1.0.faa.prot4nr
```

`cat nr_array_max1.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=M1_nr_array_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%40
#SBATCH --output=nr_max1_logs/%x_%A_%a.out
#SBATCH --error=nr_max1_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/nr/nr"
CHUNK_DIR="nr_chunks"
OUT_DIR="nr_out_max1"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Blast against swissprot database" $(date)

echo "Running chunk: $chunk_file"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.out" \
  -evalue 1e-5 \
  -max_target_seqs 1 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

And the same for max 5: 

`cat nr_array_max5.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=M5_nr_array_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%40
#SBATCH --output=nr_m5_logs/%x_%A_%a.out
#SBATCH --error=nr_m5_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/nr/nr"
CHUNK_DIR="nr_chunks"
OUT_DIR="nr_out_max5"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Blast against swissprot database" $(date)

echo "Running chunk: $chunk_file"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.out" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

Checking the outputs:

`cat chunk* > PastGeneModels_ncbi_max1.out`
`cat PastGeneModels_ncbi_max1.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PastGeneModels_vs_nr_1e-5_besthit.out`
`wc -l PastGeneModels_vs_nr_1e-5_besthit.out #59`

Lets create the for the xml file:

`cat nr_array_max5_xml.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=xml_M5_nr_array_blast
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH --array=1-139%40
#SBATCH --output=nr_xml_m5_logs/%x_%A_%a.out
#SBATCH --error=nr_xml_m5_logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu

set -euo pipefail

DB="/scratch/users/kaiku/databases/nr/nr"
CHUNK_DIR="nr_chunks"
OUT_DIR="nr_out_max5_xml"

chunk_file=$(ls ${CHUNK_DIR}/chunk_*.faa | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
base=$(basename "$chunk_file" .faa)

echo "Blast against swissprot database" $(date)

echo "Running chunk: $chunk_file"

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

blastp \
  -db "$DB" \
  -query "$chunk_file" \
  -out "${OUT_DIR}/${base}.out" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -outfmt "5 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

echo "STOP $(date)"
```

### 4. Interproscan

#### InterProScan identifies conserved protein domains and functional signatures in your proteins.

Think “Does this protein contain known functional domains"? This is what interproscan is used for. 

##### Step 1: Download InterProScan 5.59-91.0

```bash
cd /scratch/users/kaiku
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.59-91.0/interproscan-5.59-91.0-64-bit.tar.gz
tar -xzf interproscan-5.59-91.0-64-bit.tar.gz
cd interproscan-5.59-91.0
```
##### Step 2: Install Interproscan Databases

```bash
python3 setup.py -f interproscan.properties
```

##### Step 3: Chunk our large fasta output from Maker: 

```bash
cd /path/to/funcitonal_annotation_dir
mkdir ipr_chunks

awk '
/^>/ {
    if (seq % 5000 == 0) {
        file = sprintf("chunk_%03d.fa", ++chunk)
    }
    seq++
}
{ print > file }
' no_mpi_round3.2.all.maker.proteins.busco.fasta #<- maker output file
```

##### Step 4: Run the interproscan array 

`cat interproscan.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=ipr_chunk
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-19
#SBATCH --output=ipr_%A_%a.out
#SBATCH --error=ipr_%A_%a.err
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1

cd /scratch/users/kaiku/interproscan-5.59-91.0

# Unique TMP per array task
export TMPDIR=/scratch/users/kaiku/ipr_tmp_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Select chunk file
CHUNK=$(printf "chunk_%03d.fa" $SLURM_ARRAY_TASK_ID)

./interproscan.sh \
  -i /scratch/users/kaiku/feb26_pipeline/pipeline_dec_24/transcriptomics/dec25_maker2_mpi/functional_annotation/ipr_chunks/$CHUNK \
  -f XML \
  -b Past.interpro.chunk_${SLURM_ARRAY_TASK_ID} \
  -iprlookup \
  -goterms \
  -pa \
  -cpu 12 \
  -exclappl ProSitePatterns,ProSiteProfiles \
  -T $TMPDIR

echo "DONE $(date)"
```

And then we'll merged the chunk'ed xml output files using this code:

```bash 
cd /scratch/users/kaiku/interproscan-5.59-91.0

head -n 1 Past.interpro.chunk_1.xml > Past.interpro.all.xml

for f in Past.interpro.chunk_*.xml; do
    grep -v '<?xml' "$f" | \
    grep -v '<protein-matches' | \
    grep -v '</protein-matches>' >> Past.interpro.all.xml
done

echo "</protein-matches>" >> Past.interpro.all.xml
```

Nice! Also convert to gff and tsv so that we can run statistics with agat:

`cat interproscan_xml2gff.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=ipr_xml2gff
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=ipr_xml2gff.out
#SBATCH --error=ipr_xml2gff.err
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1

cd /scratch/users/kaiku/interproscan-5.59-91.0

# Unique TMP per array task
export TMPDIR=/scratch/users/kaiku/ipr_tmp_xml
mkdir -p $TMPDIR

export _JAVA_OPTIONS="-Xms8G -Xmx60G"

echo "Converting XML → GFF3"
./interproscan.sh \
   -mode convert \
   -f GFF3 \
   -i Past.interpro.all.xml \
   -o Past.interpro.all.gff3

echo "Converting XML → TSV"
./interproscan.sh \
   -mode convert \
   -f TSV \
   -i Past.interpro.all.xml \
   -o Past.interpro.all.tsv

cp Past.interpro.all.* /scratch/users/kaiku/feb26_pipeline/pipeline_dec_24/transcriptomics/dec25_maker2_mpi/functional_annotation/

echo "DONE $(date)"
```

### 5. Statistics with AGAT

#### Porites rus

`cat agat_stats.sbatch`

```bash
#!/bin/bash
#SBATCH --job-name=agat
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=stats_agat.out
#SBATCH --error=stats_agat.err
#SBATCH -p serc,spalumbi,hns
#SBATCH --nodes=1

cd /scratch/users/kaiku/interproscan-5.59-91.0

source /home/users/kaiku/miniconda3/etc/profile.d/conda.sh

conda activate agat_env

agat_sp_statistics.pl --gff Past.interpro.all.gff3 

echo "DONE $(date)"
```

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

Bye Bye.
```
