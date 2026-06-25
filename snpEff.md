# Learning snpEff to annotate our vcf of interesting SNPs

Earlier we conducted a GWAS via gemma to identify peaks associated with the top5% heat experienced by corals. We also used pcadapt to
evaluate outlier snps associated with our first two principle components from our PCA, which showed a slight grouping with a thermal outgroup. 
These snps were then converted into two different vcf files by filtering the major vcf by the positions of SNPs of interest. 

We also learned that minor allele frequencies can identify hidden trends in human population structure, and can be extremely useful for
our analysis given our fine spatial scale (< 40km) and our highly-related, mostly panmictic population. The minor allele frequency pca 
identified that there is a clear subgroup among our thermal individuals, to which we grouped together based on a PC1 threshold. 
The alleles of this subgroup (n_indiv = 19; total_possible_alleles = 38), were then filtered for private alleles at high frequency. This
represents alelles that are not present in the greater population, display a frequency of > 20%, and have a total mimimum allele count of 8.
In the end we idenitfied ~800 alleles at high frequency not present in the main population. These high frequency private alleles were 
converted into a vcf, and will be supplied to snpEff to idenitfy snp locations (intergenic, intron, upstream/downstream variants, etc.), 
and mutational changes such as synonymous, missense, stop_gained, splice_acceptor/donor variants, etc. 

snpEff is able to annotate these vcf files by mapping positions to a database. The database must first be created using the gff that we made of 
the Porites rus genome in our Gene_Annotation_Pipeline.md.

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
Great, got the SwissProt Database! Now lets do the same thing for TrEMBL (Translated EMBL Nucleotide Sequence Data Library) 

## Download required databases (TrEMBL)

`cat build_trembl_db.sbatch`

```
#!/bin/bash
#SBATCH --job-name=build_trembl_db
#SBATCH -p hns,spalumbi,serc
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu
#SBATCH --output=build_trembl_db.out
#SBATCH --error=build_trembl_db.err

set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# ============================
# CONFIGURATION
# ============================

DB_DIR=/scratch/users/kaiku/databases/trembl
UNIPROT_URL=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
FASTA_GZ=uniprot_trembl.fasta.gz
FASTA=uniprot_trembl.fasta
DB_NAME=trembl

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
# DOWNLOAD TrEMBL
# ============================

if [ ! -f ${FASTA_GZ} ]; then
    echo "Downloading TrEMBL FASTA from UniProt..."
    wget ${UNIPROT_URL}
else
    echo "TrEMBL FASTA archive already exists, skipping download."
fi

# ============================
# UNZIP FASTA
# ============================

if [ ! -f ${FASTA} ]; then
    echo "Unzipping TrEMBL FASTA..."
    gunzip -c ${FASTA_GZ} > ${FASTA}
else
    echo "Uncompressed FASTA already exists, skipping unzip."
fi

# ============================
# BUILD BLAST DATABASE
# ============================

echo "Building TrEMBL BLAST database..."
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

echo "TrEMBL BLAST database successfully built."
echo "Job finished at $(date)"
```

Awesome! Now lets get the nr database:

## Download required databases (NR)

`cat build_nr_db.sbatch`

```
#!/bin/bash
#SBATCH --job-name=build_nr_db
#SBATCH -p hns,spalumbi,serc
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kaiku@stanford.edu
#SBATCH --output=build_nr_db.out
#SBATCH --error=build_nr_db.err

set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# ============================
# CONFIGURATION
# ============================

DB_DIR=/scratch/users/kaiku/databases/nr
NCBI_FTP=https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
FASTA_GZ=nr.gz
FASTA=nr
DB_NAME=nr

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
# DOWNLOAD NR FASTA
# ============================

echo "Downloading NR FASTA from NCBI..."
#wget ${NCBI_FTP}/${FASTA_GZ}

# ============================
# UNZIP NR
# ============================

if [ ! -f ${FASTA} ]; then
    echo "Unzipping NR (this will take time)..."
    gunzip -c ${FASTA_GZ} > ${FASTA}
else
    echo "Uncompressed NR already exists, skipping unzip."
fi

# ============================
# BUILD BLAST DATABASE
# ============================

echo "Building NR BLAST database (this may take many hours)..."

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

echo "NR BLAST database successfully built."
echo "Job finished at $(date)"
```

Great work! :)
