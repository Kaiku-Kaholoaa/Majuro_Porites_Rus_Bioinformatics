# EEMS
This section is for reproducibility running EEMS.. we've had some issues!

## Installation
Installing EEMS was somewhat difficult. Although we were able to create a conda environment and get dependencies for EEMS, we found issues with the makefile compiler, which only worked on previous versions of both eigen and boost. This was later resolved thanks to this: https://github.com/dipetkov/eems/issues/14 . Basically, we needed to reinstall eigen and boost into our eems conda environment, and change the makefile paths to reflect these changes: 

```
$ cat /scratch/users/kaiku/eems/runeems_snps/src/Makefile

EIGEN_INC = /scratch/users/kaiku/conda_envs/eems_env/include/eigen3
BOOST_LIB = /scratch/users/kaiku/conda_envs/eems_env/lib
BOOST_INC = /scratch/users/kaiku/conda_envs/eems_env/include
```

With this, we were able to get EEMs working, but we then struggled with our input files, which were not created correctly. 

## Data Cleaning of Input Files
For EEMS, we need three input files and a paramters .ini file.

### datapath.fam file
The fam file in my opinion is one of the most important for data cleaning. This is because order of samples must be preserved across all input files, and it all starts with the fam file. This is essentially the foundation from which we build the rest of our files. For our analysis, we took the rare allele dataset (no_maf_no_singletons), converted them to bed files (for bed2diffs), and renamed them to prus_eems. This dataset contains 159 individuals, and so we must preserve this number and their order across all other input files. 

### datapath.diffs file
Next is the diffs file that we created using bed2diffs. After installation, we did something similar to the documentation provided by using 

`
./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode --nthreads 2
`
source: https://github.com/dipetkov/eems/tree/master/bed2diffs

### datapath.coord file
This was the problem maker for my analysis. This is because we had ~179 samples before qc, and this file had sample IDS followed by LAT LON. The first issue is that we included samples that did not pass qc, and the second issue is that sample IDS should not be in this file. It should only have LAT LONs, which emphasizes the need for correct ordering. However, both issues were resolved using an AWK command, and I'm proud to say i'm learning AWK pretty well! 

Using awk, I was able to use the fam file to make a list of IDs to keep, and then print only the lat lons (from the incorrect file) if their IDs were stored from the fam file!

``` bash
awk 'FNR==NR{keep[$2]=1; next} ($1 in keep){print $2,$3}' prus_eems.fam prus_eems_incorrect_coords.txt
```

```
# FNR==NR 
  # if total rows read == rows read in current file (essentially means "for the first file only {do action}"
# {keep[$2]=1; next} 
  # create list called 'keep', and store the second column of the fam file (IDS) with variable of 1 (not used, but needed to make the list). Then move to next row via next. 
# ($1 in keep)
  # If the ID ($1 in the incorrect_coords file) matches a key in the list called keep (from fam file),
# {print $2,$3}
  # print the lat ($2) and lon ($3) from the second file. 
```

`head prus_eems.coord`
```
7.170308 171.13354
7.170308 171.13354
7.111973 171.120796
7.111973 171.120796
```
yay! 

Then with this we can check concurrency of order using the fam file, .order file (from bed2diffs), and in our newly created .coord file. 
and trust me, they looked good! 

### datapath.ini file 
