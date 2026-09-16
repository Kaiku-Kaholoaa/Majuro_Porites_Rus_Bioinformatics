# Installation and Data Preparation and for EEMS
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

### datapath.fam and datapath.diffs files
The fam file in my opinion is one of the most important for data cleaning. This is because order of samples must be preserved across all input files, and it all starts with the fam file. This is essentially the foundation from which we build the rest of our files. For our analysis, we took the rare allele dataset (no_maf_no_singletons), converted them to bed files (for bed2diffs), and renamed them to prus_eems. This dataset contains 159 individuals, and so we must preserve this number and their order across all other input files. 

File name: prus_eems.fam

Next is the **datapath.diffs** file that we created using bed2diffs. After program installation via https://github.com/dipetkov/eems/tree/master/bed2diffs , we followed the documentation to produce our prus_eems.diffs file. Here they also recommended to remove snps with high missingess (which we did via qa/qc), and proceeding with bed2diffs_v1 instead of bed2diffs_v2:
`
./src/bed2diffs_v1 --bfile ./test/example-SNP-major-mode --nthreads 2
`
The result is our unlabeled, pairwise genotype matrix that can be used by eems. 

File name: prus_eems.diffs

### datapath.coord file
This was the problem maker for my analysis. This is because we had ~179 samples before qc, and this file had sample IDS followed by LON LAT. The first issue is that we included samples that did not pass qc, and the second issue is that sample IDS should not be in this file. It should only have LON LATs, which emphasizes the need for correct ordering. However, both issues were resolved using an AWK command, and I'm proud to say i'm learning AWK pretty well! 

Using awk, I was able to use the fam file to make a list of IDs to keep, and then print only the lat lons (from the incorrect file) if their IDs were stored from the fam file!

``` bash
awk 'FNR==NR{keep[$2]=1; next} ($1 in keep){print $3,$2}' prus_eems.fam prus_eems_incorrect_coords.txt > prus_eems.coord

#also doing one with the IDS kept so we can do some bookkeeping later:
awk 'FNR==NR{keep[$2]=1; next} ($1 in keep){print}' prus_eems.fam prus_eems_incorrect_coords.txt > prus_eems_coords_with_IDs.txt

```

```
# FNR==NR 
  # if total rows read == rows read in current file (essentially means "for the first file only {do action}"
# {keep[$2]=1; next} 
  # create list called 'keep', and store the second column of the fam file (IDS) with variable of 1 (not used, but needed to make the list). Then move to next row via next. 
# ($1 in keep)
  # If the ID ($1 in the incorrect_coords file) matches a key in the list called keep (from fam file),
# {print $3,$2}
  # print the lon ($3) and lat ($2) from the second file. 
```

`head prus_eems.coord` 
```
171.13354 7.170308
171.13354 7.170308
171.120796 7.111973
171.120796 7.111973
171.120796 7.111973
```
Yay! Also remember it needs to be in LON LAT formatting!  

Then with this we can check concurrency of order using the fam file, .order file (from bed2diffs), and in our newly created .coord file using another awk command:
```bash
awk '{print $2, NR}' prus_eems.fam > fam_order.txt
awk '{print $2, NR}' prus_eems.order > diffs_order.txt
awk '{print $1, NR}' prus_eems_coords_with_IDs.txt > sample_coords_order.txt
```
With these, we can just ensure sample size (wc -l) and order (tail) are consistent across files!

```bash
wc -l fam_order.txt diffs_order.txt sample_coords_order.txt

fam_order.txt        159
diffs_order.txt      159
sample_coords_order  159

tail fam_order.txt diffs_order.txt sample_coords_order.txt

fam_order.txt
P-rus_395_S220 156
P-rus_396_S161 157
P-rus_397_S162 158
P-rus_399_S222 159

diffs_order.txt
P-rus_395_S220 156
P-rus_396_S161 157
P-rus_397_S162 158
P-rus_399_S222 159

sample_coords_order.txt
P-rus_395_S220 156
P-rus_396_S161 157
P-rus_397_S162 158
P-rus_399_S222 159
```

Awesome! Sample input files ready, but now we need our .outer file outlining our spatial polygon (the area around and between our samples).

### datapath.outer file
This file is essentially our spatial bounds, and should encompass the locations of our collected samples (prus_eems.coord). The way I did this was I literally plotted points on google maps, extracted the lon lats, and ordered it the way the EEMS needs it, which is counter clockwise and **closed** which means the first point in the polygon is also the last. 

Here is a snippet of that file: 
```bash 
head -3 *.outer
171.057210868626811	7.226999696692798
171.039859607264788	7.189203746942331
171.021598798360600	7.150927591577632

tail -3 *.outer
171.090024025820099	7.208316356772613
171.087898211318702	7.221695887737019
171.057210868626811	7.226999696692798
```
Remember that both the .outer and .coord files need to be formatted as LON LAT, but otherwise, that's it for our input files! woot woot!

### datapath.ini file 
