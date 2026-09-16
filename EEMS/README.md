# EEMS
This section is for reproducibility running EEMS.. we've had some issues!

## Installation
Installing EEMS was somewhat difficult. Although we were able to create a conda environment and get dependencies for EEMS, we found issues with the makefile compiler, which only worked on previous versions of both eigen and boost. This was later resolved thanks to this: https://github.com/dipetkov/eems/issues/14 . Basically, we needed to reinstall eigen and boost into our eems conda environment, and change the makefile paths to reflect these changes: 

```
$ cat /scratch/users/kaiku/eems/runeems_snps/src/Makefile

EIGEN_INC = /scratch/users/kaiku/conda_envs/eems_env/include/eigen3
BOOST_LIB = /scratch/users/kaiku/conda_envs/eems_env/lib
BOOST_INC = /scratch/users/kaiku/conda_envs/eems_env/include

...
```
With this, we were able to get EEMs working, but we then struggled with our input files, which were not created correctly. 

## Data Cleaning of Input Files
For EEMS, we need three input files and a paramters .ini file.

### .fam file
The fam file in my opinion is one of the most important for data cleaning. This is because order of samples must be preserved across all input files, and it all starts with the fam file. This is essentially the foundation from which we build the rest of our files. 

### .coord file

### .diffs file

### .ini file 
