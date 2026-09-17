# Running EEMS! 
Yes! We've been waiting for this! 

## Initial Trial Run 
To run EEMs, start with the path of the source file and then supply your .ini file to the params argument. Then select a seed for testing, and finally call your input .coord, .outer, and .diffs files. 

`
/scratch/users/kaiku/eems/runeems_snps/src/runeems_snps --params prus_eems_run1.ini --seed 938271 prus_eems.coord prus_eems.outer prus_eems.diffs
`

Awesome, so EEMS is now running with default settings for your priors and hyperparameters. When it finished, you can enter your output dir and cat the eemsrun.txt file. 
`cat prus_eems_output_run1/eemsrun.txt`
```
Input parameter values:
               datapath = ./prus_eems
               mcmcpath = ./prus_eems_output_run1
               prevpath = 
               gridpath = 
               distance = euclidean
                diploid = 1
                 nIndiv = 159
                 nSites = 22953951
                 nDemes = 200
                   seed = 938271
            numMCMCIter = 2000000
            numBurnIter = 1000000
            numThinIter = 9999
              negBiSize = 10
              negBiProb = 0.670000
             qVoronoiPr = 0.250000
             mrateShape = 0.000500
             qrateShape = 0.002000
             sigmaShape = 0.001000
             qrateScale = 0.500000
             mrateScale = 2.000000
             sigmaScale = 1.000000
       mSeedsProposalS2 = 0.010000
       qSeedsProposalS2 = 0.100000
       mEffctProposalS2 = 0.100000
       qEffctProposalS2 = 0.001000
      mrateMuProposalS2 = 0.010000

Acceptance proportions:
	(31986/125086) = 26% for proposal type "qTileRate",		 with proposal variance "qEffctProposalS2"
	(6005/124681) = 4.8% for proposal type "qTileMove",		 with proposal variance "qSeedsProposalS2"
	(26033/125852) = 21% for proposal type "qBirthDeath"
	(286062/374560) = 76% for proposal type "mTileRate",		 with proposal variance "mEffctProposalS2"
	(106772/249985) = 43% for proposal type "mMeanRate",		 with proposal variance "mrateMuProposalS2"
	(98243/374763) = 26% for proposal type "mTileMove",		 with proposal variance "mSeedsProposalS2"
	(100838/375071) = 27% for proposal type "mBirthDeath"
	(36728/250002) = 15% for proposal type "degrees of freedom"

Final log prior: 25.63
Final log llike: 38304.25
```

Great, the next step is fine tuning the model so that it works properly :)

## Fine Tuning Our Model

### Understanding Acceptance Proportions
Okay, the EEMS output seems like a lot, but for now we can just focus on the acceptance proportions, which we want to be between 20% and 40%. 
These are essentially how the chain explores probable space and the likelyhood it accepts or rejects a stepping-stone movement based on genotype similarity.
Here we are aiming for a **balance between step size and acceptance between proposals,** while exploring a decent amount of probably space.

If the acceptance rate is near 0, then this reflects a stepping size that is too large and all proposals are likely to get rejected. Ex:
```
0.001 -> 79.1 REJECT 
0.001 -> 34.4 REJECT 
0.001 -> 54.3 REJECT  
```
Result: chain did not move because steps were too large. Probable space was not explored efficiently. 

On the other hand, if the acceptance is near 1 then this reflects a stepping size that is too small and all proposals are likely to get accepted. Ex:
```
0.001 -> 0.002 ACCEPT 
0.002 -> 0.003 ACCEPT 
0.003 -> 0.002 ACCEPT  
```
Result: chain moved like 0.001 between steps, and didn't explore a good range of probable space. Bad.

A good acceptance rate like 20%-40% however, allows us to explore space evenly, which is what we want :) Ex:
```
0.001 -> 10.5 ACCEPT 
10.5 -> 35.2 REJECT 
35.2 -> 22.9 ACCEPT  
```
Notice how the chain is able to move through space efficiently, given these more moderate steps! 

### Adjusting our parameters to fine tune our acceptance proportions: 
Great, now we can learn how to adjust our acceptance rates. Generally, they're created by these parameters:
```
       mSeedsProposalS2 = 0.010000 #where the migration regions are located, affects mTileMove
       qSeedsProposalS2 = 0.100000 #where the diversity regions are located, affects qTileMove
       mEffctProposalS2 = 0.100000 #how high/low migration is in one migration region, affects mTileRate
       qEffctProposalS2 = 0.001000 #how high/low diversity is in one diversity region, affects qTileRate
      mrateMuProposalS2 = 0.010000 #the overall baseline migration rate across the whole map
```
The easiest way to remember them is:

m = migration
q = local genetic diversity
Seeds = where the spatial regions are
Effct = what value each region has
ProposalS2 = how big an MCMC jump EEMS tries to make; S2 means variance

So, given our previous acceptance results, let's adjust them to get the acceptance proportions around 20-40%:

```bash
Acceptance proportions:
	(31986/125086) = 26% for proposal type "qTileRate",			#qEffctProposalS2 IS GOOD
	(6005/124681) = 4.8% for proposal type "qTileMove",		 	#qSeedsProposalS2 VERY LOW
	(26033/125852) = 21% for proposal type "qBirthDeath"		#generally ignore, parameters dont directly influence this
	(286062/374560) = 76% for proposal type "mTileRate",		#mEffctProposalS2 VERY HIGH
	(106772/249985) = 43% for proposal type "mMeanRate",		#mrateMuProposalS2 SLIGHTLY HIGH
	(98243/374763) = 26% for proposal type "mTileMove",			#mSeedsProposalS2 GOOD
	(100838/375071) = 27% for proposal type "mBirthDeath"		#generally ignore, parameters dont directly influence this
	(36728/250002) = 15% for proposal type "degrees of freedom"	#generally ignore, parameters dont directly influence this
```
So, lets adjust these parameters from round 1 by making a new round2.ini file with a new output path

`nano prus_eems_run2.ini`
```bash
datapath = ./prus_eems
mcmcpath = ./prus_eems_output_run2

nIndiv = 159
nSites = 22953951
nDemes = 200

diploid = true

numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999

mSeedsProposalS2 = 0.01
qSeedsProposalS2 = 0.01

mEffctProposalS2 = 6.0
qEffctProposalS2 = 0.001

mrateMuProposalS2 = 0.05


#previous values:			#Adjustment made in this run:
#mSeedsProposalS2 = 0.01	#KEPT SAME
#qSeedsProposalS2 = 0.10	#SLIGHTLY INCREASED

#mEffctProposalS2 = 0.10	#INCREASED A LOT
#qEffctProposalS2 = 0.001	#KEPT SAME

#mrateMuProposalS2 = 0.01	#SLIGHTLY INCREASED
```

I did this a few times, and it took about 4 iterations to get what I wanted.. and it looks good!
```bash
       mSeedsProposalS2 = 0.009000
       qSeedsProposalS2 = 0.007000
       mEffctProposalS2 = 4.500000
       qEffctProposalS2 = 0.000700
      mrateMuProposalS2 = 0.035000

Acceptance proportions:
	(35354/125265) = 28% for proposal type "qTileRate",		 with proposal variance "qEffctProposalS2"
	(33537/124585) = 27% for proposal type "qTileMove",		 with proposal variance "qSeedsProposalS2"
	(26695/125450) = 21% for proposal type "qBirthDeath"
	(113325/374568) = 30% for proposal type "mTileRate",		 with proposal variance "mEffctProposalS2"
	(66154/250273) = 26% for proposal type "mMeanRate",		 with proposal variance "mrateMuProposalS2"
	(102338/374534) = 27% for proposal type "mTileMove",		 with proposal variance "mSeedsProposalS2"
	(182368/375283) = 49% for proposal type "mBirthDeath"
	(36404/250042) = 15% for proposal type "degrees of freedom"
```
yay! now we can evaluate chains using R:

## Evaluating Chain Convergence
For our chains, we want to see the plot generally centered on the y-axis with hills and valleys. This demonstrates that the model is binding itself to the likely posterior. A bad plot would be a diagonal line that is increasing like a 1:1 line. We don't want that because it shows that the model continually increases and is not bound well by priors or parameters. 

To do this, we can load R on the server, and then make an rscript. Our goal is to read in the mcmcpilogl file which contains our log priors and log likelihoods, and then take those values to mathematically build the posterior (Posterior = prior + likelihood on the log scale). With this, we can then plot it and save it as a PDF (and it's called a traceplot!):

```nano check_chain.R```
```
# Read EEMS posterior output
chain <- read.table("prus_eems_output_run4/mcmcpilogl.txt")

# EEMS writes two columns:
# column 1 = log prior
# column 2 = log likelihood
colnames(chain) <- c("log_prior", "log_likelihood")

# Posterior = prior + likelihood on the log scale
chain$log_posterior <- chain$log_prior + chain$log_likelihood

# Look at the data
print(head(chain))
print(summary(chain))

# Save plots to a PDF
pdf("prus_run4_chain.pdf", width = 8, height = 8)

par(mfrow = c(3, 1))

plot(
    chain$log_posterior,
    type = "l",
    xlab = "Saved MCMC sample",
    ylab = "Log posterior",
    main = "Posterior trace"
)

plot(
    chain$log_likelihood,
    type = "l",
    xlab = "Saved MCMC sample",
    ylab = "Log likelihood",
    main = "Likelihood trace"
)

plot(
    chain$log_prior,
    type = "l",
    xlab = "Saved MCMC sample",
    ylab = "Log prior",
    main = "Prior trace"
)

dev.off()
```
[prus_run4_chain.pdf](https://github.com/user-attachments/files/32315102/prus_run4_chain.pdf)

Nice! The posterior traceplot looks pretty good! The only concern however is that our sample size is small at ~100. We can increase this by reducing thinning from 9999 to 999. 

## Creating Final Outputs
When running eems, we essentially need to run the model 4 times, on 4 random seeds. These are therefore 4 independent chains, and we want to plot the chains together to ensure that they are all generally exploring the same posterior space. If two chains are clumping together while others are not, this could be a potential issue. the largest red flag would be if some chains do things that no other chains do. That is bad and demonstrates that the chains can randomly move around due to chance, and not due to the data. 

The first step is to create 4 new .ini files that mimic the thinning changes and new output directories. Here are the first 15 rows of 2 of my chains, and don't forget to make the other two! And remember to keep the parameters the same!:

``` head -15  prus_eems_final_chain*.ini```
```
==> prus_eems_final_chain1.ini <==
datapath = ./prus_eems
mcmcpath = ./prus_eems_output_final_chain1

nIndiv = 159
nSites = 22953951
nDemes = 200

diploid = true

numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 999

mSeedsProposalS2 = 0.009
qSeedsProposalS2 = 0.007

==> prus_eems_final_chain2.ini <==
datapath = ./prus_eems
mcmcpath = ./prus_eems_output_final_chain2

nIndiv = 159
nSites = 22953951
nDemes = 200

diploid = true

numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 999

mSeedsProposalS2 = 0.009
qSeedsProposalS2 = 0.007
```

and finally, here is the sbatch script we can use to run all 4 chains with their different random seeds :) :



