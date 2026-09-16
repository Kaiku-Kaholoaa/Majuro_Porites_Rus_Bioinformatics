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

## Fine Tuning Our Model for Chain Convergence

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

Great, now we can learn how to adjust our acceptance rates. Generally, they're created by these parameters:
```
       mSeedsProposalS2 = 0.010000
       qSeedsProposalS2 = 0.100000
       mEffctProposalS2 = 0.100000
       qEffctProposalS2 = 0.001000
      mrateMuProposalS2 = 0.010000
```
The easiest way to remember them is:

m = migration
q = local genetic diversity
Seeds = where the spatial regions are
Effct = what value each region has
ProposalS2 = how big an MCMC jump EEMS tries to make; S2 means variance
