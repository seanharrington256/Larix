## Larix demographics

### Collaboration with Dan Turck & Dave Tank

<br>

This repository contains my code for fitting demographic models and estimating model parameters for Dan's Larix data.


### Create SFS

FastSimCoal takes the site frequency spectrum (SFS) as input. To generate the SFS, I used Isaac Overcast's `easySFS.py` script from here: [https://github.com/isaacovercast/easySFS](https://github.com/isaacovercast/easySFS ). To determine how many alleles to project down (downsample to) to account for missing data, I ran the first line to preview the sizes of datasets at various projections. The second lines generates the SFS is the projection I used:


```
	easySFS.py -i 4LAR_no_outgroups.recode_dp07mis07.recode.vcf -p LAR_Pop_file_without_outgroups_3_populations.txt --preview

	easySFS.py -i 4LAR_no_outgroups.recode_dp07mis07.recode.vcf -p LAR_Pop_file_without_outgroups_3_populations.txt --proj=8,14,22 --prefix LarK3 -o /Users/harrington/Active_Research/Larix_data_and_outs/SFS
```

The file output from this is `LarK3_MSFS.obs`.

### Divergence time

Divergence time at the root is fixed to either 20k years or 200k years. Assuming a generation time of 180 years, these times become 111 or 1,111 generations. We fit all models under each of these assumed divergence times. 



### Models


1. Make the models. Making the models and ensuring that they are specified correctly is the most tedious part of this. Models are contained in the `Models` directory. Each contains a `.tpl` and `.est` file that specify the coalescent model for FastSimCoal. We use fastsimcoal26 v2.6.0.3.  

To test that each model is correctly specified, I execute a very short test run and visualize the output `.par` file using the `ParFileInterpreter-v6.3.1.r` script. This was previously available on the site hosting FSC2, but has since been replaced by a new version. It is included here in this repo.


* Example of checking a model with a short run

A very short FSC run:
`fsc26 -t *tpl -e *est -n 1000 -L 40 -m -M -u -x  -q`

Check that it looks right - don't worry about parameter estimates, it didn't run long enough to optimize (note this uses a path to the script that exists on my computer):
`Rscript /Applications/fsc26_mac64/ParFileInterpreter-v6.3.1.r *par`


2. Get the linux executable for FSC: http://cmpg.unibe.ch/software/fastsimcoal2/ (for running this on a Linux cluster, as I did).

3. Put the `Models` directory (without any name changes or scripts below won't work) onto the cluster. The `Models` directory must contain a directory for each model to be estimated (each containing a `.tpl` and `.est` file) - the current `Models` directory does this - as well as the sfs file. 

## EDIT THIS: You will need to add the sfs `.obs` file from Dryad if running this.

4. To run multiple replicates of FSC parameter estimation, use the `Prep_FSC_reps.sh` script. Details of how to use this script are in the comments inside it. Briefly, run it from inside the Models folder made in step 6, supplying an argument specifying if you are using either multi-dimensional SFS (argument MSFS) or two-dimensional SFS (argument jointMAF). This will make 50 replicates in a directory "Reps" in the parent directory of Models, each containing the necessary tpl, est, and obs file.

5. Submit all of the jobs to the cluster as a big job array on a SLURM cluster using `FSC_Larix_Mods.slurm`.


6. For each of the models find the single run with the best likelihood. The script `Get_best_FSCacross_mods.sh` wraps the `fsc-selectbestrun.sh` script from here: [https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh](https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh) (and copied here with some slight modifications) across each model. Comments in this script further describe its usage--run it from within the `Reps` directory. To get AIC from the bestlhoods files, go inside `best_L_allMods` directory, created by last script, and run the script `Get_AIC_across_mods.R`(!!!NOTE: This AIC calculation will only work if you output 1 and exactly 1 parameter to this file for each estimated parameter otherwise the AIC values will be incorrect -- i.e., if you have a complex parameter, and you output both the estimated parameter and the complex parameter that is a transformation of the estimated parameter, you'll need a different solution here).

7. Convert the parameter estimates into useful units (e.g., years, number of migrants per generation, etc.) -- I still don't have a fully automated solution for this, so has to be done on an ad hoc basis depending on what parameters are included, etc. `Par_conv_FSC_Larix.R` will do these conversions here - see comments in the script for more info. This script will also start some prep for parametric bootstrap estimation of confidence intervals around parameter estimates.

* NOTE: I think there are other times that aren't being converted - times for when ancient migration ends, possibly others



Hunting down failed runs:
`grep CANCELLED *Larix_7145250*  | cut -d "_" -f 5 | cut -d "." -f 1 | sort > failed2_fsc.txt`

`grep 'Bad parameters'  *.err` This was only one, and also contains "CANCELLED", don't need to get this separately.



- edit sample sizes
- check that all params exist match up in the est file

# just notes to me:

models made (Xs have been checked):
numbers correspond to numbers in Dan's diagram

20k div time:
1. 20_noMK3 - X
2. 20_AncMK3 - X
3. 20_MallAsymK3 - X
4. 20_AncMintK3 - X
5. 20_SecK3 - X
6. 20_SecIntK3 - X
7. 20_AncMCoInK3 - X

200k div time:
1. 200_noMK3 - X
2. 200_AncMK3 - X
3. 200_MallAsymK3 - X
4. 200_AncMintK3 - X
5. 200_SecK3 - X
6. 200_SecIntK3 - X
7. 200_AncMCoInK3 - X
16. 200_secLGMintK3 - X

200k div time with population resize (expecting expansion)added in:
200_AncMCoInExpK3 - 200_AncMCoInK3 with population resize parameter added in - X
200_MallAsymExpK3 - 200_MallAsymK3 with population resize parameter added in - X
200_noMExpK3 - 200_noMK3 with population resize parameter added in - X
200_SecExpK3 - 200_SecK3 with population resize parameter added in - X







- use `4LAR_no_outgroups.recode_dp07mis07.recode.vcf`


