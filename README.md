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

Models are contained in the `Models` directory. Each contains a `.tpl` and `.est` file that specify the coalescent model for FastSimCoal. We use fastsimcoal26. 

To test that each model is correctly specified, I execute a very short test run and visualize the output `.par` file using the `ParFileInterpreter-v6.3.1.r` script. This was previously available on the site hosting FSC2, but has since been replaced by a new version. It is included here in this repo.


* Example of checking a model with a short run

A very short FSC run:
`fsc26 -t *tpl -e *est -n 1000 -L 40 -m -M -u -x  -q`

Check that it looks right - don't worry about parameter estimates, it didn't run long enough to optimize (note this uses a path to the script that exists on my computer):
`Rscript /Applications/fsc26_mac64/ParFileInterpreter-v6.3.1.r *par`



- FSC 2.6.0.3
- edit sample sizes
- check that all params exist match up in the est file

# just notes to me:

models made (Xs have been checked):

20k div time:
1. 20_noMK3 - X
2. 20_AncMK3 - X
3. 20_MallAsymK3 - X
4. 20_AncMintK3 - X
5. 20_SecK3 - X
6. 20_SecIntK3 - X
7. 20_AncMCoInK3 - 

200k div time:
1. 200_noMK3 -
2. 200_AncMK3 - 
3. 200_MallAsymK3 - 
4. 200_AncMintK3
5. 200_SecK3
6. 200_SecIntK3
7. 200_AncMCoInK3


Add in rules for events being above 0



- add in model 16

- fit these models, then add in models that incorporate population size change


- use `4LAR_no_outgroups.recode_dp07mis07.recode.vcf`


