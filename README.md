## Larix demographics

### Collaboration with Dan Turck & Dave Tank

<br>

This repository contains my code for fitting demographic models and estimating model parameters for Dan's Larix data.


### Models

Models are contained in the `Models` directory. Each contains a `.tpl` and `.est` file that specify the coalescent model for FastSimCoal. We use fastsimcoal26. 

To test that each model is correctly specified, I execute a very short test run and visualize the output `.par` file using the `ParFileInterpreter-v6.3.1.r` script. This was previously available on the site hosting FSC2, but has since been replaced by a new version. We include it here in this repo.




- add in complete FSC version
- edit sample sizes
- check that all params exist match up in the est file

# just notes to me:

models made (still need to be checked):

1. noMK3
2. AncMK3
3. MallAsymK3
4. AncMintK3
5. SecK3
6. SecIntK3
7. AncMCoInK3


- Need to check that the number of migration matrices are correct in all