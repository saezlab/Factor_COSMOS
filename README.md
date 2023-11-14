# NCI60_mofa_cosmos

Tutorial on how to use the MOFA-COSMOS pipeline in order to generate testable hypothesis connecting MOFA outputs with COSMOS network inference.

## TLDR:
The main pipleine script with all explanations can be found in this repository: MOFA_to_COSMOS.md (https://github.com/saezlab/Factor_COSMOS/blob/main/MOFA_to_COSMOS.md)

Sample specific analysis can be found in: Cell_line_MOFA_space.md (https://github.com/saezlab/Factor_COSMOS/blob/main/Cell_line_MOFA_space.md)

Pathway control analysis can be found in: pathway_control_analysis_moon.md (https://github.com/saezlab/Factor_COSMOS/blob/main/pathway_control_analysis_moon.md)

## If you wish to run it locally:

Clone the repository locally

open NCI60.Rproj with Rstudio (then all the file path should be properlly working)

Open the .rmd file you are interested in (ex: MOFA_to_COSMOS.rmd)

Once the required packages are installed, and MOFA is set up with python, everything should run smoothly

(If you wish to run the part of COSMOS that relied on CARNIVAL, you will need to install the cplex IBM solver too and set the path in the script)

Please cite: https://www.embopress.org/doi/full/10.15252/msb.20209730