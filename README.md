# NCI60_cosmos
Formatting NCI60 data into cosmos ready inputs and generation of testable hypothesis connecting cell-line specific TF and metabolic deregulations.

You can install Rstudio and open the NCI60.Rproj with it.

Then simply open the scritps in the Rstudio and run them.

data/cosmso/cosmos_inputs.RData contains cosmos-ready inputs for each NCI60 cell lines.

scripts/run_cosmos.R to run cosmos analysis on a specific NCI60 cell line. Just pick a cell line based on it's name.

/!\ If you do not have access to the IBM-CPLEX solver, we advise you to instead run our solver-free pipeline: scripts/net_compr_decouplRnival.R /!\

scripts/prepare_comsos_inputs.R shows how to prepare multi-omic inputs for cosmos, and how to estiamte TF activities from transcriptomic data with decoupleR (https://github.com/saezlab/decoupler)

install cosmos here: https://github.com/saezlab/cosmosR

Have fun !
