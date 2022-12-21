from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import sys

options_MOFA = sys.argv[1]

options_MOFA_dt = pd.read_csv(options_MOFA, sep=",")

likelihoods = []
for i in range(1, options_MOFA_dt['length'][0]+1): 
	likelihoods.append(options_MOFA_dt["likelihood_"+str(i)][0])


# initialise the entry point
ent = entry_point()

data_dt = pd.read_csv(options_MOFA_dt['input_file'][0], sep=",")
ent.set_data_options(
    scale_groups = options_MOFA_dt["scale_groups"][0], 
    scale_views = options_MOFA_dt["scale_views"][0]
)
ent.set_data_df(
	data_dt, 
	likelihoods = likelihoods
)
ent.set_model_options(
    factors = options_MOFA_dt["factors"][0], 
    spikeslab_weights = options_MOFA_dt["spikeslab_weights"][0],
    spikeslab_factors = options_MOFA_dt["spikeslab_factors"][0],
    ard_factors = options_MOFA_dt["ard_factors"][0],
    ard_weights = options_MOFA_dt["ard_weights"][0]
)
ent.set_train_options(
    iter = options_MOFA_dt["iter"][0], 
    convergence_mode = options_MOFA_dt["convergence_mode"][0], 
    startELBO = options_MOFA_dt["startELBO"][0], 
    freqELBO = options_MOFA_dt["freqELBO"][0], 
    dropR2 = options_MOFA_dt["dropR2"][0], 
    gpu_mode = options_MOFA_dt["gpu_mode"][0], 
    verbose = options_MOFA_dt["verbose"][0], 
    seed = options_MOFA_dt["seed"][0]
)
ent.build()
ent.run()
# Save the output
ent.save(outfile=options_MOFA_dt["output_file"][0])

