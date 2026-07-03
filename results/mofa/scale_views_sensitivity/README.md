# MOFA Scale-Views Sensitivity

This folder contains the Reviewer R1.3 sensitivity check for the NCI60 MOFA model.

The current manuscript analysis selects `results/mofa/mofa_res_10factor.hdf5`, a max-10-factor MOFA run that retained 9 active factors. To test whether the original `scale_views = False` setting materially affected the factor structure, we trained one additional model with the same settings except:

- `scale_views = True`
- `factors = 10`
- output: `results/mofa/mofa_res_10factor_scale_views_true.hdf5`

Training command:

```bash
PYTHONPATH=/private/tmp/mofapy2_deps \
/Users/dugourd/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
scripts/mofa/mofa2.py \
results/mofa/scale_views_sensitivity/options_MOFA_scale_views_true.csv
```

Comparison command:

```bash
PYTHONPATH=/private/tmp/mofapy2_deps \
/Users/dugourd/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
scripts/mofa/compare_scale_views_models.py
```

The temporary `PYTHONPATH` points to an isolated install of `mofapy2` and `h5py` used for the revision run.
