# Detecting and quantifying rare sex in natural populations (Pieszko et al., 2025)

Code accompanying the manuscript (add DOI). The directories contain:

- `sim_pipeline/`
  - Snakemake pipeline used to run simulations using models of asexuality (`scripts/model_MP.slim` and `scripts/model_CF.slim`).
- `sexsign_functions/`
  - Functions for calculating statistics and summarising tree composition.
- `stats/`
  - Scripts to calculate statistics and summarise tree composition across the space of simulated scenarios.
- `plots/`
  - Scripts to generate all plots.
- `variation_analysis/`
  - Analyses described in *Predictability and discriminability of different scenarios*.
- `validation/`
  - Validation of simulation models with analytical results (as described in Methods and SI Appendix, section A). Includes a separate Snakemake pipeline to run obligate asexual scenarios under the MP and CF models in turn.