# Code for "The graph theoretic determinants of friendship selection"

This repository contains the analysis code for Thornton, "The graph theoretic determinants of friendship selection." The corresponding data (raw, processed, and cached intermediate results) are archived separately on Zenodo at `[DOI]`.

## Pipeline

Run in this order to reproduce the full analysis from the rawest available data:

1. **`expand_data.R`** — reads `net_data.csv` (the raw per-trial platform export), computes the core set of graph theoretic features, and writes `fulldata.csv`.
2. **`expand_extra_features.R`** — reads `net_data.csv` again plus the stimulus graph pool (`graphs.csv`), computes a further set of graph theoretic features, and writes `extra_features.csv`.
3. **`merge_features.R`** — merges `fulldata.csv` and `extra_features.csv` into `alldata.csv`, the single file both of the scripts below read. Also fills the one known data quality issue (local clustering coefficient is undefined, not missing, for nodes with degree < 2) and reorders columns into one coherent scheme.
4. **`analyze_nets.R`** — all R-side analyses: zero-order, marginal-over-visual, and unique/partial linear probability models; the individual-differences PCA, its split-half reliability, and the bimodality test; produces `combined_regression.pdf`, `pca_scree_loadings.pdf`, and `pca_demographics.pdf`.
5. **`neural_network_models.ipynb`** — the main neural network analyses: zero-order, unique (leave-one-group-out), marginal-over-visual, and exact Shapley decomposition of each feature group's contribution to predicting choice. Produces `combined_brier_skill.pdf`.
6. **`linear_vs_nonlinear_ablation.ipynb`** — the linear-vs-nonlinear activation comparison and the Shapley decomposition of the nonlinear advantage across feature groups. Produces `linear_vs_nonlinear.pdf` and `shapley_nonlinear_advantage.pdf`.
7. **`conditional_logit_robustness.R`** — a robustness check on the complete/partial model: a single pooled conditional logit across every trial (McFadden-style forced choice among the trial's candidates), compared against the mean partial linear-probability coefficients. Produces `conditional_logit_results.csv`.

Steps 1-3 are slow (they iterate over every trial) and their output is large. If you only want to rerun the analyses themselves, start from `alldata.csv` directly (archived alongside the raw files on Zenodo) and skip to step 4.

## Requirements

- R 4.3.2, with `igraph`, `RColorBrewer`, `psych`, `MASS`, `diptest`, `survival`, `RJSONIO`, `lmerTest`.
- Python 3.12, with TensorFlow 2.21, pandas, numpy, scipy, scikit-learn, matplotlib.

## Notes

- `expand_data.R` and `expand_extra_features.R` both read the shared stimulus graph pool at `graphs.csv`; this file is not specific to this study and is also used by other experiments on the same platform.
- Each neural network variant is refit multiple times from different random initializations and averaged (2 repeats in most cases; 10 for the leave-one-group-out/unique contribution analysis, which showed meaningfully lower repeat-to-repeat reliability at 2). See the notebooks for the full repeat-reliability diagnostic.
- `conditional_logit_robustness.R` fits one population-average model across all subjects pooled together, not per subject — a per-subject conditional logit would risk the same quasi-complete-separation problem that motivated using linear probability models over per-subject logistic regression in the first place.
