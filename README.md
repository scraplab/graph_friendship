# Code for "The graph theoretic determinants of friendship selection"

This repository contains the analysis code for Thornton, "The graph theoretic determinants of friendship selection." The corresponding data (raw, processed, and cached intermediate results) are archived separately on Zenodo at `[DOI]`.

## Setup

### R

Requires R 4.3.2 (later 4.x versions should work; this is what was used). Install dependencies once:

```r
install.packages(c("igraph", "RColorBrewer", "psych", "MASS", "diptest",
                   "survival", "RJSONIO", "lmerTest", "pwr"))
```

### Python

Requires Python 3.12. The exact environment used:

```
tensorflow==2.21.0
scikit-learn==1.9.1
numpy==2.5.3
pandas==3.0.5
scipy==1.18.1
matplotlib
jupyter
```

With conda:

```bash
conda create -n graph_friendship python=3.12
conda activate graph_friendship
pip install tensorflow==2.21.0 scikit-learn==1.9.1 numpy==2.5.3 pandas==3.0.5 scipy==1.18.1 matplotlib jupyter
```

Exact pinned versions matter less for the R-side scripts than for the neural network notebooks; TensorFlow's numerics are not bit-identical across versions, so an exact match is the safest way to reproduce the cached Brier skill scores precisely, though results should be qualitatively identical with any reasonably recent 2.x release.

### Data and working directory

None of the scripts take a data-directory argument; each one has a hardcoded `setwd()` (R) or assumes the working directory (Python) at the top, pointing at the path this project happened to live at. **Before running anything**, either:

- edit that line in each script to point at wherever you've put the data (simplest), or
- place the files from the Zenodo data deposit directly alongside this code and remove/comment out the `setwd()` calls.

`generate_networks.R` and `expand_all_features.R` additionally read a shared stimulus-graph file (`graphs.csv`) via a relative path (`../../graphs.csv`), reflecting this project's original folder depth (two levels under a shared `net/` directory with other experiments on the same platform). Adjust that path too if your layout differs.

## Pipeline

Run in this order to reproduce the full analysis from the rawest available data:

1. **`generate_networks.R`** — generates the stimulus graph pool: samples 100,000 random geometric graphs, filters to the 10,000 meeting the study's realism criteria (connected, density < .25, assortativity > 0, modularity > .1), and computes each graph's node- and network-level properties. Produces the `graphs.csv` stimulus pool (written in 10 chunk files to avoid an import timeout; these were concatenated into the single `graphs.csv` the rest of the pipeline reads). You do not need to rerun this unless you want a fresh stimulus pool — `graphs.csv` is static platform infrastructure, not something regenerated per analysis.
2. **`expand_all_features.R`** — reads `net_data.csv` (the raw per-trial platform export) and `graphs.csv`, computes every graph theoretic feature used in the study, and writes `alldata.csv`. This consolidates what used to be three separate scripts (feature computation, a second round of additional features, and a merge step) into one; see the comment at the top of the script for exactly how it relates to that earlier version. **Note:** the script's own header comment still says "NOT YET RUN" as of this repository snapshot — confirm it has in fact been verified against the original three-script output before relying on that comment being stale.
3. **`analyze_nets.R`** — all R-side analyses: zero-order, marginal-over-visual, and unique/partial linear probability models; the individual-differences PCA, its split-half reliability, and the bimodality test. Produces `combined_regression.pdf`, `pca_scree_loadings.pdf`, and `pca_demographics.pdf`. Runtime: the per-subject coefficient loop is the slow part (tens of minutes for 8,614 subjects).
4. **`neural_network_models.ipynb`** — the main neural network analyses: zero-order, unique (leave-one-group-out), marginal-over-visual, and exact Shapley decomposition of each feature group's contribution to predicting choice. Produces `combined_brier_skill.pdf`. Runtime: by far the slowest step in the whole pipeline (hours), since it fits on the order of 100+ model variants.
5. **`linear_vs_nonlinear_ablation.ipynb`** — the linear-vs-nonlinear activation comparison and the Shapley decomposition of the nonlinear advantage across feature groups. Produces `linear_vs_nonlinear.pdf` and `shapley_nonlinear_advantage.pdf`. Also slow (252 model fits for the nonlinear-advantage Shapley decomposition alone).
6. **`conditional_logit_robustness.R`** — a robustness check on the complete/partial model: a single pooled conditional logit across every trial (McFadden-style forced choice among that trial's candidates), compared against the mean partial linear-probability coefficients. Produces `conditional_logit_results.csv` and `conditional_logit_diagnostics.pdf` (a side-by-side coefficient comparison and a concordance scatterplot, colored by feature category using the same palette as `combined_regression.pdf` — currently for internal inspection, not yet a manuscript figure). Runtime: a few minutes.

Steps 1-2 are slow and their output is large. If you only want to rerun the analyses themselves, start from `alldata.csv` directly (archived alongside the raw files on Zenodo) and skip to step 3.

### Replication planning (not part of the original paper's pipeline)

- **`power_analysis_prereg.R`** and **`power_analysis_prereg.py`** — power analyses for every analysis specified in the replication preregistration (see `preregistration_draft.html` in the manuscript folder), used to determine what sample size the planned replication needs. These read the same cached results as the main pipeline but are a separate, later addition supporting the *replication*, not the original study. `power_analysis_prereg.R` covers profile replication (Wasserstein distance + permutation, and Pearson correlation), the sign-reversal checks, and the PCA generalization tests; `power_analysis_prereg.py` covers the frozen-model transfer tests. The Wasserstein+permutation simulation in the R script is the slow part (tens of minutes); everything else runs in well under a minute.

## Notes

- `generate_networks.R` and `expand_all_features.R` both read the shared stimulus graph pool at `graphs.csv`; this file is not specific to this study and is also used by other experiments on the same platform.
- Each neural network variant is refit multiple times from different random initializations and averaged (2 repeats in most cases; 10 for the leave-one-group-out/unique contribution analysis, which showed meaningfully lower repeat-to-repeat reliability at 2). See the notebooks for the full repeat-reliability diagnostic.
- `conditional_logit_robustness.R` fits one population-average model across all subjects pooled together, not per subject — a per-subject conditional logit would risk the same quasi-complete-separation problem that motivated using linear probability models over per-subject logistic regression in the first place.
- If you cloned or copied this repository from somewhere that also has a `.git` folder synced via Dropbox or a similar service, make sure you're working against only one copy — two working directories sharing the same underlying `.git` state (via file sync rather than a real git remote) can corrupt the index if both are touched around the same time. Push this repository to GitHub and treat that as the single source of truth going forward.
