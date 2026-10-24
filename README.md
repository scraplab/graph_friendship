# Code for "The graph theoretic determinants of friendship selection"

This repository contains the analysis code for Thornton, "The graph theoretic determinants of friendship selection." The corresponding data (raw, processed, and cached intermediate results) are archived separately on Zenodo at [doi.org/10.5281/zenodo.23628825](https://doi.org/10.5281/zenodo.23628825).

**`task_frontend/`** has the participant-facing task code (the HTML/CSS/JS that presented the network choice task on MySocialBrain.org), shared for transparency rather than as a runnable copy -- see its own README.

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
2. **`expand_all_features.R`** — reads `net_data.csv` (the raw per-trial platform export) and `graphs.csv`, computes every graph theoretic feature used in the study, and writes `alldata.csv`. This consolidates what used to be three separate scripts (feature computation, a second round of additional features, and a merge step) into one; see the comment at the top of the script for exactly how it relates to that earlier version, and for how its output was verified against the original three-script pipeline on a reconstructed sample.
3. **`analyze_nets.R`** — all R-side analyses: zero-order, marginal-over-visual, and unique/partial linear probability models; the individual-differences PCA, its split-half reliability, and the bimodality test. Produces `combined_regression.pdf`, `pca_scree_loadings.pdf`, and `pca_demographics.pdf`. Runtime: the per-subject coefficient loop is the slow part (tens of minutes for 8,614 subjects).
4. **`neural_network_models.ipynb`** — the main neural network analyses: zero-order, unique (leave-one-group-out), marginal-over-visual, and exact Shapley decomposition of each feature group's contribution to predicting choice. Produces `combined_brier_skill.pdf`. Runtime: by far the slowest step in the whole pipeline (hours), since it fits on the order of 100+ model variants.
5. **`linear_vs_nonlinear_ablation.ipynb`** — the linear-vs-nonlinear activation comparison and the Shapley decomposition of the nonlinear advantage across feature groups. Produces `linear_vs_nonlinear.pdf` and `shapley_nonlinear_advantage.pdf`. Also slow (252 model fits for the nonlinear-advantage Shapley decomposition alone).
6. **`conditional_logit_robustness.R`** — a robustness check on the complete/partial model: a single pooled conditional logit across every trial (McFadden-style forced choice among that trial's candidates), compared against the mean partial linear-probability coefficients. Produces `conditional_logit_results.csv` and `conditional_logit_diagnostics.pdf` (a side-by-side coefficient comparison and a concordance scatterplot, colored by feature category using the same palette as `combined_regression.pdf` — currently for internal inspection, not yet a manuscript figure). Runtime: a few minutes.

Steps 1-2 are slow and their output is large. If you only want to rerun the analyses themselves, start from `alldata.csv` directly (archived alongside the raw files on Zenodo) and skip to step 3.

### Replication planning (not part of the original paper's pipeline)

- **`power_analysis_prereg.R`** and **`power_analysis_prereg.py`** — power analyses for every analysis specified in the replication preregistration (see `preregistration_draft.html` in the manuscript folder), used to determine what sample size the planned replication needs. These read the same cached results as the main pipeline but are a separate, later addition supporting the *replication*, not the original study. `power_analysis_prereg.R` covers profile replication (Wasserstein distance + permutation, and Pearson correlation), the sign-reversal checks, and the PCA generalization tests; `power_analysis_prereg.py` covers the frozen-model transfer tests. Both print a "HEADLINE" summary at the end specifically for N = 300, the replication's planned sample size, pulled from the same grid each script already computes (which also spans a wider range of N, in case the planned size changes). The Wasserstein+permutation simulation in the R script is the slow part (tens of minutes); everything else runs in well under a minute.
- **`replication_analysis_rq1_rq3.R`** and **`replication_analysis_rq4.py`** — the actual confirmatory analyses to run once the preregistration is registered and the replication cohort's own `alldata.csv` exists (point `REPLICATION_ALLDATA_PATH` at it near the top of each script). These implement exactly what `preregistration_draft.html`'s "Statistical models" section specifies: RQ1 (profile replication: Wasserstein distance + permutation, and Pearson correlation, a true two-sample version of the test in the power script above), RQ2 (two-sided sign-reversal t-tests), RQ3 (out-of-sample PCA reconstruction + noise ceiling, and Tucker congruence on a fresh fit) in the R script; RQ4 (frozen complete/linear/Shapley-coalition model transfer, with Holm-Bonferroni correction across the six Shapley tests) in the Python script. Both were tested end-to-end on a real, disjoint held-out subsample of the original study's own data standing in for a genuine replication cohort (not just unit-tested in isolation), which caught and fixed two bugs: a column-schema mismatch when combining the three RQs' result tables in the R script, and a subtler one in the Python script -- the Shapley decomposition's "full coalition" (all six substantive groups) must be the exact same fitted model as the complete model, not an independently retrained copy of the same architecture, or the Shapley shares' defining property (they sum exactly to the complete-vs-chance gap) only holds approximately instead of exactly.
- **`train_frozen_models.py`** — produces the frozen model artifacts `replication_analysis_rq4.py`'s `load_frozen_models()` reads: the complete GELU model, the complete linear model, the `StandardScaler` they were trained with, and the 63 non-trivial feature-group coalition models needed for the Shapley decomposition (the 64th, all six substantive groups, is the complete model itself -- not a separate fit; see that script's docstring). This is not a new analysis -- the data split, scaling, architecture, and `BASE_SEED+r` training procedure are copied verbatim from `neural_network_models.ipynb` (the GELU models) and `linear_vs_nonlinear_ablation.ipynb` (the one linear-activation model); the only addition is saving each fitted model to disk, since neither notebook otherwise persists anything but per-subject Brier skill scores. Resumable (skips any model whose files already exist).

  Already run once on the full N = 8,614 sample (took about 15.5 hours on a CPU-only machine). As a sanity check, the freshly refit complete model's Brier skill (GELU: .0554, linear: .0452, averaged across repeats) landed close to the values already cached from the original notebook runs (GELU: .0561-.0563, linear: .0450) -- consistent with the same procedure, not a new or different model.

  **The resulting artifacts (31MB, 263 files, including `bss_sanity_check.csv` with the full per-fit numbers) are archived on Zenodo at [doi.org/10.5281/zenodo.23628825](https://doi.org/10.5281/zenodo.23628825) alongside the raw/processed data, not in this git repository** -- they're training output, not code. Place the Zenodo deposit's `frozen_models/` folder directly alongside this code (same convention as `alldata.csv` -- see Data and working directory above) before running `replication_analysis_rq4.py`, or point its `ARTIFACT_DIR` at wherever you put it.

## Notes

- `generate_networks.R` and `expand_all_features.R` both read the shared stimulus graph pool at `graphs.csv`; this file is not specific to this study and is also used by other experiments on the same platform.
- Each neural network variant is refit multiple times from different random initializations and averaged (2 repeats in most cases; 10 for the leave-one-group-out/unique contribution analysis, which showed meaningfully lower repeat-to-repeat reliability at 2). See the notebooks for the full repeat-reliability diagnostic.
- `conditional_logit_robustness.R` fits one population-average model across all subjects pooled together, not per subject — a per-subject conditional logit would risk the same quasi-complete-separation problem that motivated using linear probability models over per-subject logistic regression in the first place.
- If you cloned or copied this repository from somewhere that also has a `.git` folder synced via Dropbox or a similar service, make sure you're working against only one copy — two working directories sharing the same underlying `.git` state (via file sync rather than a real git remote) can corrupt the index if both are touched around the same time. Push this repository to GitHub and treat that as the single source of truth going forward.
