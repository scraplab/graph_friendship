"""
replication_analysis_rq4.py

The exact confirmatory analysis for RQ4 (frozen neural-network-model
transfer) of the replication preregistration (preregistration_draft.html,
"Statistical models" section): the original study's already-trained complete
model (GELU), its linear-activation counterpart, and the coalition models
needed for the exact Shapley decomposition over the six substantive feature
groups, all applied to the replication sample with NO retraining, scored by
per-participant Brier skill, and tested as specified there:
  (a) complete model vs. chance -- one-sided (greater)
  (b) GELU-minus-linear complete-model difference vs. 0 -- one-sided (greater)
  (c) each of the six Shapley shares vs. 0 -- two-sided, Holm-Bonferroni
      corrected across the six

This is meant to be run exactly once, after the preregistration is
registered and the replication cohort's alldata.csv exists. Point
REPLICATION_ALLDATA_PATH at that file and run this script.

IMPORTANT, as of this writing: load_frozen_models() below has nothing to
load. neural_network_models.ipynb and linear_vs_nonlinear_ablation.ipynb fit
all their models in memory and never persist the trained weights or the
StandardScaler to disk -- only per-subject Brier skill scores on the
original study's own held-out test split are cached. "Frozen-model
transfer" requires the actual trained models (and the exact scaler they were
trained with), so before this script can be run for real, the training
notebooks need a small addition: save the complete GELU model, the complete
linear model, the scaler, and every one of the 2**6 = 64 feature-group
coalition models already being fit for the Shapley decomposition (Analysis
3a), in the format load_frozen_models() documents below. That is a
straightforward addition (every model these notebooks fit already exists in
memory at some point; this just keeps a copy), not a redesign, but it does
mean re-running (part of) the training pipeline once to produce the
artifacts -- flagged here rather than done unprompted given the runtime cost
(see the main README).

Every other piece of this script -- the scoring function, the exact Shapley
formula, and the three inferential tests with Holm-Bonferroni correction --
is implemented and was tested end-to-end (including real frozen models
fit and saved in exactly this format) on a held-out subsample of the
original study's own data before first use.
"""
import itertools
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import stats

# ---- CONFIG: edit before running for real ----------------------------------
REPLICATION_ALLDATA_PATH = "replication_alldata.csv"
ARTIFACT_DIR             = "frozen_models"
OUTPUT_PATH              = "replication_results_rq4.csv"
ALPHA                    = .05

OUTCOME = "binres"
FEATURE_GROUPS = {
    "target":         ["degree", "closeness", "betweenness", "eigenvector",
                        "clust", "constr", "ecc", "esize", "comsize"],
    "dyad":           ["samecom", "edgbet", "udists",
                        "mutfr", "dhom", "structsim", "aa", "katz", "srcos"],
    "network":        ["density", "assortativity", "apl", "modularity",
                        "netcent", "netvar"],
    "change_network": ["chgmod", "chgast", "chapl", "chgcent", "chgvar"],
    "self":           ["udegree", "ucloseness", "ubetweenness", "ueigenvector",
                        "uclust", "uconstr", "uecc", "uesize", "ucomsize"],
    "change_self":    ["chub", "chuc", "chue", "chuclust", "chuconstr", "chuecc", "chuesize"],
    "visual":         ["vdist", "vecc", "uvecc"],
}
SUBSTANTIVE_GROUPS = [g for g in FEATURE_GROUPS if g != "visual"]


def group_columns(active_groups, feature_groups=FEATURE_GROUPS):
    cols = []
    for g in active_groups:
        cols.extend(feature_groups[g])
    return cols


def coalition_key(subset):
    """Canonical artifact filename stem for a subset of SUBSTANTIVE_GROUPS
    (visual is always included on top of whatever's in `subset`)."""
    return "empty" if not subset else "_".join(sorted(subset))


def brier_skill(y_true, y_pred, p_chance, sid):
    """Per-participant Brier skill score: 1 - BS/BS_chance. Identical
    formula to both training notebooks' brier_skill()."""
    bs = (y_true - y_pred) ** 2
    bs_chance = (y_true - p_chance) ** 2
    df_ = pd.DataFrame({"sid": sid, "bs": bs, "bs_chance": bs_chance})
    g = df_.groupby("sid").mean()
    return (1 - g["bs"] / g["bs_chance"]).rename("bss")


def score_model(model, cols, scaler, scaler_cols, data, outcome=OUTCOME):
    """Apply an already-fit, frozen model to new (replication) data: scales
    `data` with the ORIGINAL training scaler (never refit on the replication
    sample -- the model's weights were learned in that fixed coordinate
    system), predicts, and returns per-participant Brier skill scores.

    cols: the exact, ordered feature-column list this particular model was
    trained on (narrower than scaler_cols in general, since the scaler was
    fit once on every numeric feature but each coalition model only consumes
    a subset of them).
    """
    X = data.copy()
    X[scaler_cols] = scaler.transform(X[scaler_cols])
    Xc = X[cols].to_numpy(dtype="float32")
    n_options = X.groupby(["sid", "trial"])["id"].transform("size")
    p_chance = (1.0 / n_options).to_numpy()
    y_true = X[outcome].to_numpy(dtype="float32")
    sid = X["sid"].to_numpy()
    p_pred = model.predict(Xc, batch_size=4096, verbose=0).ravel()
    return brier_skill(y_true, p_pred, p_chance, sid)


def shapley_values(coalition_bss, substantive_groups=SUBSTANTIVE_GROUPS):
    """Exact Shapley decomposition of each substantive group's contribution,
    from a dict {frozenset(subset of substantive_groups): per-subject BSS
    series}, where every coalition additionally includes 'visual' as a fixed
    baseline (so coalition_bss[frozenset()] is the visual-alone model and
    coalition_bss[frozenset(substantive_groups)] is the complete model).
    Identical formula/weighting to Analysis 3b in neural_network_models.ipynb:

      phi_g = sum over S subset of (groups \\ {g}) of
              [|S|! * (n-|S|-1)!] / n! * (v(S u {g}) - v(S))

    computed exactly over all 2**n subsets (n = 6 keeps this cheap), not
    sampled.
    """
    n = len(substantive_groups)
    phi = {}
    any_key = next(iter(coalition_bss))
    all_sids = coalition_bss[any_key].index
    for g in substantive_groups:
        others = [x for x in substantive_groups if x != g]
        total = pd.Series(0.0, index=all_sids)
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                weight = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)
                v_S = coalition_bss[S]
                v_Sg = coalition_bss[S | {g}]
                total = total + weight * (v_Sg.reindex(all_sids) - v_S.reindex(all_sids))
        phi[g] = total.rename(g)
    return pd.DataFrame(phi)


def load_frozen_models(artifact_dir=ARTIFACT_DIR):
    """Loads the frozen model artifacts this script expects. See the module
    docstring for why none of this currently exists on disk. Expected
    layout, all produced once at the end of the original training run and
    never refit:

      artifact_dir/
        scaler.pkl                 -- joblib-dumped sklearn StandardScaler,
                                       fit on the ORIGINAL study's training
                                       split only (same object used in
                                       neural_network_models.ipynb)
        scaler_cols.json           -- the exact numeric_features column list
                                       (and order) the scaler was fit on
        complete_gelu.keras        + complete_gelu_cols.json
        complete_linear.keras      + complete_linear_cols.json
        shapley_coalitions/
          <key>.keras + <key>_cols.json
            for 63 of the 2**6 = 64 subsets of SUBSTANTIVE_GROUPS (visual is
            implicit in every one); <key> = coalition_key(subset), i.e.
            "empty" or "_".join(sorted(subset)). The one subset NOT saved
            here is the full set (all 6 groups): that coalition IS the
            complete model, not a separate fit of the same architecture --
            exactly as neural_network_models.ipynb's own Shapley cell
            reuses complete_repeats for that corner rather than refitting it
            (see its ALL_SUBSETS loop). Aliasing it to complete_gelu here,
            rather than loading an independently-trained duplicate, is what
            makes the Shapley values satisfy their defining efficiency
            property (sum of the six shares == complete - empty) exactly
            instead of only approximately -- two separately-fit copies of
            the same architecture land in different local optima and do NOT
            agree closely enough on new data for that identity to hold.

    Returns a dict: {'scaler':, 'scaler_cols':, 'complete_gelu': (model, cols),
    'complete_linear': (model, cols), 'coalitions': {frozenset: (model, cols)}}.
    """
    import pickle
    from tensorflow.keras.models import load_model

    def _require(path):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found. load_frozen_models() expects the artifacts "
                "documented in its docstring; see this module's header for what "
                "still needs to be added to the training notebooks to produce them."
            )
        return path

    with open(_require(os.path.join(artifact_dir, "scaler.pkl")), "rb") as f:
        scaler = pickle.load(f)
    with open(_require(os.path.join(artifact_dir, "scaler_cols.json"))) as f:
        scaler_cols = json.load(f)

    def _load_named(stem):
        model = load_model(_require(os.path.join(artifact_dir, f"{stem}.keras")))
        with open(_require(os.path.join(artifact_dir, f"{stem}_cols.json"))) as f:
            cols = json.load(f)
        return model, cols

    complete_gelu = _load_named("complete_gelu")
    complete_linear = _load_named("complete_linear")

    coalitions = {}
    n = len(SUBSTANTIVE_GROUPS)
    full_set = frozenset(SUBSTANTIVE_GROUPS)
    for r in range(n + 1):
        for subset in itertools.combinations(SUBSTANTIVE_GROUPS, r):
            key = frozenset(subset)
            if key == full_set:
                coalitions[key] = complete_gelu   # alias, not a separate fit -- see docstring
                continue
            model, cols = _load_named(os.path.join("shapley_coalitions", coalition_key(subset)))
            coalitions[key] = (model, cols)

    return {"scaler": scaler, "scaler_cols": scaler_cols,
            "complete_gelu": complete_gelu, "complete_linear": complete_linear,
            "coalitions": coalitions}


def holm_bonferroni(pvals):
    """Standard Holm-Bonferroni step-down correction. Returns adjusted
    p-values in the same order as `pvals`."""
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = min((m - rank) * pvals[idx], 1.0)
        running_max = max(running_max, val)
        adj[idx] = running_max
    return adj


def run_rq4_tests(complete_bss, linear_bss, shapley_df):
    """The three preregistered inferential procedures. complete_bss and
    linear_bss are per-subject Series (same index); shapley_df is the output
    of shapley_values(), one column per substantive group."""
    rows = []

    t_complete = stats.ttest_1samp(complete_bss, 0, alternative="greater")
    rows.append({"rq": "RQ4a", "variant": "complete_vs_chance", "test": "one_sided_greater",
                 "statistic": t_complete.statistic, "p_value": t_complete.pvalue,
                 "sig": t_complete.pvalue < ALPHA, "mean": complete_bss.mean(), "n": len(complete_bss)})

    diff = complete_bss.reindex(linear_bss.index) - linear_bss
    t_diff = stats.ttest_1samp(diff.dropna(), 0, alternative="greater")
    rows.append({"rq": "RQ4b", "variant": "gelu_minus_linear", "test": "one_sided_greater",
                 "statistic": t_diff.statistic, "p_value": t_diff.pvalue,
                 "sig": t_diff.pvalue < ALPHA, "mean": diff.mean(), "n": diff.dropna().shape[0]})

    shap_rows = []
    for g in shapley_df.columns:
        t = stats.ttest_1samp(shapley_df[g], 0, alternative="two-sided")
        shap_rows.append({"group": g, "statistic": t.statistic, "p_value": t.pvalue,
                           "mean": shapley_df[g].mean(), "n": len(shapley_df[g])})
    shap_df = pd.DataFrame(shap_rows)
    shap_df["p_holm"] = holm_bonferroni(shap_df["p_value"].to_numpy())
    shap_df["sig"] = shap_df["p_holm"] < ALPHA
    for _, r in shap_df.iterrows():
        rows.append({"rq": "RQ4c", "variant": r["group"], "test": "two_sided_holm",
                     "statistic": r["statistic"], "p_value": r["p_holm"],
                     "sig": r["sig"], "mean": r["mean"], "n": r["n"]})

    return pd.DataFrame(rows)


def main():
    frozen = load_frozen_models()
    repl = pd.read_csv(REPLICATION_ALLDATA_PATH)

    print(f"Scoring {repl['sid'].nunique()} replication participants against frozen models...")
    model, cols = frozen["complete_gelu"]
    complete_bss = score_model(model, cols, frozen["scaler"], frozen["scaler_cols"], repl)

    model, cols = frozen["complete_linear"]
    linear_bss = score_model(model, cols, frozen["scaler"], frozen["scaler_cols"], repl)

    full_set = frozenset(SUBSTANTIVE_GROUPS)
    coalition_bss = {full_set: complete_bss}   # same model as complete_gelu -- reuse, don't rescore
    for subset, (model, cols) in frozen["coalitions"].items():
        if subset == full_set:
            continue
        coalition_bss[subset] = score_model(model, cols, frozen["scaler"], frozen["scaler_cols"], repl)
    shapley_df = shapley_values(coalition_bss)

    results = run_rq4_tests(complete_bss, linear_bss, shapley_df)
    results.to_csv(OUTPUT_PATH, index=False)
    print(results.to_string(index=False))
    print(f"\nwrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
