"""
replication_analysis_rq3.py

The exact confirmatory analysis for RQ3 (frozen neural-network-model
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

(The preregistration originally also specified a sign-reversal RQ, on
mutual friend count and target distance; it was cut before registration,
and the RQ numbering was shifted down accordingly -- the old RQ4 is this
script's RQ3.)

This is meant to be run exactly once, after the preregistration is
registered and the replication cohort's alldata.csv exists. Point
REPLICATION_ALLDATA_PATH at that file and run this script.

The frozen model artifacts load_frozen_models() reads are produced by
train_frozen_models.py, which reproduces the exact data split, scaling,
architecture, and BASE_SEED+r training procedure used in
neural_network_models.ipynb (the GELU models) and linear_vs_nonlinear_
ablation.ipynb (the one linear-activation model), with model-saving added
-- see that script's header for details. Each named model (complete_gelu,
complete_linear, and every coalition) is saved as N_REPEATS=2 separate
fits, matching both notebooks' own repeat-and-average convention exactly:
scoring averages each repeat's per-subject Brier skill, not raw predicted
probabilities (see score_model_repeats()).

Every piece of this script -- the scoring function, the exact Shapley
formula, and the three inferential tests with Holm-Bonferroni correction --
was tested end-to-end (including real frozen models fit and saved in
exactly this format) on a held-out subsample of the original study's own
data before first use.
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
OUTPUT_PATH              = "replication_results_rq3.csv"
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
    """Apply one already-fit, frozen model to new (replication) data: scales
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


def score_model_repeats(repeats, scaler, scaler_cols, data, outcome=OUTCOME):
    """Scores every repeat of one variant separately, then averages the
    per-subject Brier skill scores across repeats -- matching
    neural_network_models.ipynb's own `pd.concat(reps, axis=1).mean(axis=1)`
    convention exactly (average the scores, not the raw predicted
    probabilities, and not until each repeat has gone through its own
    per-subject groupby)."""
    per_repeat = [score_model(model, cols, scaler, scaler_cols, data, outcome) for model, cols in repeats]
    return pd.concat(per_repeat, axis=1).mean(axis=1)


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


N_REPEATS = 2   # must match train_frozen_models.py / both notebooks' N_REPEATS


def load_frozen_models(artifact_dir=ARTIFACT_DIR):
    """Loads the frozen model artifacts train_frozen_models.py produces.
    Expected layout, all produced once and never refit:

      artifact_dir/
        scaler.pkl                 -- joblib-dumped sklearn StandardScaler,
                                       fit on the ORIGINAL study's training
                                       split only (same object used in
                                       neural_network_models.ipynb)
        scaler_cols.json           -- the exact numeric_features column list
                                       (and order) the scaler was fit on
        complete_gelu_r{0,1}.keras       + complete_gelu_r{0,1}_cols.json
        complete_linear_r{0,1}.keras     + complete_linear_r{0,1}_cols.json
        shapley_coalitions/
          <key>_r{0,1}.keras + <key>_r{0,1}_cols.json
            for 63 of the 2**6 = 64 subsets of SUBSTANTIVE_GROUPS (visual is
            implicit in every one); <key> = coalition_key(subset), i.e.
            "empty" or "_".join(sorted(subset)). Each named model has
            N_REPEATS=2 separately-fit copies, matching both training
            notebooks' own repeat-and-average convention (see
            score_model_repeats()). The one subset NOT saved under
            shapley_coalitions/ is the full set (all 6 groups): that
            coalition IS the complete model, not a separate fit of the same
            architecture -- exactly as neural_network_models.ipynb's own
            Shapley cell reuses complete_repeats for that corner rather than
            refitting it (see its ALL_SUBSETS loop). Aliasing it to
            complete_gelu here, rather than loading an independently-trained
            duplicate, is what makes the Shapley values satisfy their
            defining efficiency property (sum of the six shares == complete
            - empty) exactly instead of only approximately -- two
            separately-fit copies of the same architecture land in
            different local optima and do NOT agree closely enough on new
            data for that identity to hold.

    Returns a dict: {'scaler':, 'scaler_cols':, 'complete_gelu': [(model, cols), ...],
    'complete_linear': [(model, cols), ...], 'coalitions': {frozenset: [(model, cols), ...]}},
    each list having N_REPEATS entries.
    """
    import pickle
    from tensorflow.keras.models import load_model

    def _require(path):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found. load_frozen_models() expects the artifacts "
                "train_frozen_models.py produces; see this module's docstring "
                "and that script's header for the expected layout."
            )
        return path

    with open(_require(os.path.join(artifact_dir, "scaler.pkl")), "rb") as f:
        scaler = pickle.load(f)
    with open(_require(os.path.join(artifact_dir, "scaler_cols.json"))) as f:
        scaler_cols = json.load(f)

    def _load_repeats(stem):
        repeats = []
        for r in range(N_REPEATS):
            model = load_model(_require(os.path.join(artifact_dir, f"{stem}_r{r}.keras")))
            with open(_require(os.path.join(artifact_dir, f"{stem}_r{r}_cols.json"))) as f:
                cols = json.load(f)
            repeats.append((model, cols))
        return repeats

    complete_gelu = _load_repeats("complete_gelu")
    complete_linear = _load_repeats("complete_linear")

    coalitions = {}
    n = len(SUBSTANTIVE_GROUPS)
    full_set = frozenset(SUBSTANTIVE_GROUPS)
    for r_size in range(n + 1):
        for subset in itertools.combinations(SUBSTANTIVE_GROUPS, r_size):
            key = frozenset(subset)
            if key == full_set:
                coalitions[key] = complete_gelu   # alias, not a separate fit -- see docstring
                continue
            coalitions[key] = _load_repeats(os.path.join("shapley_coalitions", coalition_key(subset)))

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


def run_rq3_tests(complete_bss, linear_bss, shapley_df):
    """The three preregistered inferential procedures. complete_bss and
    linear_bss are per-subject Series (same index); shapley_df is the output
    of shapley_values(), one column per substantive group."""
    rows = []

    t_complete = stats.ttest_1samp(complete_bss, 0, alternative="greater")
    rows.append({"rq": "RQ3a", "variant": "complete_vs_chance", "test": "one_sided_greater",
                 "statistic": t_complete.statistic, "p_value": t_complete.pvalue,
                 "sig": t_complete.pvalue < ALPHA, "mean": complete_bss.mean(), "n": len(complete_bss)})

    diff = complete_bss.reindex(linear_bss.index) - linear_bss
    t_diff = stats.ttest_1samp(diff.dropna(), 0, alternative="greater")
    rows.append({"rq": "RQ3b", "variant": "gelu_minus_linear", "test": "one_sided_greater",
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
        rows.append({"rq": "RQ3c", "variant": r["group"], "test": "two_sided_holm",
                     "statistic": r["statistic"], "p_value": r["p_holm"],
                     "sig": r["sig"], "mean": r["mean"], "n": r["n"]})

    return pd.DataFrame(rows)


def main():
    frozen = load_frozen_models()
    repl = pd.read_csv(REPLICATION_ALLDATA_PATH)

    print(f"Scoring {repl['sid'].nunique()} replication participants against frozen models...")
    complete_bss = score_model_repeats(frozen["complete_gelu"], frozen["scaler"], frozen["scaler_cols"], repl)
    linear_bss = score_model_repeats(frozen["complete_linear"], frozen["scaler"], frozen["scaler_cols"], repl)

    full_set = frozenset(SUBSTANTIVE_GROUPS)
    coalition_bss = {full_set: complete_bss}   # same model as complete_gelu -- reuse, don't rescore
    for subset, repeats in frozen["coalitions"].items():
        if subset == full_set:
            continue
        coalition_bss[subset] = score_model_repeats(repeats, frozen["scaler"], frozen["scaler_cols"], repl)
    shapley_df = shapley_values(coalition_bss)

    results = run_rq3_tests(complete_bss, linear_bss, shapley_df)
    results.to_csv(OUTPUT_PATH, index=False)
    print(results.to_string(index=False))
    print(f"\nwrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
