"""
train_frozen_models.py

Produces the frozen model artifacts replication_analysis_rq4.py's
load_frozen_models() expects: the complete GELU model, the complete linear
model, the StandardScaler they were trained with, and the 63 non-trivial
feature-group coalition models needed for the exact Shapley decomposition
(the 64th, all six substantive groups, IS the complete model -- see
replication_analysis_rq4.py's docstring).

This is NOT a new analysis. Every model fit here, and the exact procedure
used to fit it (data load, the participant-level 70/10/20 split, z-scoring,
FEATURE_GROUPS, build_model, fit_and_score, the BASE_SEED+r seed convention,
N_LAYERS, epochs/batch_size/patience), is copied verbatim from
neural_network_models.ipynb (the gelu models) and linear_vs_nonlinear_
ablation.ipynb (the linear complete model -- the one intentional difference
there is activation='linear'). Nothing about the architecture, split, or
training procedure is changed; the only addition is saving each fitted
model (and the scaler) to disk, since neither notebook currently persists
anything but per-subject Brier skill scores.

Resumable: every model's save path is checked before fitting, so an
interrupted run can simply be restarted and will skip whatever already
completed. This matters here because the full set (1 complete-gelu + 1
complete-linear + 63 coalitions, x2 repeats each = 130 fits on the full
N=8,614 sample) is expected to take on the order of hours on this
CPU-only machine, consistent with neural_network_models.ipynb's own
documented runtime for a comparable fit count.

Per-subject Brier skill scores for every fit are also written out
(frozen_models/bss_sanity_check.csv) purely as a sanity check against the
existing cached_results (nn_results_cache.npz, linear_vs_nonlinear_
results.csv) -- NOT to replace them. Those cached files back the already-
written manuscript and are left untouched; see the chat response this
script's commit refers back to for the comparison.
"""
import itertools
import json
import os
import pickle

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from tensorflow.keras import Model
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.layers import BatchNormalization, Dense, Input

ARTIFACT_DIR = "frozen_models"
os.makedirs(os.path.join(ARTIFACT_DIR, "shapley_coalitions"), exist_ok=True)

np.random.seed(1)
tf.random.set_seed(1)

# ---- load data -- identical to both notebooks -------------------------------
print("loading alldata.csv...")
df = pd.read_csv("alldata.csv")
print(f"  {len(df)} rows, {df.sid.nunique()} participants")

# ---- participant-level partition -- identical to both notebooks -------------
usid = np.asarray(df.sid.unique()).copy()
np.random.shuffle(usid)
n = usid.size
n_train = int(round(0.70 * n))
n_val = int(round(0.10 * n))
utrain, uval, utest = np.split(usid, [n_train, n_train + n_val])
traindat = df.loc[df.sid.isin(utrain)].copy()
valdat = df.loc[df.sid.isin(uval)].copy()
testdat = df.loc[df.sid.isin(utest)].copy()
print(f"  train={len(utrain)} val={len(uval)} test={len(utest)} participants")

OUTCOME = "binres"

# ---- z-score normalization -- identical column list to both notebooks -------
numeric_features = [
    "age",
    "degree", "closeness", "betweenness", "eigenvector",
    "clust", "constr", "ecc", "esize", "comsize",
    "udegree", "ucloseness", "ubetweenness", "ueigenvector",
    "uclust", "uconstr", "uecc", "uesize", "ucomsize",
    "edgbet", "udists",
    "mutfr", "dhom", "structsim", "aa", "katz", "srcos",
    "density", "assortativity", "apl", "modularity",
    "netcent", "netvar",
    "chgmod", "chgast", "chapl", "chub", "chuc", "chue",
    "chgcent", "chgvar", "chuclust", "chuconstr", "chuecc", "chuesize",
    "vdist", "vecc", "uvecc",
]
scaler = StandardScaler()
scaler.fit(traindat[numeric_features])
for d in (traindat, valdat, testdat):
    d[numeric_features] = scaler.transform(d[numeric_features])

with open(os.path.join(ARTIFACT_DIR, "scaler.pkl"), "wb") as f:
    pickle.dump(scaler, f)
with open(os.path.join(ARTIFACT_DIR, "scaler_cols.json"), "w") as f:
    json.dump(numeric_features, f)

# ---- feature-group registry -- identical to both notebooks ------------------
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
COMPLETE_GROUPS = list(FEATURE_GROUPS.keys())


def group_columns(active_groups, feature_groups=FEATURE_GROUPS):
    cols = []
    for g in active_groups:
        cols.extend(feature_groups[g])
    return cols


def next_pow2(n):
    return 1 if n <= 1 else 2 ** (n - 1).bit_length()


def coalition_key(subset):
    return "empty" if not subset else "_".join(sorted(subset))


# ---- model builder -- identical to both notebooks ----------------------------
def build_model(active_groups, feature_groups=FEATURE_GROUPS, n_layers=None,
                 width=None, width_rule="pow2", activation="gelu",
                 dropout=0.0, use_batchnorm=True, seed=None, name=None):
    cols = group_columns(active_groups, feature_groups)
    n_features = len(cols)
    if n_layers is None:
        n_layers = len(active_groups)
    if width is None:
        width = next_pow2(n_features) if width_rule == "pow2" else n_features

    inputs = Input(shape=(n_features,), name="features")
    x = inputs
    for i in range(n_layers):
        x = Dense(width, activation=activation, kernel_initializer=GlorotUniform(seed=seed), name=f"dense_{i}")(x)
        if use_batchnorm:
            x = BatchNormalization(name=f"bn_{i}")(x)
    outputs = Dense(1, activation="sigmoid", name="p_choice")(x)
    model = Model(inputs=inputs, outputs=outputs, name=name or "model_" + "_".join(active_groups))
    model.compile(optimizer="adam", loss="binary_crossentropy")
    return model, cols


# ---- fit + evaluate helper -- identical to both notebooks --------------------
TEST_N_OPTIONS = testdat.groupby(["sid", "trial"])["id"].transform("size")
TEST_P_CHANCE = (1.0 / TEST_N_OPTIONS).to_numpy()
TEST_Y_TRUE = testdat[OUTCOME].to_numpy(dtype="float32")
TEST_SID = testdat["sid"].to_numpy()


def brier_skill(y_true, y_pred, p_chance, sid):
    bs = (y_true - y_pred) ** 2
    bs_chance = (y_true - p_chance) ** 2
    df_ = pd.DataFrame({"sid": sid, "bs": bs, "bs_chance": bs_chance})
    g = df_.groupby("sid").mean()
    return (1 - g["bs"] / g["bs_chance"]).rename("bss")


def fit_and_score(active_groups, seed, n_layers=None, epochs=100, batch_size=1024,
                   patience=10, verbose=0, **build_kwargs):
    """Identical to both notebooks' fit_and_score (plot_loss dropped --
    nothing here needs the training-curve plots)."""
    tf.random.set_seed(seed)
    model, cols = build_model(active_groups, n_layers=n_layers, seed=seed, **build_kwargs)

    X_train = traindat[cols].to_numpy(dtype="float32")
    y_train = traindat[OUTCOME].to_numpy(dtype="float32")
    X_val = valdat[cols].to_numpy(dtype="float32")
    y_val = valdat[OUTCOME].to_numpy(dtype="float32")
    X_test = testdat[cols].to_numpy(dtype="float32")

    early_stop = EarlyStopping(monitor="val_loss", patience=patience, restore_best_weights=True)
    model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=epochs,
              batch_size=batch_size, callbacks=[early_stop], verbose=verbose)

    p_test = model.predict(X_test, batch_size=4096, verbose=0).ravel()
    bss = brier_skill(TEST_Y_TRUE, p_test, TEST_P_CHANCE, TEST_SID)
    return model, cols, bss


N_REPEATS = 2
BASE_SEED = 100
N_LAYERS = len(FEATURE_GROUPS)

bss_records = []


def fit_save_repeats(active_groups, save_stem, activation="gelu"):
    """Fits N_REPEATS copies of one variant (same BASE_SEED+r convention as
    both notebooks) and saves each repeat's model + column list, skipping
    any repeat whose files already exist (resumable)."""
    for r in range(N_REPEATS):
        model_path = os.path.join(ARTIFACT_DIR, f"{save_stem}_r{r}.keras")
        cols_path = os.path.join(ARTIFACT_DIR, f"{save_stem}_r{r}_cols.json")
        if os.path.exists(model_path) and os.path.exists(cols_path):
            print(f"  [skip, exists] {save_stem} r{r}")
            continue
        model, cols, bss = fit_and_score(active_groups, seed=BASE_SEED + r, n_layers=N_LAYERS,
                                          activation=activation)
        os.makedirs(os.path.dirname(model_path), exist_ok=True)
        model.save(model_path)
        with open(cols_path, "w") as f:
            json.dump(cols, f)
        bss_records.append({"variant": save_stem, "repeat": r, "mean_bss": bss.mean()})
        print(f"  {save_stem} r{r}: mean participant Brier skill = {bss.mean():.4f}")


print("\n=== complete GELU model (2 repeats) ===")
fit_save_repeats(COMPLETE_GROUPS, "complete_gelu", activation="gelu")

print("\n=== complete linear model (2 repeats) ===")
fit_save_repeats(COMPLETE_GROUPS, "complete_linear", activation="linear")

print("\n=== 63 Shapley coalition models (2 repeats each; visual always included) ===")
full_set = frozenset(SUBSTANTIVE_GROUPS)
n_done = 0
for r_size in range(len(SUBSTANTIVE_GROUPS) + 1):
    for subset in itertools.combinations(SUBSTANTIVE_GROUPS, r_size):
        if frozenset(subset) == full_set:
            continue   # the complete model IS this coalition -- not a separate fit
        active_groups = ["visual"] + list(subset)
        stem = os.path.join("shapley_coalitions", coalition_key(subset))
        fit_save_repeats(active_groups, stem, activation="gelu")
        n_done += 1
        if n_done % 10 == 0:
            print(f"  --- {n_done}/63 coalitions done ---")

pd.DataFrame(bss_records).to_csv(os.path.join(ARTIFACT_DIR, "bss_sanity_check.csv"), index=False)
print(f"\ndone. {n_done} coalitions + complete_gelu + complete_linear saved under {ARTIFACT_DIR}/")
print("see bss_sanity_check.csv to compare against the existing cached_results.")
