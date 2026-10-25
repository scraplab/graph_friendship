"""
power_analysis_prereg.py

Power analysis for RQ3 of the replication preregistration (preregistration_
draft.html): frozen neural-network-model transfer to new data. Covers the
complete model vs. chance, the nonlinear-vs-linear activation advantage, and
each of the six feature groups' Shapley contribution vs. null. All three are
evaluated entirely by bootstrap resampling of the real per-participant
values already on disk (no retraining, matching the preregistration's own
"apply the already-trained model to new data" design) -- this gives exact,
finite-sample power estimates rather than relying on a normal/t
approximation, and is cheap since it only resamples existing arrays.

RQ1-RQ2 (the R-side analyses -- profile replication, PCA generalization)
are covered in power_analysis_prereg.R.

(The preregistration originally also specified a sign-reversal RQ, on
mutual friend count and target distance; it was cut before registration,
and the RQ numbering was shifted down accordingly -- the old RQ4 is this
script's RQ3.)
"""
import numpy as np
import pandas as pd
from scipy import stats

rng = np.random.default_rng(42)

data = np.load(r"nn_results_cache.npz", allow_pickle=True)
complete_bss = pd.DataFrame.from_records(data['complete_bss']).set_index('sid')['value']

lvn = pd.read_csv(r"linear_vs_nonlinear_results.csv", index_col=0)
diff_gelu_linear = lvn['gelu'] - lvn['linear']

shap_nl = pd.read_csv(r"shapley_nonlinear_results.csv", index_col=0)
# NOTE: RQ3(c) in the preregistration is the MAIN Shapley decomposition
# (zero-order/unique/marginal/Shapley "which feature groups matter" result,
# cached in nn_results_cache.npz), evaluated against null on new data -- not
# the nonlinear-advantage Shapley decomposition from this session's add-on
# analysis. Both are included below since both are plausible readings of
# "each feature group's Shapley contribution"; keep whichever matches your
# final RQ3(c) wording, or both if you want to preregister both.
shapley_bss = pd.DataFrame.from_records(data['shapley_bss'])
if 'sid' in shapley_bss.columns:
    shapley_bss = shapley_bss.set_index('sid')

Ns = [15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500]
N_REPS = 1000
TARGET_N = 300   # the replication's actual planned sample size

def bootstrap_power(x, Ns, n_reps, alt='two-sided'):
    x = x.dropna().to_numpy()
    out = []
    for N in Ns:
        sig = 0
        for _ in range(n_reps):
            samp = rng.choice(x, size=min(N, len(x)), replace=(N > len(x)))
            t = stats.ttest_1samp(samp, 0, alternative=alt)
            if t.pvalue < .05:
                sig += 1
        out.append(sig / n_reps)
    return out

print("RQ3(a): complete model vs. chance (one-sided)")
power_complete = bootstrap_power(complete_bss, Ns, N_REPS, alt='greater')
for n, p in zip(Ns, power_complete):
    print(f"  N={n:5d}  power={p:.3f}")

print("\nRQ3(b): GELU (nonlinear) vs. linear complete model (one-sided)")
power_diff = bootstrap_power(diff_gelu_linear, Ns, N_REPS, alt='greater')
for n, p in zip(Ns, power_diff):
    print(f"  N={n:5d}  power={p:.3f}")

results = {'N': Ns, 'power_complete_vs_chance': power_complete, 'power_gelu_vs_linear': power_diff}

print("\nRQ3(c): each Shapley-decomposed feature group vs. null (two-sided)")
for g in shapley_bss.columns:
    print(f" -- {g} --")
    power_g = bootstrap_power(shapley_bss[g], Ns, N_REPS, alt='two-sided')
    for n, p in zip(Ns, power_g):
        print(f"  N={n:5d}  power={p:.3f}")
    results[f'power_shapley_{g}'] = power_g

out = pd.DataFrame(results)
out.to_csv("power_prereg_rq3_results.csv", index=False)
print("\nwrote power_prereg_rq3_results.csv")

print(f"\n=== HEADLINE: power at the planned N = {TARGET_N} ===")
row = out.loc[out['N'] == TARGET_N].iloc[0]
print(f"  RQ3a complete vs. chance:        power = {row['power_complete_vs_chance']:.4f}")
print(f"  RQ3b GELU vs. linear:            power = {row['power_gelu_vs_linear']:.4f}")
for g in shapley_bss.columns:
    print(f"  RQ3c shapley[{g}]: power = {row[f'power_shapley_{g}']:.4f}")
