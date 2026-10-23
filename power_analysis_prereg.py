"""
power_analysis_prereg.py

Power analysis for RQ4 of the replication preregistration (preregistration_
draft.html): frozen neural-network-model transfer to new data. Covers the
complete model vs. chance, the nonlinear-vs-linear activation advantage, and
each of the six feature groups' Shapley contribution vs. null. All three are
evaluated entirely by bootstrap resampling of the real per-participant
values already on disk (no retraining, matching the preregistration's own
"apply the already-trained model to new data" design) -- this gives exact,
finite-sample power estimates rather than relying on a normal/t
approximation, and is cheap since it only resamples existing arrays.

RQ1-RQ3 (the R-side analyses -- profile replication, sign reversals, PCA
generalization) are covered in power_analysis_prereg.R.
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
# NOTE: RQ4(c) in the preregistration is the MAIN Shapley decomposition
# (zero-order/unique/marginal/Shapley "which feature groups matter" result,
# cached in nn_results_cache.npz), evaluated against null on new data -- not
# the nonlinear-advantage Shapley decomposition from this session's add-on
# analysis. Both are included below since both are plausible readings of
# "each feature group's Shapley contribution"; keep whichever matches your
# final RQ4(c) wording, or both if you want to preregister both.
shapley_bss = pd.DataFrame.from_records(data['shapley_bss'])
if 'sid' in shapley_bss.columns:
    shapley_bss = shapley_bss.set_index('sid')

Ns = [15, 20, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500]
N_REPS = 1000

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

print("RQ4(a): complete model vs. chance (one-sided)")
power_complete = bootstrap_power(complete_bss, Ns, N_REPS, alt='greater')
for n, p in zip(Ns, power_complete):
    print(f"  N={n:5d}  power={p:.3f}")

print("\nRQ4(b): GELU (nonlinear) vs. linear complete model (one-sided)")
power_diff = bootstrap_power(diff_gelu_linear, Ns, N_REPS, alt='greater')
for n, p in zip(Ns, power_diff):
    print(f"  N={n:5d}  power={p:.3f}")

results = {'N': Ns, 'power_complete_vs_chance': power_complete, 'power_gelu_vs_linear': power_diff}

print("\nRQ4(c): each Shapley-decomposed feature group vs. null (two-sided)")
for g in shapley_bss.columns:
    print(f" -- {g} --")
    power_g = bootstrap_power(shapley_bss[g], Ns, N_REPS, alt='two-sided')
    for n, p in zip(Ns, power_g):
        print(f"  N={n:5d}  power={p:.3f}")
    results[f'power_shapley_{g}'] = power_g

out = pd.DataFrame(results)
out.to_csv("power_prereg_rq4_results.csv", index=False)
print("\nwrote power_prereg_rq4_results.csv")
