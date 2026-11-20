"""
Fig. 1F: discriminability d' between 50 s and 100 s activations (STD of actin intensity),
with a 95% bias-corrected and accelerated (BCa) bootstrap confidence interval.

    d' = (mean_100 - mean_50) / sqrt((var_50 + var_100) / 2)

The two groups are resampled independently (10,000 resamples).
Published result: d' = 1.54, 95% CI [1.07, 1.97], n = 18 per group.

Usage:  python fig1F_dprime_bootstrap_CI.py [path/to/Fig1F_STD_50s_100s.csv]
The CSV has two columns: STD_50s, STD_100s (one experiment per row).
"""
import sys
import numpy as np
from scipy import stats

path = sys.argv[1] if len(sys.argv) > 1 else "../../data/Fig1/Fig1F_STD_50s_100s.csv"
data = np.genfromtxt(path, delimiter=",", names=True)
x50, x100 = data["STD_50s"], data["STD_100s"]


def dprime(a, b, axis=-1):
    return (b.mean(axis) - a.mean(axis)) / np.sqrt(
        (a.var(axis, ddof=1) + b.var(axis, ddof=1)) / 2
    )


d = dprime(x50, x100)
res = stats.bootstrap(
    (x50, x100),
    dprime,
    n_resamples=10_000,
    method="BCa",
    confidence_level=0.95,
    vectorized=True,
    random_state=70,  # fixed seed for reproducibility
)
lo, hi = res.confidence_interval
print(f"n = {len(x50)} / {len(x100)}")
print(f"d' = {d:.2f}  (95% BCa bootstrap CI [{lo:.2f}, {hi:.2f}])")
