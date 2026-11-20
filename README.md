# Intrinsic fluctuations constrain mechanical information in active cytoskeletal networks

Source data and analysis code for:

> X. Xiao, S. Chen, V. Yadav, I. Linsmeier, M. P. Murrell, "Intrinsic fluctuations constrain mechanical information in active cytoskeletal networks." *Proc. Natl. Acad. Sci. U.S.A.* (2026). https://doi.org/10.1073/pnas.2609991123

Each `data/` file contains the values plotted in the corresponding figure panel. Raw microscopy images and PIV fields are not included.

## Data (`data/`)

| Figure | File | Contents |
|---|---|---|
| Fig. 1D | `Fig1/Fig1D_Rrp_vs_RL.csv` | PIV-defined response radius R_rp vs. laser radius R_L (n = 18) |
| Fig. 1E | `Fig1/Fig1E_cumulative_strain.csv` | Cumulative strain (ROI mean ± SD) vs. time; selected response states S1, S1′, S1″ (d′ ≥ 2) |
| Fig. 1F | `Fig1/Fig1F_STD_50s_100s.csv` | Temporal STD of actin intensity, 50 s vs. 100 s activation (n = 18 per group) |
| Fig. 1G | `Fig1/Fig1G_Rrp_crosslinker.csv` | Response radius for β = 0 (n = 11) and β = 1/50 (n = 9) |
| Fig. 1H | `Fig1/Fig1H_D2max_crosslinker.csv` | Maximum ROI-averaged non-affinity D²_max for β = 0 (n = 11) and β = 1/50 (n = 9) |
| Fig. 2E | `Fig2/Fig2E_matrixA_normalizedSTD.csv` | Matrix A: normalized temporal STD at the 8 × 4 lattice sites |
| Fig. 2G–H | `Fig2/Fig2G-H_KL_all_8191_configurations.csv` | KL divergence of all 2¹³ − 1 = 8,191 candidate configurations, ranked (rank 1 = decoded "YALE", D_KL = 0.487; rank 2 = 0.511; decoding margin 0.024) |
| Fig. 3G | `Fig3/Fig3G_line_profiles.csv`, `Fig3/Fig3G_bundling_parameter_summary.csv` | Row-averaged actin line profiles at t = 10 s and 180 s; bundling parameter λ (0.058, 0.131) |
| Fig. 3I | `Fig3/Fig3I_KL_vs_lambda_points.csv`, `Fig3/Fig3I_MI_vs_lambda_points.csv` | D_KL and I(X₁; X₂) vs. λ for each dual-spot experiment (column `IncludedInFit` marks points used in the regression) |
| Fig. 3I | `Fig3/Fig3I_KL_individual_fits.csv`, `Fig3/Fig3I_MI_individual_fits.csv`, `Fig3/Fig3I_representative_experiment_summary.csv` | Per-experiment linear fits (slope, R², P) and the representative experiment (exp. 10) |

Notes:
- Fig. 2G: site numbering is row-major on the 8 × 4 lattice (site = 4 × (row − 1) + column).
- Fig. 3: experiment 9 was excluded because its network ruptured during acquisition.
- Units: lengths in µm, time in s, KL divergence and mutual information in nats.

## Code (`code/`)

MATLAB scripts unless noted. Replace `PATH_TO_DATA` with the local data directory before running.

| Figure | Script | Analysis |
|---|---|---|
| Fig. 1D | `Fig1/Fig1D_response_radius_linear_fit.m` | Linear fit of R_rp vs. R_L; intercept R* |
| Fig. 1E | `Fig1/Fig1E_cumulative_strain_states.m` | Cumulative strain and state selection (d′ ≥ 2) |
| Fig. 1F | `Fig1/Fig1F_dprime.m` | d′ and box plot |
| Fig. 1F | `Fig1/Fig1F_dprime_bootstrap_CI.py` (Python) | d′ with 95% BCa bootstrap CI (independent resampling of groups, 10,000 resamples) |
| Fig. 1G | `Fig1/Fig1G_crosslinker_response_radius.m` | Response radius vs. crosslinker |
| Fig. 1H | `Fig1/Fig1H_nonaffinity_D2min.m` | Falk–Langer non-affinity D²_min (Eqs. 10–13) |
| Main text | `Fig1/information_density.m` | Geometric information-density estimate (Eqs. 8–9) |
| Fig. 2D | `Fig2/Fig2D_candidate_sites_inward_flux.m` | Candidate active sites from net inward PIV flux |
| Fig. 2E–H | `Fig2/Fig2E-H_KL_decoding.m` | Normalized STD matrix, enumeration of 8,191 configurations, KL decoding (Eqs. 5, 19–21) |
| SI | `Fig2/decoding_noise_robustness.m` | Decoding with added synthetic noise (Eq. 22) |
| Fig. 3G | `Fig3/Fig3G_line_scan_bundling_parameter.m` | Line scans and bundling parameter λ |
| Fig. 3 | `Fig3/bundling_parameter_lambda.m` | λ = σ_I / μ_I for all frames (Eq. 14) |
| Fig. 3 | `Fig3/KL_divergence_single_vs_dual.m` | KL divergence between single- and dual-spot distributions (Eqs. 15–17) |
| Fig. 3 | `Fig3/mutual_information.m` | Spatial mutual information I(X₁; X₂) (Eq. 18) |
| Fig. 3I | `Fig3/Fig3I_KL_MI_vs_lambda_fits.m` | Regression of D_KL and I(X₁; X₂) on λ, per experiment |

## Contact

Corresponding author: Michael P. Murrell (michael.murrell@yale.edu)

## License

MIT License (see `LICENSE`).
