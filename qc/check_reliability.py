#!/usr/bin/env python3
"""
Statistical Reliability Assessment for RNA-seq DEG Analysis
GSE254408 - SH-SY5Y neurons treated with PFAS compounds
Control vs PFOA, n=3 per group
CPM normalization + Welch t-test + BH FDR correction
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import gzip
import numpy as np
from scipy import stats
from scipy.stats import nct
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

SEP = "=" * 78

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', 'geo_GSE254408')
DATA_FILE = os.path.join(DATA_DIR, 'GSE254408_all_samples_raw_counts.txt.gz')
MAPPING_FILE = os.path.join(DATA_DIR, 'ensembl_to_symbol.tsv')

CTRL_COLS = [8, 15, 22]
PFOA_COLS = [11, 18, 25]

print(SEP)
print("STATISTICAL RELIABILITY ASSESSMENT - GSE254408 Control vs PFOA")
print(SEP)

# Load ensembl mapping
ens2sym = {}
with open(MAPPING_FILE, "r") as f:
    f.readline()
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            ens2sym[parts[0]] = parts[1]

# Load raw counts
gene_ids = []
ctrl_raw = []
pfoa_raw = []

with gzip.open(DATA_FILE, "rt") as f:
    f.readline()
    f.readline()
    for line in f:
        parts = line.strip().split("\t")
        gene_id = parts[0]
        c_vals = [int(parts[i]) for i in CTRL_COLS]
        p_vals = [int(parts[i]) for i in PFOA_COLS]
        gene_ids.append(gene_id)
        ctrl_raw.append(c_vals)
        pfoa_raw.append(p_vals)

ctrl_raw = np.array(ctrl_raw, dtype=np.float64)
pfoa_raw = np.array(pfoa_raw, dtype=np.float64)
n_total_genes = len(gene_ids)

print(f"\nTotal genes in count matrix: {n_total_genes}")
print(f"Control samples: {ctrl_raw.shape[1]}, PFOA samples: {pfoa_raw.shape[1]}")

def compute_cpm(counts):
    lib_sizes = counts.sum(axis=0)
    return counts / lib_sizes * 1e6

all_counts = np.hstack([ctrl_raw, pfoa_raw])
cpm_all = compute_cpm(all_counts)
cpm_ctrl = cpm_all[:, :3]
cpm_pfoa = cpm_all[:, 3:]

cpm_above1 = (cpm_all >= 1).sum(axis=1)
keep = cpm_above1 >= 3

gene_ids_filt = [gene_ids[i] for i in range(n_total_genes) if keep[i]]
cpm_ctrl_filt = cpm_ctrl[keep]
cpm_pfoa_filt = cpm_pfoa[keep]
n_genes = len(gene_ids_filt)

print(f"Genes after CPM>=1 filter: {n_genes}")

log2_ctrl = np.log2(cpm_ctrl_filt + 1)
log2_pfoa = np.log2(cpm_pfoa_filt + 1)

pvals = np.zeros(n_genes)
tvals = np.zeros(n_genes)
log2fc = np.zeros(n_genes)

for i in range(n_genes):
    c = log2_ctrl[i]
    p = log2_pfoa[i]
    t_stat, p_val = stats.ttest_ind(p, c, equal_var=False)
    tvals[i] = t_stat
    pvals[i] = p_val
    log2fc[i] = np.mean(p) - np.mean(c)

_, fdr_bh, _, _ = multipletests(pvals, method="fdr_bh")

n_nom05 = int(np.sum(pvals < 0.05))
n_nom01 = int(np.sum(pvals < 0.01))
n_fdr05 = int(np.sum(fdr_bh < 0.05))
n_fdr10 = int(np.sum(fdr_bh < 0.10))
n_fdr25 = int(np.sum(fdr_bh < 0.25))

nom_mask = pvals < 0.05
n_up = int(np.sum(nom_mask & (log2fc > 0)))
n_down = int(np.sum(nom_mask & (log2fc < 0)))

print(f"\n--- DEG Counts ---")
print(f"Nominal p < 0.05:  {n_nom05} ({n_up} UP, {n_down} DOWN)")
print(f"Nominal p < 0.01:  {n_nom01}")
print(f"BH FDR < 0.25:     {n_fdr25}")
print(f"BH FDR < 0.10:     {n_fdr10}")
print(f"BH FDR < 0.05:     {n_fdr05}")

# =========================================================================
# 1. FALSE DISCOVERY RATE CHECK
# =========================================================================

print("\n" + SEP)
print("1. FALSE DISCOVERY RATE CHECK - Storey pi0 & P-value Distribution")
print(SEP)

print("\n--- P-value Distribution (deciles) ---")
decile_edges = np.linspace(0, 1, 11)
counts_hist, _ = np.histogram(pvals, bins=decile_edges)
expected_per_bin = n_genes / 10

print("         Bin  Observed  Expected  Ratio  Bar")
print("-" * 65)
for i in range(10):
    lo = decile_edges[i]
    hi = decile_edges[i + 1]
    obs = counts_hist[i]
    ratio = obs / expected_per_bin
    bar = "#" * int(round(ratio * 20))
    print(f"  [{lo:.1f}-{hi:.1f})  {obs:8d}  {expected_per_bin:8.0f}  {ratio:6.2f}  {bar}")

lambdas = np.arange(0.05, 0.96, 0.05)
pi0_lambda = np.zeros(len(lambdas))
for j, lam in enumerate(lambdas):
    pi0_lambda[j] = np.sum(pvals > lam) / (n_genes * (1 - lam))

pi0_est = np.sum(pvals > 0.5) / (n_genes * 0.5)
pi0_est = min(pi0_est, 1.0)

pi0_high = np.mean(pi0_lambda[lambdas >= 0.5])
pi0_high = min(pi0_high, 1.0)

print(f"\n--- Storey pi0 Estimation ---")
print(f"pi0 at lambda=0.5:      {pi0_est:.4f}")
print(f"pi0 (avg, lambda>=0.5): {pi0_high:.4f}")
print(f"  -> Estimated fraction of TRUE nulls: {pi0_high:.1%}")
print(f"  -> Estimated fraction of TRUE alternatives: {1 - pi0_high:.1%}")
print(f"  -> Estimated number of truly DE genes: ~{int(n_genes * (1 - pi0_high))}")

expected_fp = n_genes * pi0_high * 0.05
observed_deg = n_nom05
estimated_tp = max(observed_deg - expected_fp, 0)

print(f"\n--- False Positive Estimate at p<0.05 ---")
print(f"Tested genes:              {n_genes}")
print(f"Expected FP (pi0*m*alpha): {expected_fp:.0f}")
print(f"Observed DEGs:             {observed_deg}")
print(f"Estimated true positives:  {estimated_tp:.0f}")
print(f"Estimated FDR among nominal DEGs: {expected_fp / observed_deg:.1%}")
print(f"Estimated true positive rate:     {estimated_tp / observed_deg:.1%}")

print(f"\n--- pi0(lambda) across thresholds ---")
print("  lambda      pi0")
for j in range(0, len(lambdas), 2):
    print(f"    {lambdas[j]:.2f}    {pi0_lambda[j]:.4f}")

# =========================================================================
# 2. PERMUTATION-BASED FDR
# =========================================================================

print("\n" + SEP)
print("2. PERMUTATION-BASED FDR (1000 label shuffles)")
print(SEP)

np.random.seed(42)
n_perm = 1000
perm_deg_counts = np.zeros(n_perm)

log2_combined = np.hstack([log2_ctrl, log2_pfoa])

for perm_i in range(n_perm):
    idx = np.random.permutation(6)
    pc = log2_combined[:, idx[:3]]
    pp = log2_combined[:, idx[3:]]

    n1 = n2 = 3
    m1 = pc.mean(axis=1)
    m2 = pp.mean(axis=1)
    v1 = pc.var(axis=1, ddof=1)
    v2 = pp.var(axis=1, ddof=1)

    se = np.sqrt(v1 / n1 + v2 / n2)
    se[se == 0] = np.inf
    t_perm = (m2 - m1) / se

    num_df = (v1 / n1 + v2 / n2) ** 2
    den_df = (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
    den_df[den_df == 0] = 1
    df = num_df / den_df
    df = np.maximum(df, 1)

    p_perm = 2 * stats.t.sf(np.abs(t_perm), df)
    perm_deg_counts[perm_i] = np.sum(p_perm < 0.05)

mean_null_degs = perm_deg_counts.mean()
median_null_degs = np.median(perm_deg_counts)
std_null_degs = perm_deg_counts.std()
pct95_null = np.percentile(perm_deg_counts, 95)
pct99_null = np.percentile(perm_deg_counts, 99)

print(f"\nNull distribution of DEGs at p<0.05 (from {n_perm} permutations):")
print(f"  Mean:              {mean_null_degs:.1f}")
print(f"  Median:            {median_null_degs:.1f}")
print(f"  Std Dev:           {std_null_degs:.1f}")
print(f"  95th percentile:   {pct95_null:.1f}")
print(f"  99th percentile:   {pct99_null:.1f}")
print(f"  Min:               {perm_deg_counts.min():.0f}")
print(f"  Max:               {perm_deg_counts.max():.0f}")

print(f"\nComparison to observed:")
print(f"  Observed DEGs at p<0.05:     {observed_deg}")
print(f"  Expected under null (mean):  {mean_null_degs:.0f}")
print(f"  Ratio (observed/expected):   {observed_deg / mean_null_degs:.2f}x")
print(f"  Permutation-based FDR:       {mean_null_degs / observed_deg:.1%}")

n_exceed = int(np.sum(perm_deg_counts >= observed_deg))
print(f"  Null permutations >= observed: {n_exceed}/{n_perm} ({n_exceed/n_perm:.1%})")

print(f"\n--- Permutation DEG count distribution ---")
perm_edges = np.linspace(perm_deg_counts.min(), max(perm_deg_counts.max(), observed_deg), 11)
perm_hist, _ = np.histogram(perm_deg_counts, bins=perm_edges)
for i in range(len(perm_hist)):
    lo = perm_edges[i]
    hi = perm_edges[i + 1]
    bar = "#" * max(1, int(perm_hist[i] / n_perm * 100))
    print(f"  [{lo:6.0f}-{hi:6.0f})  {perm_hist[i]:4d}  {bar}")

# =========================================================================
# 3. REPRODUCIBILITY CHECK
# =========================================================================

print("\n" + SEP)
print("3. REPRODUCIBILITY CHECK - Split-sample validation")
print(SEP)

split_results = []

print("\nFor each split: discover DEGs with 2 vs 2, check direction in held-out pair")
print("                         Split   Disc DEGs   Dir Match    Frac")
print("-" * 70)

for pfoa_out in range(3):
    for ctrl_out in range(3):
        disc_pfoa_idx = [j for j in range(3) if j != pfoa_out]
        disc_ctrl_idx = [j for j in range(3) if j != ctrl_out]

        disc_pfoa = log2_pfoa[:, disc_pfoa_idx]
        disc_ctrl = log2_ctrl[:, disc_ctrl_idx]

        disc_pvals = np.zeros(n_genes)
        disc_fc = np.zeros(n_genes)
        for i in range(n_genes):
            t_s, p_v = stats.ttest_ind(disc_pfoa[i], disc_ctrl[i], equal_var=False)
            disc_pvals[i] = p_v
            disc_fc[i] = np.mean(disc_pfoa[i]) - np.mean(disc_ctrl[i])

        disc_mask = disc_pvals < 0.05
        n_disc = int(disc_mask.sum())

        if n_disc == 0:
            split_results.append((0, 0, np.nan))
            continue

        val_pfoa_val = log2_pfoa[:, pfoa_out]
        val_ctrl_val = log2_ctrl[:, ctrl_out]
        val_fc = val_pfoa_val - val_ctrl_val

        dir_match = int(np.sum(np.sign(disc_fc[disc_mask]) == np.sign(val_fc[disc_mask])))
        frac = dir_match / n_disc

        split_name = f"out: PFOA_{pfoa_out+1}, Ctrl_{ctrl_out+1}"
        print(f"  {split_name:>28s}  {n_disc:10d}  {dir_match:10d}  {frac:6.1%}")
        split_results.append((n_disc, dir_match, frac))

fracs_list = [r[2] for r in split_results if not np.isnan(r[2])]
mean_frac = np.mean(fracs_list)
print(f"\n  Average direction concordance: {mean_frac:.1%}")
print(f"  Expected by chance:           50.0%")
print(f"  Above chance:                 {mean_frac - 0.5:+.1%}")

print(f"\n--- Per-replicate direction consistency (full-data DEGs) ---")
nom_degs_mask = pvals < 0.05
n_degs = int(nom_degs_mask.sum())

consistent_all = 0
consistent_2of3 = 0
for i in np.where(nom_degs_mask)[0]:
    pfoa_vals = log2_pfoa[i]
    ctrl_vals = log2_ctrl[i]
    diffs = []
    for pi_idx in range(3):
        for ci_idx in range(3):
            diffs.append(pfoa_vals[pi_idx] - ctrl_vals[ci_idx])
    diffs = np.array(diffs)
    expected_sign = np.sign(log2fc[i])
    n_agree = np.sum(np.sign(diffs) == expected_sign)
    if n_agree == 9:
        consistent_all += 1
    if n_agree >= 6:
        consistent_2of3 += 1

print(f"  Nominal DEGs: {n_degs}")
print(f"  All 9 pairwise comparisons agree in direction: {consistent_all} ({consistent_all/n_degs:.1%})")
print(f"  >= 6/9 pairwise comparisons agree:             {consistent_2of3} ({consistent_2of3/n_degs:.1%})")

# =========================================================================
# 4. EFFECT SIZE CHECK
# =========================================================================

print("\n" + SEP)
print("4. EFFECT SIZE CHECK - log2FC Distribution & Biological Significance")
print(SEP)

up_fc = log2fc[nom_degs_mask & (log2fc > 0)]
down_fc = log2fc[nom_degs_mask & (log2fc < 0)]

all_fc = log2fc[nom_degs_mask]
abs_fc = np.abs(all_fc)

print(f"\n--- Effect sizes among nominal p<0.05 DEGs ---")
print(f"\n  ALL {n_degs} DEGs:")
print(f"    Mean |log2FC|:   {np.mean(abs_fc):.3f}")
print(f"    Median |log2FC|: {np.median(abs_fc):.3f}")
print(f"    Std |log2FC|:    {np.std(abs_fc):.3f}")

print(f"\n  {len(up_fc)} UP-regulated DEGs:")
print(f"    Mean log2FC:     {np.mean(up_fc):.3f}")
print(f"    Median log2FC:   {np.median(up_fc):.3f}")
print(f"    Range:           [{np.min(up_fc):.3f}, {np.max(up_fc):.3f}]")

print(f"\n  {len(down_fc)} DOWN-regulated DEGs:")
print(f"    Mean log2FC:     {np.mean(down_fc):.3f}")
print(f"    Median log2FC:   {np.median(down_fc):.3f}")
print(f"    Range:           [{np.min(down_fc):.3f}, {np.max(down_fc):.3f}]")

print(f"\n--- |log2FC| distribution among nominal DEGs ---")
fc_bins = [0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, np.inf]
fc_labels = ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.5", "0.5-0.75",
             "0.75-1.0", "1.0-1.5", "1.5-2.0", ">2.0"]
fc_hist, _ = np.histogram(abs_fc, bins=fc_bins)

print("  |log2FC| range   Count    Pct  Cumul  Bar")
print("-" * 65)
cumul = 0
for i in range(len(fc_hist)):
    pct = fc_hist[i] / len(all_fc) * 100
    cumul += pct
    bar = "#" * max(1, int(pct / 2))
    print(f"  {fc_labels[i]:>12s}  {fc_hist[i]:6d}  {pct:5.1f}%  {cumul:5.1f}%  {bar}")

print(f"\n--- Cohen d (standardized effect size) ---")
cohens_d = np.zeros(n_genes)
for i in range(n_genes):
    pooled_std = np.sqrt((np.var(log2_ctrl[i], ddof=1) + np.var(log2_pfoa[i], ddof=1)) / 2)
    if pooled_std > 0:
        cohens_d[i] = (np.mean(log2_pfoa[i]) - np.mean(log2_ctrl[i])) / pooled_std
    else:
        cohens_d[i] = 0

d_degs = np.abs(cohens_d[nom_degs_mask])
print(f"  Among {n_degs} nominal DEGs:")
print(f"    Mean |d|:   {np.mean(d_degs):.2f}")
print(f"    Median |d|: {np.median(d_degs):.2f}")
n_large = int(np.sum(d_degs >= 0.8))
n_vlarge = int(np.sum(d_degs >= 1.2))
n_huge = int(np.sum(d_degs >= 2.0))
print(f"    |d| >= 0.8 (large):   {n_large} ({n_large/n_degs:.1%})")
print(f"    |d| >= 1.2 (v.large): {n_vlarge} ({n_vlarge/n_degs:.1%})")
print(f"    |d| >= 2.0 (huge):    {n_huge} ({n_huge/n_degs:.1%})")

print(f"\n--- Statistical Power at n=3 per group ---")
for d in [0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]:
    n_per = 3
    df_power = 2 * n_per - 2
    ncp_val = d * np.sqrt(n_per / 2)
    t_crit = stats.t.ppf(1 - 0.025, df_power)
    power = 1 - nct.cdf(t_crit, df_power, ncp_val) + nct.cdf(-t_crit, df_power, ncp_val)
    print(f"    d={d:.1f}: power = {power:.1%}")

# =========================================================================
# 5. SUMMARY
# =========================================================================

print("\n" + SEP)
print("5. OVERALL RELIABILITY ASSESSMENT")
print(SEP)

n_fc05 = int(np.sum(abs_fc >= 0.5))

print(f"""
Key Findings:
=============

1. P-VALUE DISTRIBUTION / pi0:
   - Estimated pi0 = {pi0_high:.3f} (fraction of true null hypotheses)
   - ~{(1-pi0_high)*100:.1f}% of tested genes ({int(n_genes*(1-pi0_high))}) may be truly DE
   - Among {observed_deg} nominal DEGs, estimated ~{estimated_tp:.0f} true positives
   - Estimated FDR among nominal DEGs: {expected_fp/observed_deg:.1%}

2. PERMUTATION FDR:
   - Under random label shuffles, we expect {mean_null_degs:.0f} +/- {std_null_degs:.0f} DEGs
   - We observe {observed_deg}, which is {observed_deg/mean_null_degs:.1f}x the null expectation
   - Permutation-based FDR: {mean_null_degs/observed_deg:.1%}
   - {n_exceed}/{n_perm} permutations produced >= {observed_deg} DEGs

3. REPRODUCIBILITY:
   - Average direction concordance in leave-one-out: {mean_frac:.1%} (vs 50% by chance)
   - {consistent_all}/{n_degs} DEGs ({consistent_all/n_degs:.1%}) show fully consistent direction
     across all pairwise replicate comparisons

4. EFFECT SIZES:
   - Median |log2FC| = {np.median(abs_fc):.3f} among nominal DEGs
   - {n_fc05}/{n_degs} ({n_fc05/n_degs:.1%}) have |log2FC| >= 0.5
   - At n=3, even a large effect (d=2.0) has only ~50% power at alpha=0.05
""")

print("BOTTOM LINE:")
if observed_deg / mean_null_degs > 2.0 and mean_frac > 0.6:
    print("  There IS detectable signal above noise, but reliability is LIMITED.")
    perm_fdr_pct = f"{mean_null_degs/observed_deg:.0%}"
    conc_pct = f"{mean_frac:.0%}"
    print(f"  The {observed_deg} nominal DEGs contain an estimated {perm_fdr_pct} false")
    print(f"  positives. Direction concordance of {conc_pct} across splits exceeds chance.")
    print("  Individual gene-level calls should be treated as HYPOTHESES, not conclusions.")
elif observed_deg / mean_null_degs > 1.5:
    perm_fdr_pct = f"{mean_null_degs/observed_deg:.0%}"
    print("  There is MODEST signal above noise.")
    print(f"  The permutation-based FDR of {perm_fdr_pct} means a large fraction")
    print("  of nominal DEGs are likely false positives.")
    print("  Focus on the strongest effects (largest |log2FC|) and pathway-level patterns.")
else:
    print("  The signal is WEAK and barely distinguishable from noise.")
    print(f"  The observed DEG count ({observed_deg}) is close to what we expect under the null")
    print(f"  ({mean_null_degs:.0f}), suggesting most are false positives.")
    print("  Gene-level results are unreliable; aggregate/pathway analysis may still be informative.")

print("""
RECOMMENDATIONS:
  1. Do NOT rely on individual gene FDR-corrected results at n=3
  2. Aggregate analyses (GSEA, pathway enrichment) are more robust
  3. Focus on genes with |log2FC| > 0.5 and consistent replicate behavior
  4. Validation (qPCR, independent dataset) is essential for any specific gene
  5. Report effect sizes alongside p-values; n=3 is severely underpowered
""")

print(SEP)
print("Analysis complete.")
print(SEP)
