"""
GSE288358: Genome-wide bivalent enrichment analysis
Placental EWAS (n=151, Illumina EPIC 850K) vs PFAS
Tests whether PFAS-associated differentially methylated probes (DMPs)
are enriched at bivalent domain genes (Court & Arnaud 2017).

Approach:
1. Load beta matrix (151 samples x 365K probes)
2. Load PFAS concentrations + covariates from series matrix
3. For each probe, run linear regression: beta ~ PFAS + covariates
4. Identify significant DMPs (nominal p < 0.05)
5. Map probes to genes (EPIC manifest)
6. Test Fisher exact enrichment of bivalent genes among DMPs
7. Separate tests for HYPER and HYPO methylated probes
"""

import numpy as np
import pickle
import gzip
import os
import sys
from scipy import stats

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# =====================================================
# STEP 1: Load probe-to-gene mapping
# =====================================================
print("=" * 80)
print("GSE288358: GENOME-WIDE BIVALENT ENRICHMENT ANALYSIS")
print("=" * 80)

print("\n[1] Loading probe-to-gene mapping...")
with open('probe_to_gene.pkl', 'rb') as f:
    probe_genes = pickle.load(f)
print(f"  {len(probe_genes)} probes with gene annotations")

# =====================================================
# STEP 2: Load bivalent gene list
# =====================================================
print("\n[2] Loading bivalent domain gene list...")
bivalent_genes = set()
with open(os.path.join('..', 'bivalent_domains_court_arnaud.txt')) as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            bivalent_genes.add(line.upper())
print(f"  {len(bivalent_genes)} bivalent genes loaded")

# =====================================================
# STEP 3: Parse metadata from series matrix
# =====================================================
print("\n[3] Parsing metadata...")
sample_ids = []
pfas_data = {}
covariates = {}

with gzip.open('GSE288358_series_matrix.txt.gz', 'rt', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if line.startswith('!Sample_geo_accession'):
            sample_ids = [p.strip('"') for p in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            values = [p.strip('"') for p in line.split('\t')[1:]]
            if values and ':' in values[0]:
                key = values[0].split(':')[0].strip()

                if key in ['pfhxs', 'pfoa', 'pfna', 'pfos', 'pfda']:
                    vals = []
                    for v in values:
                        try:
                            vals.append(float(v.split(':')[1].strip()))
                        except:
                            vals.append(np.nan)
                    pfas_data[key] = np.array(vals)

                elif key == 'matage_years':
                    vals = [float(v.split(':')[1].strip()) if ':' in v else np.nan for v in values]
                    covariates['matage'] = np.array(vals)

                elif key == 'gestationalage_weeks':
                    vals = [float(v.split(':')[1].strip()) if ':' in v else np.nan for v in values]
                    covariates['gestage'] = np.array(vals)

                elif key == 'Sex':
                    vals = [1.0 if 'Male' in v else 0.0 for v in values]
                    covariates['sex'] = np.array(vals)

                elif key == 'matbmi':
                    vals = [1.0 if '25' in v else 0.0 for v in values]
                    covariates['bmi_cat'] = np.array(vals)

                elif key == 'batch':
                    batch_vals = [v.split(':')[1].strip() for v in values]
                    unique_b = sorted(set(batch_vals))
                    bmap = {b: i for i, b in enumerate(unique_b)}
                    covariates['batch'] = np.array([float(bmap[b]) for b in batch_vals])

                elif key in ['trophoblasts', 'stromal', 'hofbauer',
                             'endothelial', 'nrbc', 'syncytiotrophoblast']:
                    vals = []
                    for v in values:
                        try:
                            vals.append(float(v.split(':')[1].strip()))
                        except:
                            vals.append(np.nan)
                    covariates['cell_' + key] = np.array(vals)

n_samples = len(sample_ids)
print(f"  N={n_samples}, PFAS: {list(pfas_data.keys())}")
print(f"  Covariates: {list(covariates.keys())}")

# Build covariate matrix
cov_names = ['matage', 'gestage', 'sex', 'bmi_cat', 'batch',
             'cell_trophoblasts', 'cell_stromal', 'cell_hofbauer',
             'cell_endothelial', 'cell_nrbc']
cov_matrix = np.column_stack([covariates[c] for c in cov_names])
cov_valid = ~np.any(np.isnan(cov_matrix), axis=1)
print(f"  Samples with complete covariates: {np.sum(cov_valid)}/{n_samples}")

# =====================================================
# STEP 4: Read full beta matrix
# =====================================================
print("\n[4] Reading full beta matrix (this takes a while)...")

with gzip.open('GSE288358_Beta_Matrix.csv.gz', 'rt') as f:
    header = f.readline().strip().split(',')
    probe_ids = [h.strip('"') for h in header[1:]]  # skip row ID column

    n_probes = len(probe_ids)
    print(f"  {n_probes} probes in matrix")

    # Read all data into numpy array
    beta_matrix = np.zeros((n_samples, n_probes), dtype=np.float32)

    for row_idx, line in enumerate(f):
        parts = line.strip().split(',')
        for col_idx in range(n_probes):
            try:
                beta_matrix[row_idx, col_idx] = float(parts[col_idx + 1])
            except:
                beta_matrix[row_idx, col_idx] = np.nan

        if (row_idx + 1) % 50 == 0:
            print(f"    Read {row_idx + 1}/{n_samples} samples...")

    print(f"  Matrix shape: {beta_matrix.shape}")
    print(f"  NaN fraction: {np.isnan(beta_matrix).mean():.4f}")

# =====================================================
# STEP 5: Run EWAS for PFOA (primary analysis)
# =====================================================
print("\n[5] Running genome-wide association: beta ~ PFOA + covariates")

pfoa = pfas_data['pfoa']
mask_pfoa = ~np.isnan(pfoa) & cov_valid
n_valid = np.sum(mask_pfoa)
print(f"  Valid samples for PFOA analysis: {n_valid}")

# Prepare design matrix (once): intercept + PFOA + covariates
X_base = np.column_stack([np.ones(n_valid), pfoa[mask_pfoa], cov_matrix[mask_pfoa]])
k = X_base.shape[1]
df = n_valid - k

# Precompute (X'X)^-1 X' for OLS (same design matrix for all probes)
try:
    XtX_inv = np.linalg.inv(X_base.T @ X_base)
    XtX_inv_Xt = XtX_inv @ X_base.T  # k x n matrix
    hat_diag_se_factor = np.sqrt(np.diag(XtX_inv))  # for SE computation
    print(f"  Design matrix: {X_base.shape}, df={df}")
except np.linalg.LinAlgError:
    print("  ERROR: Singular design matrix!")
    sys.exit(1)

# Vectorized OLS: run all probes at once using matrix multiplication
print("  Computing OLS for all probes (vectorized)...")

Y = beta_matrix[mask_pfoa, :].astype(np.float64)  # n_valid x n_probes

# Beta coefficients: (k x n) @ (n x p) = k x p
betas_all = XtX_inv_Xt @ Y  # k x n_probes

# Residuals
Y_hat = X_base @ betas_all  # n_valid x n_probes
residuals = Y - Y_hat

# MSE for each probe
mse = np.sum(residuals ** 2, axis=0) / df  # n_probes

# SE of PFOA coefficient (index 1)
se_pfoa = hat_diag_se_factor[1] * np.sqrt(mse)  # n_probes

# t-statistics for PFOA
beta_pfoa = betas_all[1, :]  # PFOA coefficients
t_stats = beta_pfoa / se_pfoa

# p-values (two-sided)
pvals = 2 * stats.t.sf(np.abs(t_stats), df)

# Handle probes with all-NaN or zero variance
bad_probes = np.isnan(t_stats) | np.isinf(t_stats) | (se_pfoa == 0)
pvals[bad_probes] = np.nan
t_stats[bad_probes] = np.nan

valid_pvals = ~np.isnan(pvals)
n_tested = np.sum(valid_pvals)
print(f"  Tested {n_tested} probes")

# Also run PFOS
print("\n  Running OLS for PFOS...")
pfos = pfas_data['pfos']
mask_pfos = ~np.isnan(pfos) & cov_valid
n_valid_s = np.sum(mask_pfos)

X_s = np.column_stack([np.ones(n_valid_s), pfos[mask_pfos], cov_matrix[mask_pfos]])
XtX_inv_s = np.linalg.inv(X_s.T @ X_s)
XtX_inv_Xt_s = XtX_inv_s @ X_s.T
hat_se_s = np.sqrt(np.diag(XtX_inv_s))
df_s = n_valid_s - X_s.shape[1]

Y_s = beta_matrix[mask_pfos, :].astype(np.float64)
betas_s = XtX_inv_Xt_s @ Y_s
resid_s = Y_s - X_s @ betas_s
mse_s = np.sum(resid_s ** 2, axis=0) / df_s
se_s = hat_se_s[1] * np.sqrt(mse_s)
beta_pfos = betas_s[1, :]
t_s = beta_pfos / se_s
pvals_pfos = 2 * stats.t.sf(np.abs(t_s), df_s)
bad_s = np.isnan(t_s) | np.isinf(t_s) | (se_s == 0)
pvals_pfos[bad_s] = np.nan

# Also run PFHxS
print("  Running OLS for PFHxS...")
pfhxs = pfas_data['pfhxs']
mask_pfhxs = ~np.isnan(pfhxs) & cov_valid
n_valid_h = np.sum(mask_pfhxs)

X_h = np.column_stack([np.ones(n_valid_h), pfhxs[mask_pfhxs], cov_matrix[mask_pfhxs]])
XtX_inv_h = np.linalg.inv(X_h.T @ X_h)
XtX_inv_Xt_h = XtX_inv_h @ X_h.T
hat_se_h = np.sqrt(np.diag(XtX_inv_h))
df_h = n_valid_h - X_h.shape[1]

Y_h = beta_matrix[mask_pfhxs, :].astype(np.float64)
betas_h = XtX_inv_Xt_h @ Y_h
resid_h = Y_h - X_h @ betas_h
mse_h = np.sum(resid_h ** 2, axis=0) / df_h
se_h = hat_se_h[1] * np.sqrt(mse_h)
beta_pfhxs = betas_h[1, :]
t_h = beta_pfhxs / se_h
pvals_pfhxs = 2 * stats.t.sf(np.abs(t_h), df_h)
bad_h = np.isnan(t_h) | np.isinf(t_h) | (se_h == 0)
pvals_pfhxs[bad_h] = np.nan

# =====================================================
# STEP 6: Compute BH-FDR
# =====================================================
print("\n[6] Computing BH-FDR correction...")

def benjamini_hochberg(pvals):
    n = len(pvals)
    valid = ~np.isnan(pvals)
    result = np.full(n, np.nan)
    pv = pvals[valid]
    m = len(pv)
    if m == 0:
        return result
    idx = np.argsort(pv)
    sorted_p = pv[idx]
    q = sorted_p * m / (np.arange(1, m + 1))
    # Enforce monotonicity
    for i in range(m - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q = np.minimum(q, 1.0)
    # Map back
    orig_idx = np.where(valid)[0]
    for i, si in enumerate(idx):
        result[orig_idx[si]] = q[i]
    return result

qvals_pfoa = benjamini_hochberg(pvals)
qvals_pfos = benjamini_hochberg(pvals_pfos)
qvals_pfhxs = benjamini_hochberg(pvals_pfhxs)

for name, pv, qv in [('PFOA', pvals, qvals_pfoa), ('PFOS', pvals_pfos, qvals_pfos),
                       ('PFHxS', pvals_pfhxs, qvals_pfhxs)]:
    n_nom = np.sum(pv[~np.isnan(pv)] < 0.05)
    n_fdr10 = np.sum(qv[~np.isnan(qv)] < 0.10)
    n_fdr05 = np.sum(qv[~np.isnan(qv)] < 0.05)
    print(f"  {name}: nominal p<0.05 = {n_nom}, FDR<0.10 = {n_fdr10}, FDR<0.05 = {n_fdr05}")

# =====================================================
# STEP 7: Map probes to genes and classify
# =====================================================
print("\n[7] Mapping probes to genes for enrichment analysis...")

# For each PFAS, identify significant probes and map to genes
def get_enrichment_data(pv, bcoef, pfas_name, threshold=0.05):
    """Classify probes as sig/nonsig, hyper/hypo, bivalent/not."""
    sig_mask = pv < threshold
    sig_mask[np.isnan(pv)] = False

    # Map probes to gene-level (a probe is bivalent if ANY of its annotated genes is bivalent)
    # Background: all tested probes with gene annotation
    bg_bivalent = 0
    bg_nonbivalent = 0
    sig_bivalent = 0
    sig_nonbivalent = 0
    sig_hyper_biv = 0
    sig_hyper_nonbiv = 0
    sig_hypo_biv = 0
    sig_hypo_nonbiv = 0

    # Gene-level analysis
    sig_genes_all = set()
    sig_genes_hyper = set()
    sig_genes_hypo = set()
    bg_genes = set()

    for i, pid in enumerate(probe_ids):
        if np.isnan(pv[i]):
            continue

        # Get genes for this probe
        genes = probe_genes.get(pid, set())
        if not genes:
            continue  # Skip probes without gene annotation for gene-level analysis

        genes_upper = set(g.upper() for g in genes)
        is_bivalent = bool(genes_upper & bivalent_genes)

        bg_genes.update(genes_upper)
        if is_bivalent:
            bg_bivalent += 1
        else:
            bg_nonbivalent += 1

        if sig_mask[i]:
            sig_genes_all.update(genes_upper)
            if is_bivalent:
                sig_bivalent += 1
            else:
                sig_nonbivalent += 1

            if bcoef[i] > 0:
                sig_genes_hyper.update(genes_upper)
                if is_bivalent:
                    sig_hyper_biv += 1
                else:
                    sig_hyper_nonbiv += 1
            else:
                sig_genes_hypo.update(genes_upper)
                if is_bivalent:
                    sig_hypo_biv += 1
                else:
                    sig_hypo_nonbiv += 1

    return {
        'pfas': pfas_name,
        'bg_bivalent': bg_bivalent,
        'bg_nonbivalent': bg_nonbivalent,
        'sig_bivalent': sig_bivalent,
        'sig_nonbivalent': sig_nonbivalent,
        'sig_hyper_biv': sig_hyper_biv,
        'sig_hyper_nonbiv': sig_hyper_nonbiv,
        'sig_hypo_biv': sig_hypo_biv,
        'sig_hypo_nonbiv': sig_hypo_nonbiv,
        'sig_genes_all': sig_genes_all,
        'sig_genes_hyper': sig_genes_hyper,
        'sig_genes_hypo': sig_genes_hypo,
        'bg_genes': bg_genes,
    }

res_pfoa = get_enrichment_data(pvals, beta_pfoa, 'PFOA')
res_pfos = get_enrichment_data(pvals_pfos, beta_pfos, 'PFOS')
res_pfhxs = get_enrichment_data(pvals_pfhxs, beta_pfhxs, 'PFHxS')

# =====================================================
# STEP 8: Fisher exact tests
# =====================================================
print("\n[8] Bivalent enrichment (Fisher exact test)")
print("=" * 80)

def fisher_test(a, b, c, d, label):
    """2x2: a=sig+biv, b=sig+nonbiv, c=bg_biv-a, d=bg_nonbiv-b"""
    table = np.array([[a, b], [c - a, d - b]])
    if any(x < 0 for x in table.flatten()):
        print(f"  {label}: INVALID TABLE")
        return None
    odds, p = stats.fisher_exact(table, alternative='greater')
    total_sig = a + b
    total_bg = c + d
    frac_sig = a / total_sig if total_sig > 0 else 0
    frac_bg = c / total_bg if total_bg > 0 else 0
    enrichment = frac_sig / frac_bg if frac_bg > 0 else 0
    print(f"  {label}:")
    print(f"    Sig DMPs: {a} bivalent / {total_sig} total ({frac_sig:.1%})")
    print(f"    Background: {c} bivalent / {total_bg} total ({frac_bg:.1%})")
    print(f"    Enrichment: {enrichment:.2f}x, Fisher p = {p:.4e}")
    return {'label': label, 'enrichment': enrichment, 'p': p,
            'sig_biv': a, 'sig_total': total_sig, 'bg_frac': frac_bg}

results_table = []

for res in [res_pfoa, res_pfos, res_pfhxs]:
    name = res['pfas']
    print(f"\n--- {name} ---")

    r = fisher_test(res['sig_bivalent'], res['sig_nonbivalent'],
                    res['bg_bivalent'], res['bg_nonbivalent'],
                    f"{name} ALL DMPs")
    if r: results_table.append(r)

    r = fisher_test(res['sig_hyper_biv'], res['sig_hyper_nonbiv'],
                    res['bg_bivalent'], res['bg_nonbivalent'],
                    f"{name} HYPER (gain methylation)")
    if r: results_table.append(r)

    r = fisher_test(res['sig_hypo_biv'], res['sig_hypo_nonbiv'],
                    res['bg_bivalent'], res['bg_nonbivalent'],
                    f"{name} HYPO (lose methylation)")
    if r: results_table.append(r)

# =====================================================
# STEP 9: Gene-level enrichment
# =====================================================
print("\n\n[9] Gene-level bivalent enrichment")
print("=" * 80)

for res in [res_pfoa, res_pfos, res_pfhxs]:
    name = res['pfas']
    bg = res['bg_genes']
    sig_all = res['sig_genes_all']
    sig_hyper = res['sig_genes_hyper']
    sig_hypo = res['sig_genes_hypo']

    bg_biv_genes = bg & bivalent_genes
    bg_nonbiv_genes = bg - bivalent_genes

    print(f"\n--- {name} (gene-level) ---")
    print(f"  Background: {len(bg)} genes ({len(bg_biv_genes)} bivalent = {len(bg_biv_genes)/len(bg):.1%})")

    for label, gene_set in [("ALL DMPs", sig_all), ("HYPER", sig_hyper), ("HYPO", sig_hypo)]:
        n_biv = len(gene_set & bivalent_genes)
        n_nonbiv = len(gene_set - bivalent_genes)
        n_total = n_biv + n_nonbiv

        if n_total == 0:
            print(f"  {label}: no genes")
            continue

        frac = n_biv / n_total
        bg_frac = len(bg_biv_genes) / len(bg)
        enrichment = frac / bg_frac if bg_frac > 0 else 0

        # Fisher: gene-level
        table = [[n_biv, n_nonbiv],
                 [len(bg_biv_genes) - n_biv, len(bg_nonbiv_genes) - n_nonbiv]]
        if all(x >= 0 for row in table for x in row):
            odds, p = stats.fisher_exact(table, alternative='greater')
        else:
            p = np.nan

        print(f"  {label}: {n_biv}/{n_total} bivalent ({frac:.1%}), "
              f"enrichment {enrichment:.2f}x, Fisher p = {p:.4e}")

# =====================================================
# STEP 10: HOX gene analysis
# =====================================================
print("\n\n[10] HOX gene analysis")
print("=" * 80)

hox_genes = set()
for g in probe_genes.values():
    for gene in g:
        if gene.upper().startswith('HOX'):
            hox_genes.add(gene.upper())

print(f"  HOX genes in EPIC array: {len(hox_genes)}")

for name, pv, bcoef in [('PFOA', pvals, beta_pfoa), ('PFOS', pvals_pfos, beta_pfos),
                          ('PFHxS', pvals_pfhxs, beta_pfhxs)]:
    print(f"\n  --- {name} ---")
    hox_results = []
    for i, pid in enumerate(probe_ids):
        if np.isnan(pv[i]):
            continue
        genes = probe_genes.get(pid, set())
        genes_upper = set(g.upper() for g in genes)
        hox_hit = genes_upper & hox_genes
        if hox_hit and pv[i] < 0.05:
            direction = "HYPER" if bcoef[i] > 0 else "HYPO"
            hox_results.append((pid, sorted(hox_hit), pv[i], bcoef[i], direction))

    hox_results.sort(key=lambda x: x[2])
    if hox_results:
        for pid, genes, p, b, d in hox_results[:20]:
            print(f"    {pid} -> {','.join(genes)}: p={p:.4e}, beta={b:+.4f} ({d})")
        if len(hox_results) > 20:
            print(f"    ... and {len(hox_results) - 20} more")
        print(f"    Total HOX DMPs: {len(hox_results)}")
    else:
        print(f"    No significant HOX DMPs")

# =====================================================
# STEP 11: Summary
# =====================================================
print("\n\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\nProbe-level enrichment (nominal p < 0.05):")
print(f"{'Analysis':<35s} {'Enrich':>8s} {'Fisher p':>12s} {'Sig biv/total':>15s}")
print("-" * 72)
for r in results_table:
    print(f"{r['label']:<35s} {r['enrichment']:>8.2f}x {r['p']:>12.4e} "
          f"{r['sig_biv']}/{r['sig_total']}")

print("\n\nInterpretation for Bivalency Erosion Hypothesis:")
print("  If bivalent genes enriched among HYPER DMPs -> ")
print("    gains of DNA methylation at bivalent promoters (consistent with model)")
print("  If bivalent genes enriched among HYPO DMPs ->")
print("    losses of DNA methylation at bivalent promoters (less expected)")
print("  If no enrichment -> bivalent genes not preferentially targeted")

# Save results
import json
output = {
    'n_probes': n_probes,
    'n_samples': n_samples,
    'probe_level': [{k: v for k, v in r.items()} for r in results_table],
}
with open('ewas_bivalent_results.json', 'w') as f:
    json.dump(output, f, indent=2)
print("\nResults saved to ewas_bivalent_results.json")

print("\nDONE.")
