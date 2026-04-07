"""
GSE288358 v2: PFAS vs DNA Methylation with FDR correction + covariates
Linear regression: beta ~ PFAS + matage + gestage + sex + bmi + batch + cell_types
Benjamini-Hochberg FDR correction
"""
import urllib.request, json, time, pickle, sys, gzip, os
import numpy as np
from scipy import stats

os.chdir(os.path.dirname(os.path.abspath(__file__)))

TARGET_GENES = [
    'PHF8','KDM7B','KDM2B',
    'HOXA1','HOXA3','HOXA5','HOXA7','HOXA9','HOXA10','HOXA11','HOXA13',
    'EZH2','SUZ12','EED',
    'PRNP',
    'HSPA5','HSPA8','HSPH1','DNAJB1','BAG3',
    'TTR',
    'ATP7A','ATP7B','SLC31A1','MT3',
    'FTH1','FTL','IREB2',
    'TFEB','MAP1LC3A','MAP1LC3B','LAMP1','LAMP2','CTSD','NPC1',
    'DHCR7','DHCR24','HMGCS1','HMGCR',
    'SOD3','DDIT3','XBP1','ERN1',
]

PROMOTER_MARGIN = 1500

# =====================================================
# STEP 1: Get gene coordinates from NCBI (or load cache)
# =====================================================
COORD_CACHE = 'gene_coords_cache.pkl'
if os.path.exists(COORD_CACHE):
    print("Loading cached gene coordinates...")
    with open(COORD_CACHE, 'rb') as f:
        gene_coords = pickle.load(f)
    print("  {} genes loaded from cache".format(len(gene_coords)))
else:
    print("Getting gene coordinates from NCBI...")
    gene_coords = {}
    for gene in TARGET_GENES:
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}[sym]+AND+Homo+sapiens[orgn]&retmode=json".format(gene)
            with urllib.request.urlopen(url, timeout=10) as resp:
                data = json.loads(resp.read())
            ids = data['esearchresult']['idlist']
            if not ids:
                print("  {}: NOT FOUND".format(gene))
                continue
            gene_id = ids[0]
            url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}&retmode=json".format(gene_id)
            with urllib.request.urlopen(url2, timeout=10) as resp:
                data2 = json.loads(resp.read())
            info = data2['result'][gene_id]
            gi_list = info.get('genomicinfo', [])
            if gi_list:
                gi = gi_list[0]
                chrom = info.get('chromosome', '')
                start = int(gi['chrstart']) - PROMOTER_MARGIN
                end = int(gi['chrstop']) + PROMOTER_MARGIN
                if start > end:
                    start, end = end, start
                gene_coords[gene] = (chrom, start, end)
                print("  {}: chr{}:{}-{}".format(gene, chrom, start, end))
            time.sleep(0.35)
        except Exception as e:
            print("  {}: ERROR {}".format(gene, e))
            time.sleep(0.5)
    with open(COORD_CACHE, 'wb') as f:
        pickle.dump(gene_coords, f)
    print("  Cached to {}".format(COORD_CACHE))

print("\nGot coordinates for {}/{} genes".format(len(gene_coords), len(TARGET_GENES)))

# =====================================================
# STEP 2: Map probes to genes
# =====================================================
print("\nLoading probe coordinate index...")
with open('probe_coords.pkl', 'rb') as f:
    probes = pickle.load(f)
print("  {} probes loaded".format(len(probes)))

gene_probes = {g: [] for g in gene_coords}
for probe_id, (pchr, ppos) in probes.items():
    for gene, (gchr, gstart, gend) in gene_coords.items():
        if pchr == gchr and gstart <= ppos <= gend:
            gene_probes[gene].append(probe_id)

total = sum(len(v) for v in gene_probes.values())
print("Mapped {} probes to {} genes".format(total, sum(1 for v in gene_probes.values() if v)))

# =====================================================
# STEP 3: Parse ALL metadata (PFAS + covariates)
# =====================================================
print("\nParsing metadata (PFAS + covariates)...")
sample_ids = []
pfas_data = {}
covariates = {}

with open('GSE288358_series_matrix.txt', encoding='utf-8') as f:
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
                    vals = []
                    for v in values:
                        try:
                            vals.append(float(v.split(':')[1].strip()))
                        except:
                            vals.append(np.nan)
                    covariates['matage'] = np.array(vals)

                elif key == 'gestationalage_weeks':
                    vals = []
                    for v in values:
                        try:
                            vals.append(float(v.split(':')[1].strip()))
                        except:
                            vals.append(np.nan)
                    covariates['gestage'] = np.array(vals)

                elif key == 'Sex':
                    # Male=1, Female=0
                    vals = []
                    for v in values:
                        s = v.split(':')[1].strip()
                        vals.append(1.0 if s == 'Male' else 0.0)
                    covariates['sex'] = np.array(vals)

                elif key == 'matbmi':
                    # category: 19.5-24.9 -> 0 (normal), 25-35 -> 1 (overweight)
                    vals = []
                    for v in values:
                        s = v.split(':')[1].strip()
                        vals.append(1.0 if '25' in s else 0.0)
                    covariates['bmi_cat'] = np.array(vals)

                elif key == 'batch':
                    # Convert batch to numeric categories
                    batch_vals = [v.split(':')[1].strip() for v in values]
                    unique_batches = sorted(set(batch_vals))
                    batch_map = {b: i for i, b in enumerate(unique_batches)}
                    covariates['batch'] = np.array([float(batch_map[b]) for b in batch_vals])

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
print("  N={}, PFAS: {}".format(n_samples, list(pfas_data.keys())))
print("  Covariates: {}".format(list(covariates.keys())))

# Create PFAS sum (z-score standardized)
pfas_names = list(pfas_data.keys())
pfas_z = {}
for pf in pfas_names:
    v = pfas_data[pf].copy()
    mask = ~np.isnan(v)
    v[mask] = (v[mask] - np.nanmean(v)) / np.nanstd(v)
    pfas_z[pf] = v

pfas_sum = np.zeros(n_samples)
valid_count = np.zeros(n_samples)
for pf in pfas_names:
    m = ~np.isnan(pfas_z[pf])
    pfas_sum[m] += pfas_z[pf][m]
    valid_count[m] += 1
pfas_sum = np.where(valid_count > 0, pfas_sum / valid_count, np.nan)
pfas_data['pfas_sum'] = pfas_sum
print("  Added PFAS_sum (mean of z-scores)")

# =====================================================
# STEP 4: Read beta values
# =====================================================
print("\nReading beta matrix for target probes...")

all_targets = set()
for pl in gene_probes.values():
    all_targets.update(pl)
print("  Looking for {} probes...".format(len(all_targets)))

with gzip.open('GSE288358_Beta_Matrix.csv.gz', 'rt') as f:
    header = [h.strip('"') for h in f.readline().strip().split(',')]

    col_map = {}
    for i, h in enumerate(header):
        if h in all_targets:
            col_map[h] = i

    print("  Found {} probes in matrix".format(len(col_map)))

    if not col_map:
        print("  ERROR: No target probes found!")
        sys.exit(1)

    beta = {p: [] for p in col_map}
    for line in f:
        parts = line.strip().split(',')
        for probe, idx in col_map.items():
            try:
                beta[probe].append(float(parts[idx]))
            except:
                beta[probe].append(np.nan)

for p in beta:
    beta[p] = np.array(beta[p])

print("  Read {} samples".format(len(beta[list(beta.keys())[0]])))

# =====================================================
# STEP 5: Build covariate matrix
# =====================================================
print("\nBuilding covariate matrix...")

# Use: matage, gestage, sex, bmi_cat, batch, cell types (drop syncytio — reference)
cov_names = ['matage', 'gestage', 'sex', 'bmi_cat', 'batch',
             'cell_trophoblasts', 'cell_stromal', 'cell_hofbauer',
             'cell_endothelial', 'cell_nrbc']

cov_matrix = np.column_stack([covariates[c] for c in cov_names])
print("  Covariates: {} ({} columns)".format(cov_names, cov_matrix.shape[1]))

# Check for NaNs in covariates
cov_valid = ~np.any(np.isnan(cov_matrix), axis=1)
print("  Samples with complete covariates: {}/{}".format(np.sum(cov_valid), n_samples))

# =====================================================
# STEP 6: Linear regression with covariates + FDR
# =====================================================
print("\n" + "=" * 80)
print("ANALYSIS v2: LINEAR REGRESSION + BENJAMINI-HOCHBERG FDR")
print("Model: beta_methylation ~ PFAS + matage + gestage + sex + bmi + batch + cell_types")
print("=" * 80)

def ols_regression(y, X):
    """OLS regression. Returns beta_pfas, t_stat, p_value for first predictor (PFAS)."""
    # Add intercept
    n = len(y)
    X_full = np.column_stack([np.ones(n), X])

    try:
        # beta = (X'X)^-1 X'y
        XtX = X_full.T @ X_full
        Xty = X_full.T @ y
        betas = np.linalg.solve(XtX, Xty)

        # Residuals
        y_hat = X_full @ betas
        residuals = y - y_hat

        # Degrees of freedom
        k = X_full.shape[1]
        df = n - k
        if df <= 0:
            return np.nan, np.nan, np.nan

        # MSE and SE of coefficients
        mse = np.sum(residuals**2) / df
        cov_beta = mse * np.linalg.inv(XtX)
        se = np.sqrt(np.diag(cov_beta))

        # t-stat for PFAS coefficient (index 1, after intercept)
        t_stat = betas[1] / se[1]
        p_value = 2 * stats.t.sf(np.abs(t_stat), df)

        return betas[1], t_stat, p_value
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan

def benjamini_hochberg(pvals):
    """Compute BH-adjusted p-values (FDR q-values)."""
    n = len(pvals)
    if n == 0:
        return np.array([])

    # Handle NaN
    valid = ~np.isnan(pvals)
    idx_valid = np.where(valid)[0]
    pvals_valid = pvals[valid]

    n_valid = len(pvals_valid)
    if n_valid == 0:
        return np.full(n, np.nan)

    # Sort
    sorted_idx = np.argsort(pvals_valid)
    sorted_pvals = pvals_valid[sorted_idx]

    # BH adjustment
    qvals = np.zeros(n_valid)
    for i in range(n_valid):
        rank = i + 1
        qvals[i] = sorted_pvals[i] * n_valid / rank

    # Enforce monotonicity (from end to start)
    for i in range(n_valid - 2, -1, -1):
        qvals[i] = min(qvals[i], qvals[i + 1])

    # Cap at 1
    qvals = np.minimum(qvals, 1.0)

    # Map back
    result = np.full(n, np.nan)
    for i, si in enumerate(sorted_idx):
        result[idx_valid[si]] = qvals[i]

    return result

all_results = []

for gene, probes_list in sorted(gene_probes.items()):
    for probe in probes_list:
        if probe not in beta:
            continue
        b = beta[probe]

        for pf_name, pf_vals in pfas_data.items():
            # Valid mask: no NaN in beta, PFAS, or covariates
            mask = ~np.isnan(b) & ~np.isnan(pf_vals) & cov_valid
            n_valid = np.sum(mask)
            if n_valid < 50:
                continue

            # Prepare X: [PFAS, covariates]
            X = np.column_stack([pf_vals[mask], cov_matrix[mask]])
            y = b[mask]

            beta_coef, t_stat, pval = ols_regression(y, X)

            if np.isnan(pval):
                continue

            # Also run simple Spearman for comparison
            rho, sp_pval = stats.spearmanr(b[mask], pf_vals[mask])

            all_results.append({
                'gene': gene,
                'probe': probe,
                'pfas': pf_name,
                'beta_coef': beta_coef,  # regression coefficient for PFAS
                't_stat': t_stat,
                'pval_adj': pval,  # p-value from adjusted model
                'pval_unadj': sp_pval,  # unadjusted Spearman p-value
                'rho': rho,  # Spearman rho (for direction)
                'n': n_valid,
                'mean_beta': float(np.mean(b[mask])),
                'direction': 'HYPER' if beta_coef > 0 else 'HYPO'
            })

print("\nTotal tests: {}".format(len(all_results)))

# Apply BH FDR correction
pvals = np.array([r['pval_adj'] for r in all_results])
qvals = benjamini_hochberg(pvals)
for i, r in enumerate(all_results):
    r['qval'] = qvals[i]

# Also FDR on unadjusted
pvals_unadj = np.array([r['pval_unadj'] for r in all_results])
qvals_unadj = benjamini_hochberg(pvals_unadj)
for i, r in enumerate(all_results):
    r['qval_unadj'] = qvals_unadj[i]

# Summary
sig_nominal = [r for r in all_results if r['pval_adj'] < 0.05]
sig_fdr10 = [r for r in all_results if r['qval'] < 0.10]
sig_fdr05 = [r for r in all_results if r['qval'] < 0.05]

print("\n--- ADJUSTED MODEL (with covariates) ---")
print("Nominal p < 0.05: {}".format(len(sig_nominal)))
print("FDR q < 0.10: {}".format(len(sig_fdr10)))
print("FDR q < 0.05: {}".format(len(sig_fdr05)))

sig_nominal_unadj = [r for r in all_results if r['pval_unadj'] < 0.05]
sig_fdr10_unadj = [r for r in all_results if r['qval_unadj'] < 0.10]
sig_fdr05_unadj = [r for r in all_results if r['qval_unadj'] < 0.05]

print("\n--- UNADJUSTED MODEL (Spearman, for comparison with v1) ---")
print("Nominal p < 0.05: {}".format(len(sig_nominal_unadj)))
print("FDR q < 0.10: {}".format(len(sig_fdr10_unadj)))
print("FDR q < 0.05: {}".format(len(sig_fdr05_unadj)))

# =====================================================
# STEP 7: Detailed results
# =====================================================
print("\n" + "=" * 80)
print("RESULTS: ADJUSTED MODEL (nominal p < 0.05)")
print("=" * 80)

sig_nominal.sort(key=lambda x: x['pval_adj'])

header_fmt = "{:<12s} {:<14s} {:<10s} {:>10s} {:>10s} {:>8s} {:>6s} {:>6s}"
print("\n" + header_fmt.format("Gene", "Probe", "PFAS", "p(adj)", "q(FDR)", "Rho", "Dir", "Beta"))
print("-" * 82)

for r in sig_nominal[:80]:
    fdr_flag = ""
    if r['qval'] < 0.05:
        fdr_flag = " ***"
    elif r['qval'] < 0.10:
        fdr_flag = " **"
    elif r['qval'] < 0.20:
        fdr_flag = " *"

    row_fmt = "{:<12s} {:<14s} {:<10s} {:>10.4f} {:>10.4f} {:>+8.3f} {:>6s} {:>6.3f}{}"
    print(row_fmt.format(
        r['gene'], r['probe'], r['pfas'],
        r['pval_adj'], r['qval'], r['rho'],
        r['direction'], r['mean_beta'], fdr_flag
    ))

# =====================================================
# STEP 8: Summary by gene (adjusted model)
# =====================================================
print("\n" + "=" * 80)
print("SUMMARY BY GENE (adjusted model, nominal p < 0.05)")
print("=" * 80)

gene_sig = {}
for r in sig_nominal:
    g = r['gene']
    if g not in gene_sig:
        gene_sig[g] = []
    gene_sig[g].append(r)

for gene in sorted(gene_sig.keys()):
    hits = gene_sig[gene]
    pfas_set = sorted(set(h['pfas'] for h in hits))
    hyper = sum(1 for h in hits if h['direction'] == 'HYPER')
    hypo = len(hits) - hyper
    best = sorted(hits, key=lambda x: x['pval_adj'])[0]
    best_q = best['qval']

    fdr_note = ""
    if best_q < 0.05:
        fdr_note = " [FDR<0.05 !!!]"
    elif best_q < 0.10:
        fdr_note = " [FDR<0.10 !!]"
    elif best_q < 0.20:
        fdr_note = " [FDR<0.20 !]"

    print("\n  {}:{}".format(gene, fdr_note))
    print("    {} sig CpGs, PFAS: {}, {} hyper / {} hypo".format(
        len(hits), pfas_set, hyper, hypo))
    print("    Best: {} ~ {}: coef={:+.5f}, p={:.6f}, q={:.4f} ({})".format(
        best['probe'], best['pfas'], best['beta_coef'], best['pval_adj'],
        best['qval'], best['direction']))

# =====================================================
# STEP 9: Comparison v1 vs v2
# =====================================================
print("\n" + "=" * 80)
print("COMPARISON: UNADJUSTED (v1) vs ADJUSTED (v2)")
print("=" * 80)

# Check which genes survive adjustment
v1_genes = set(r['gene'] for r in all_results if r['pval_unadj'] < 0.05)
v2_genes = set(r['gene'] for r in sig_nominal)
fdr_genes = set(r['gene'] for r in all_results if r['qval'] < 0.10)

print("\nGenes significant in unadjusted (v1): {}".format(sorted(v1_genes)))
print("Genes significant in adjusted (v2):   {}".format(sorted(v2_genes)))
print("Genes surviving FDR<0.10:              {}".format(sorted(fdr_genes)))
print("\nLost after adjustment: {}".format(sorted(v1_genes - v2_genes)))
print("Gained after adjustment: {}".format(sorted(v2_genes - v1_genes)))
print("Stable (both): {}".format(sorted(v1_genes & v2_genes)))

# =====================================================
# STEP 10: PFAS_sum results
# =====================================================
print("\n" + "=" * 80)
print("PFAS_SUM RESULTS (mean z-score of 5 PFAS)")
print("=" * 80)

pfas_sum_results = [r for r in all_results if r['pfas'] == 'pfas_sum' and r['pval_adj'] < 0.05]
pfas_sum_results.sort(key=lambda x: x['pval_adj'])

if pfas_sum_results:
    for r in pfas_sum_results:
        fdr_flag = " [FDR<0.10]" if r['qval'] < 0.10 else ""
        print("  {:<12s} {:<14s} coef={:+.5f} p={:.4f} q={:.4f} {}{}".format(
            r['gene'], r['probe'], r['beta_coef'],
            r['pval_adj'], r['qval'], r['direction'], fdr_flag))
else:
    print("  No significant results for PFAS_sum at nominal p<0.05")

# =====================================================
# STEP 11: Final interpretation
# =====================================================
print("\n" + "=" * 80)
print("INTERPRETATION v2")
print("=" * 80)

hyper_genes_adj = sorted(g for g in gene_sig
    if sum(1 for h in gene_sig[g] if h['direction'] == 'HYPER') >
       sum(1 for h in gene_sig[g] if h['direction'] == 'HYPO'))
hypo_genes_adj = sorted(g for g in gene_sig
    if sum(1 for h in gene_sig[g] if h['direction'] == 'HYPO') >
       sum(1 for h in gene_sig[g] if h['direction'] == 'HYPER'))

print("\nHYPERMETHYLATED (silenced) after covariate adjustment:")
for g in hyper_genes_adj:
    best_q = min(r['qval'] for r in gene_sig[g])
    print("  -> {} (best FDR q={:.4f})".format(g, best_q))

print("\nHYPOMETHYLATED (activated) after covariate adjustment:")
for g in hypo_genes_adj:
    best_q = min(r['qval'] for r in gene_sig[g])
    print("  -> {} (best FDR q={:.4f})".format(g, best_q))

print("\nKEY QUESTION: Which v1 findings SURVIVE covariate adjustment + FDR?")
print("If LAMP1/EZH2/SLC31A1 survive -> result is ROBUST, not confounded")
print("If they don't survive -> effect was driven by batch/age/BMI/cell composition")

# Save results
import json as json_mod
output = {
    'n_tests': len(all_results),
    'n_nominal_adj': len(sig_nominal),
    'n_fdr10': len(sig_fdr10),
    'n_fdr05': len(sig_fdr05),
    'n_nominal_unadj': len(sig_nominal_unadj),
    'genes_v1': sorted(v1_genes),
    'genes_v2': sorted(v2_genes),
    'genes_fdr10': sorted(fdr_genes),
    'top_results': [{k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
                     for k, v in r.items()} for r in sig_nominal[:30]]
}
with open('analysis_v2_results.json', 'w') as f:
    json_mod.dump(output, f, indent=2)
print("\nResults saved to analysis_v2_results.json")

print("\nDONE.")
