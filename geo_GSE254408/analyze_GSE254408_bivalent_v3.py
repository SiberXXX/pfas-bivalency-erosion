import sys, io, gzip, os
import numpy as np
from scipy import stats
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# =====================================================
# GSE254408 v3: Raw counts (all 21 samples), CPM
# Focus on: individual PFAS, combined, directional enrichment
# =====================================================

os.chdir(os.path.dirname(os.path.abspath(__file__)))

ensembl_to_sym = {}
with open('ensembl_to_symbol.tsv', 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            ensembl_to_sym[parts[0]] = parts[1]

# Load raw counts
gene_ids = []
all_counts = []
with gzip.open('GSE254408_all_samples_raw_counts.txt.gz', 'rt') as f:
    f.readline()  # comment
    header = f.readline().strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        gene_ids.append(parts[0])
        all_counts.append([int(parts[i]) for i in range(6, 27)])

count_matrix = np.array(all_counts, dtype=np.float64)
n_genes, n_samples = count_matrix.shape

# 0-based indices into count_matrix (columns 6..26 → 0..20)
groups = {
    'Control': [2, 9, 16],   # Control_1(col8→2), Control_2(col15→9), Control_3(col22→16)
    'PFOA': [5, 12, 19],
    'PFOS': [6, 13, 20],
    'PFDA': [3, 10, 17],
    'PFDS': [4, 11, 18],
    'FTOH': [0, 7, 14],
    'FTS': [1, 8, 15],
}

gene_symbols = [ensembl_to_sym.get(gid, '') for gid in gene_ids]

# CPM
lib_sizes = count_matrix.sum(axis=0)
cpm = count_matrix / lib_sizes * 1e6

print(f'Genes: {n_genes}, Samples: {n_samples}')
print(f'Mapped to symbols: {sum(1 for s in gene_symbols if s)}')
print(f'Library sizes: {lib_sizes.min():.0f} - {lib_sizes.max():.0f}')

# Load bivalent
bivalent = set()
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bivalent_domains_court_arnaud.txt'), 'r') as f:
    in_list = False
    for line in f:
        line = line.strip()
        if '===== UNIQUE GENE LIST' in line:
            in_list = True
            continue
        if in_list and line and not line.startswith('#') and not line.startswith('='):
            if '\t' in line:
                break
            bivalent.add(line.upper())

print(f'Bivalent genes: {len(bivalent)}')

from statsmodels.stats.multitest import multipletests

def run_deg(grp_idx, ctrl_idx):
    ctrl_cpm = cpm[:, ctrl_idx]
    grp_cpm = cpm[:, grp_idx]
    ctrl_mean = ctrl_cpm.mean(axis=1)
    grp_mean = grp_cpm.mean(axis=1)
    expressed = (ctrl_mean > 1) | (grp_mean > 1)

    results = []
    for i in range(n_genes):
        if not expressed[i] or not gene_symbols[i]:
            continue
        c = ctrl_cpm[i, :]
        g = grp_cpm[i, :]
        if c.std() == 0 and g.std() == 0:
            continue
        t, p = stats.ttest_ind(g, c, equal_var=False)
        if np.isnan(p):
            continue
        fc = (grp_mean[i] + 0.1) / (ctrl_mean[i] + 0.1)
        log2fc = np.log2(fc) if fc > 0 else 0
        results.append((gene_symbols[i], log2fc, p))

    results.sort(key=lambda x: x[2])
    pvals = [r[2] for r in results]
    if not pvals:
        return results, pvals, np.array([])
    _, padj, _, _ = multipletests(pvals, method='fdr_bh')
    return results, pvals, padj

def fisher_test(deg_set, bg_set, biv_set, alt='greater'):
    deg_biv = deg_set & biv_set
    a = len(deg_biv)
    b = len(deg_set) - a
    c = len(bg_set & biv_set) - a
    d = len(bg_set - biv_set) - b
    if d < 0: d = 0  # edge case
    ort, fp = stats.fisher_exact([[a, b], [c, d]], alternative=alt)
    return a, len(deg_set), ort, fp

ctrl_idx = groups['Control']

# =====================================================
# INDIVIDUAL PFAS
# =====================================================
for pfas in ['PFOA', 'PFOS', 'PFDA', 'PFDS', 'FTOH', 'FTS']:
    print(f'\n{"=" * 60}')
    print(f'{pfas} vs Control (raw counts, CPM, n=3 vs n=3)')
    print(f'{"=" * 60}')

    results, pvals, padj = run_deg(groups[pfas], ctrl_idx)
    n_tested = len(results)
    n_fdr05 = sum(1 for p in padj if p < 0.05) if len(padj) > 0 else 0
    n_nom = sum(1 for p in pvals if p < 0.05)
    print(f'Tested: {n_tested}, FDR<0.05: {n_fdr05}, nom p<0.05: {n_nom}')

    bg = set(r[0].upper() for r in results)
    bg_biv_n = len(bg & bivalent)
    bg_ratio = bg_biv_n / len(bg)
    print(f'Background: {len(bg)} tested, {bg_biv_n} bivalent ({100*bg_ratio:.1f}%)')

    # Nominal DEGs
    nom_all = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05)
    nom_up = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05 and r[1] > 0)
    nom_dn = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05 and r[1] < 0)

    # Enrichment
    a_all, n_all, or_all, fp_all = fisher_test(nom_all, bg, bivalent)
    a_up, n_up, or_up, fp_up = fisher_test(nom_up, bg, bivalent)
    a_dn, n_dn, or_dn, fp_dn = fisher_test(nom_dn, bg, bivalent)

    sig = lambda p: '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
    print(f'  ALL:  {n_all:5d} DEGs, {a_all:3d} biv ({100*a_all/max(n_all,1):.1f}%)  OR={or_all:.2f} p={fp_all:.4e} {sig(fp_all)}')
    print(f'  UP:   {n_up:5d} DEGs, {a_up:3d} biv ({100*a_up/max(n_up,1):.1f}%)  OR={or_up:.2f} p={fp_up:.4e} {sig(fp_up)}')
    print(f'  DOWN: {n_dn:5d} DEGs, {a_dn:3d} biv ({100*a_dn/max(n_dn,1):.1f}%)  OR={or_dn:.2f} p={fp_dn:.4e} {sig(fp_dn)}')

# =====================================================
# COMBINED: All PFAS vs Control (n=18 vs n=3)
# =====================================================
print(f'\n{"=" * 60}')
print('COMBINED: All 18 PFAS vs 3 Control')
print(f'{"=" * 60}')

all_pfas_idx = []
for pfas in ['PFOA', 'PFOS', 'PFDA', 'PFDS', 'FTOH', 'FTS']:
    all_pfas_idx.extend(groups[pfas])

results, pvals, padj = run_deg(all_pfas_idx, ctrl_idx)
n_tested = len(results)
n_fdr05 = sum(1 for p in padj if p < 0.05)
n_fdr10 = sum(1 for p in padj if p < 0.10)
n_nom = sum(1 for p in pvals if p < 0.05)
print(f'Tested: {n_tested}, FDR<0.05: {n_fdr05}, FDR<0.10: {n_fdr10}, nom p<0.05: {n_nom}')

bg = set(r[0].upper() for r in results)
bg_biv_n = len(bg & bivalent)
bg_ratio = bg_biv_n / len(bg)
print(f'Background: {len(bg)} tested, {bg_biv_n} bivalent ({100*bg_ratio:.1f}%)')

# FDR<0.05
fdr05_all = set(r[0].upper() for r, pa in zip(results, padj) if pa < 0.05)
fdr05_up = set(r[0].upper() for r, pa in zip(results, padj) if pa < 0.05 and r[1] > 0)
fdr05_dn = set(r[0].upper() for r, pa in zip(results, padj) if pa < 0.05 and r[1] < 0)

a_all, n_all, or_all, fp_all = fisher_test(fdr05_all, bg, bivalent)
a_up, n_up, or_up, fp_up = fisher_test(fdr05_up, bg, bivalent)
a_dn, n_dn, or_dn, fp_dn = fisher_test(fdr05_dn, bg, bivalent)

sig = lambda p: '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
print(f'\nFDR<0.05:')
print(f'  ALL:  {n_all:5d} DEGs, {a_all:3d} biv ({100*a_all/max(n_all,1):.1f}%)  OR={or_all:.2f} p={fp_all:.4e} {sig(fp_all)}')
print(f'  UP:   {n_up:5d} DEGs, {a_up:3d} biv ({100*a_up/max(n_up,1):.1f}%)  OR={or_up:.2f} p={fp_up:.4e} {sig(fp_up)}')
print(f'  DOWN: {n_dn:5d} DEGs, {a_dn:3d} biv ({100*a_dn/max(n_dn,1):.1f}%)  OR={or_dn:.2f} p={fp_dn:.4e} {sig(fp_dn)}')

# Nominal p<0.05
nom_all = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05)
nom_up = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05 and r[1] > 0)
nom_dn = set(r[0].upper() for r, p in zip(results, pvals) if p < 0.05 and r[1] < 0)

a_all, n_all, or_all, fp_all = fisher_test(nom_all, bg, bivalent)
a_up, n_up, or_up, fp_up = fisher_test(nom_up, bg, bivalent)
a_dn, n_dn, or_dn, fp_dn = fisher_test(nom_dn, bg, bivalent)

print(f'\nNominal p<0.05:')
print(f'  ALL:  {n_all:5d} DEGs, {a_all:3d} biv ({100*a_all/max(n_all,1):.1f}%)  OR={or_all:.2f} p={fp_all:.4e} {sig(fp_all)}')
print(f'  UP:   {n_up:5d} DEGs, {a_up:3d} biv ({100*a_up/max(n_up,1):.1f}%)  OR={or_up:.2f} p={fp_up:.4e} {sig(fp_up)}')
print(f'  DOWN: {n_dn:5d} DEGs, {a_dn:3d} biv ({100*a_dn/max(n_dn,1):.1f}%)  OR={or_dn:.2f} p={fp_dn:.4e} {sig(fp_dn)}')

# =====================================================
# KEY BIVALENT DEGs
# =====================================================
print(f'\n{"=" * 60}')
print('KEY BIVALENT DEGs (Combined, FDR<0.05)')
print(f'{"=" * 60}')

fdr05_list = [(r[0], r[1], pvals[i], padj[i]) for i, r in enumerate(results) if padj[i] < 0.05]
biv_fdr05 = [(g, lfc, p, pa) for g, lfc, p, pa in fdr05_list if g.upper() in bivalent]
print(f'Total FDR<0.05: {len(fdr05_list)}, Bivalent: {len(biv_fdr05)}')
print(f'\nTop 30 bivalent DEGs by p-value:')
for g, lfc, p, pa in sorted(biv_fdr05, key=lambda x: x[2])[:30]:
    dr = 'UP' if lfc > 0 else 'DN'
    print(f'  {g:15s} log2FC={lfc:+.3f} ({dr}) padj={pa:.4e}')

# HOX
print(f'\nHOX genes (any threshold):')
hox_all = [(r[0], r[1], pvals[i], padj[i]) for i, r in enumerate(results)
           if r[0].upper().startswith('HOX') and pvals[i] < 0.05]
if hox_all:
    for g, lfc, p, pa in sorted(hox_all, key=lambda x: x[2]):
        tag = ' [BIV]' if g.upper() in bivalent else ''
        fdr_tag = ' FDR<0.05' if pa < 0.05 else ''
        print(f'  {g:12s} log2FC={lfc:+.3f} p={p:.4e} padj={pa:.4e}{tag}{fdr_tag}')
else:
    print('  None at p<0.05')

# Epigenetic machinery
epi_genes = {'EZH2', 'EZH1', 'SUZ12', 'EED', 'RING1', 'RNF2', 'BMI1',
             'KDM2A', 'KDM2B', 'KDM3A', 'KDM3B', 'KDM4A', 'KDM4B', 'KDM4C',
             'KDM5A', 'KDM5B', 'KDM5C', 'KDM5D', 'KDM6A', 'KDM6B', 'KDM7A',
             'PHF2', 'PHF8', 'JMJD1C', 'JMJD6',
             'SETD1A', 'SETD1B', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D',
             'DNMT1', 'DNMT3A', 'DNMT3B', 'TET1', 'TET2', 'TET3',
             'HDAC1', 'HDAC2', 'HDAC3', 'HDAC4', 'HDAC5', 'HDAC6'}

print(f'\nEpigenetic machinery (combined, all with p<0.10):')
epi_hits = [(r[0], r[1], pvals[i], padj[i]) for i, r in enumerate(results)
            if r[0].upper() in epi_genes and pvals[i] < 0.10]
for g, lfc, p, pa in sorted(epi_hits, key=lambda x: x[2]):
    dr = 'UP' if lfc > 0 else 'DN'
    fdr_tag = 'FDR<0.05' if pa < 0.05 else 'FDR<0.10' if pa < 0.10 else 'nom'
    print(f'  {g:12s} log2FC={lfc:+.3f} ({dr}) p={p:.4e} padj={pa:.4e} [{fdr_tag}]')

# =====================================================
# SUMMARY TABLE
# =====================================================
print(f'\n{"=" * 60}')
print('SUMMARY: Bivalent Enrichment Across All Analyses')
print(f'{"=" * 60}')
print(f'{"Dataset":<30s} {"DEGs":<8s} {"Biv%":<8s} {"BG%":<8s} {"Enrich":<8s} {"Fisher p":<12s}')
print(f'{"-"*30} {"-"*8} {"-"*8} {"-"*8} {"-"*8} {"-"*12}')

for pfas in ['PFOA', 'PFOS', 'PFDA', 'PFDS', 'FTOH', 'FTS']:
    r, pv, pa = run_deg(groups[pfas], ctrl_idx)
    bg2 = set(x[0].upper() for x in r)
    bg2_r = len(bg2 & bivalent) / len(bg2)
    dn = set(x[0].upper() for x, p in zip(r, pv) if p < 0.05 and x[1] < 0)
    dn_b = dn & bivalent
    dn_r = len(dn_b) / max(len(dn), 1)
    _, fp = fisher_test(dn, bg2, bivalent)[-2:]
    sig_s = sig(fp)
    print(f'{pfas+" DOWN nom":<30s} {len(dn):<8d} {100*dn_r:<8.1f} {100*bg2_r:<8.1f} {dn_r/bg2_r if bg2_r>0 else 0:<8.2f} {fp:<12.4e} {sig_s}')

# Combined
r, pv, pa = run_deg(all_pfas_idx, ctrl_idx)
bg2 = set(x[0].upper() for x in r)
bg2_r = len(bg2 & bivalent) / len(bg2)
fdr05_dn2 = set(x[0].upper() for x, p in zip(r, pa) if p < 0.05 and x[1] < 0)
fdr05_dn_b = fdr05_dn2 & bivalent
fdr05_dn_r = len(fdr05_dn_b) / max(len(fdr05_dn2), 1)
_, fp = fisher_test(fdr05_dn2, bg2, bivalent)[-2:]
print(f'{"Combined DN FDR<0.05":<30s} {len(fdr05_dn2):<8d} {100*fdr05_dn_r:<8.1f} {100*bg2_r:<8.1f} {fdr05_dn_r/bg2_r if bg2_r>0 else 0:<8.2f} {fp:<12.4e} {sig(fp)}')

print('\nDone.')
