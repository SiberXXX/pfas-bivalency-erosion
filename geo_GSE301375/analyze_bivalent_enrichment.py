import sys, io, gzip, os
import numpy as np
from scipy import stats
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# =====================================================
# FULL DEG ANALYSIS: GSE301375 (Mouse cerebellum PFHxA)
# Then intersect DEGs with bivalent domains
# =====================================================

os.chdir(os.path.dirname(os.path.abspath(__file__)))

metadata = {
    'GSM9081768_201': ('High', 'M'), 'GSM9081769_205': ('High', 'F'),
    'GSM9081776_241': ('High', 'M'), 'GSM9081777_245': ('High', 'F'),
    'GSM9081780_261': ('High', 'M'), 'GSM9081781_266': ('High', 'F'),
    'GSM9081788_301': ('High', 'M'), 'GSM9081789_305': ('High', 'F'),
    'GSM9081796_351': ('High', 'M'), 'GSM9081797_354': ('High', 'F'),
    'GSM9081798_361': ('High', 'M'), 'GSM9081799_364': ('High', 'F'),
    'GSM9081772_221': ('Control', 'M'), 'GSM9081773_224': ('Control', 'F'),
    'GSM9081774_231': ('Control', 'M'), 'GSM9081775_234': ('Control', 'F'),
    'GSM9081784_281': ('Control', 'M'), 'GSM9081785_283': ('Control', 'F'),
    'GSM9081790_311': ('Control', 'M'), 'GSM9081791_314': ('Control', 'F'),
    'GSM9081792_331': ('Control', 'M'), 'GSM9081793_333': ('Control', 'F'),
    'GSM9081770_211': ('Low', 'M'), 'GSM9081771_213': ('Low', 'F'),
    'GSM9081778_251': ('Low', 'M'), 'GSM9081779_255': ('Low', 'F'),
    'GSM9081782_271': ('Low', 'M'), 'GSM9081783_274': ('Low', 'F'),
    'GSM9081786_291': ('Low', 'M'), 'GSM9081787_297': ('Low', 'F'),
    'GSM9081794_341': ('Low', 'M'), 'GSM9081795_344': ('Low', 'F'),
    'GSM9081800_371': ('Low', 'M'), 'GSM9081801_374': ('Low', 'F'),
}

counts_dir = 'GSE301375_counts'
all_genes = None
data = {}

for fname in sorted(os.listdir(counts_dir)):
    if fname.endswith('_countsU.txt.gz'):
        key = fname.replace('_countsU.txt.gz', '')
        if key not in metadata:
            continue
        genes = []
        values = []
        with gzip.open(os.path.join(counts_dir, fname), 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    genes.append(parts[0])
                    values.append(int(parts[1]))
        if all_genes is None:
            all_genes = genes
        data[key] = np.array(values)

gene_names = all_genes
n_genes = len(gene_names)
ctrl_keys = [k for k, v in metadata.items() if v[0] == 'Control']
high_keys = [k for k, v in metadata.items() if v[0] == 'High']

ctrl_matrix = np.column_stack([data[k] for k in ctrl_keys if k in data])
high_matrix = np.column_stack([data[k] for k in high_keys if k in data])

def to_cpm(mat):
    lib_sizes = mat.sum(axis=0)
    return mat / lib_sizes * 1e6

ctrl_cpm = to_cpm(ctrl_matrix)
high_cpm = to_cpm(high_matrix)
ctrl_mean = ctrl_cpm.mean(axis=1)
high_mean = high_cpm.mean(axis=1)

expressed = (ctrl_mean > 1) | (high_mean > 1)
print(f'Expressed genes (CPM>1): {expressed.sum()} / {n_genes}')

results = []
for i in range(n_genes):
    if not expressed[i]:
        continue
    c = ctrl_cpm[i, :]
    h = high_cpm[i, :]
    if c.std() == 0 and h.std() == 0:
        continue
    t, p = stats.ttest_ind(h, c, equal_var=False)
    fc = (high_mean[i] + 0.1) / (ctrl_mean[i] + 0.1)
    log2fc = np.log2(fc)
    results.append((gene_names[i], log2fc, p, ctrl_mean[i], high_mean[i]))

results.sort(key=lambda x: x[2])
print(f'Tested genes: {len(results)}')

from statsmodels.stats.multitest import multipletests
pvals = [r[2] for r in results]
reject, padj, _, _ = multipletests(pvals, method='fdr_bh')

degs = [(results[i][0], results[i][1], pvals[i], padj[i]) for i in range(len(results)) if padj[i] < 0.05]
print(f'DEGs (FDR<0.05): {len(degs)}')
degs_up = [d for d in degs if d[1] > 0]
degs_down = [d for d in degs if d[1] < 0]
print(f'  Up: {len(degs_up)}, Down: {len(degs_down)}')

degs_01 = [(results[i][0], results[i][1], pvals[i], padj[i]) for i in range(len(results)) if padj[i] < 0.1]
print(f'DEGs (FDR<0.10): {len(degs_01)}')

degs_nom = [(results[i][0], results[i][1], pvals[i], padj[i]) for i in range(len(results)) if pvals[i] < 0.05]
print(f'DEGs (nominal p<0.05): {len(degs_nom)}')

# Save DEG lists
with open('GSE301375_degs_nominal.txt', 'w') as f:
    f.write('gene\tlog2fc\tpval\tpadj\n')
    for g, lfc, p, pa in sorted(degs_nom, key=lambda x: x[2]):
        f.write(f'{g}\t{lfc:.4f}\t{p:.6e}\t{pa:.6e}\n')

with open('GSE301375_degs_fdr05.txt', 'w') as f:
    f.write('gene\tlog2fc\tpval\tpadj\n')
    for g, lfc, p, pa in sorted(degs, key=lambda x: x[2]):
        f.write(f'{g}\t{lfc:.4f}\t{p:.6e}\t{pa:.6e}\n')

# =====================================================
# LOAD BIVALENT DOMAINS
# =====================================================

bivalent = set()
with open(os.path.join('..', 'bivalent_domains_court_arnaud.txt'), 'r') as f:
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

def mouse_to_human(gene):
    return gene.upper()

# =====================================================
# ENRICHMENT ANALYSIS
# =====================================================
print()
print('=' * 60)
print('BIVALENT ENRICHMENT ANALYSIS')
print('=' * 60)

# Background
all_expressed = [results[i][0] for i in range(len(results))]
all_expressed_human = set(mouse_to_human(g) for g in all_expressed)
bg_bivalent = all_expressed_human & bivalent
bg_total = len(all_expressed_human)
bg_biv_count = len(bg_bivalent)
bg_ratio = bg_biv_count / bg_total

print(f'Background: {bg_total} expressed genes, {bg_biv_count} bivalent ({100*bg_ratio:.1f}%)')

# DEGs nominal
deg_genes_human = set(mouse_to_human(g) for g, _, _, _ in degs_nom)
deg_bivalent = deg_genes_human & bivalent
deg_biv_count = len(deg_bivalent)
deg_total = len(deg_genes_human)
deg_ratio = deg_biv_count / deg_total if deg_total > 0 else 0

print(f'DEGs (nominal p<0.05): {deg_total} genes, {deg_biv_count} bivalent ({100*deg_ratio:.1f}%)')
enrichment = deg_ratio / bg_ratio if bg_ratio > 0 else 0
print(f'Enrichment: {enrichment:.2f}x')

a = deg_biv_count
b = deg_total - deg_biv_count
c = bg_biv_count - deg_biv_count
d = bg_total - bg_biv_count - b
oddsratio, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
print(f'Fisher exact: OR={oddsratio:.2f}, p={fisher_p:.4e}')

# FDR<0.05
if degs:
    deg05_human = set(mouse_to_human(g) for g, _, _, _ in degs)
    deg05_biv = deg05_human & bivalent
    deg05_ratio = len(deg05_biv) / max(len(deg05_human), 1)
    print(f'\nDEGs (FDR<0.05): {len(deg05_human)} genes, {len(deg05_biv)} bivalent ({100*deg05_ratio:.1f}%)')
    enrichment05 = deg05_ratio / bg_ratio if bg_ratio > 0 else 0
    print(f'Enrichment: {enrichment05:.2f}x')
    a2 = len(deg05_biv)
    b2 = len(deg05_human) - a2
    c2 = bg_biv_count - a2
    d2 = bg_total - bg_biv_count - b2
    or2, fp2 = stats.fisher_exact([[a2, b2], [c2, d2]], alternative='greater')
    print(f'Fisher exact: OR={or2:.2f}, p={fp2:.4e}')

# Bivalent DEGs list
print(f'\nBivalent DEGs (nominal p<0.05, top 30 by p-value):')
biv_degs = [(g, lfc, p, pa) for g, lfc, p, pa in sorted(degs_nom, key=lambda x: x[2]) if mouse_to_human(g) in bivalent]
for g, lfc, p, pa in biv_degs[:30]:
    direction = 'UP' if lfc > 0 else 'DOWN'
    print(f'  {g:15s} log2FC={lfc:+.3f} ({direction})  p={p:.4e}  padj={pa:.4e}')
if len(biv_degs) > 30:
    print(f'  ... and {len(biv_degs)-30} more')

# HOX genes specifically
print(f'\nHOX gene DEGs:')
hox_degs = [(g, lfc, p, pa) for g, lfc, p, pa in degs_nom if g.upper().startswith('HOX')]
if hox_degs:
    for g, lfc, p, pa in sorted(hox_degs, key=lambda x: x[2]):
        biv_tag = ' [BIVALENT]' if mouse_to_human(g) in bivalent else ''
        print(f'  {g:15s} log2FC={lfc:+.3f}  p={p:.4e}  padj={pa:.4e}{biv_tag}')
else:
    print('  None at p<0.05')

# Also check: are DOWN-regulated DEGs more enriched in bivalent than UP?
degs_down_nom = [(g, lfc, p, pa) for g, lfc, p, pa in degs_nom if lfc < 0]
degs_up_nom = [(g, lfc, p, pa) for g, lfc, p, pa in degs_nom if lfc > 0]
down_human = set(mouse_to_human(g) for g, _, _, _ in degs_down_nom)
up_human = set(mouse_to_human(g) for g, _, _, _ in degs_up_nom)
down_biv = down_human & bivalent
up_biv = up_human & bivalent
print(f'\nDirectional enrichment:')
print(f'  UP DEGs: {len(up_human)}, bivalent {len(up_biv)} ({100*len(up_biv)/max(len(up_human),1):.1f}%)')
print(f'  DOWN DEGs: {len(down_human)}, bivalent {len(down_biv)} ({100*len(down_biv)/max(len(down_human),1):.1f}%)')
print(f'  Background: {100*bg_ratio:.1f}%')

# =====================================================
# ALSO: Haimbaugh 2022 DEGs
# =====================================================
print(f'\n{"=" * 60}')
print('HAIMBAUGH 2022 (Zebrafish) - BIVALENT INTERSECTION')
print(f'{"=" * 60}')

haim_file = os.path.join('..', 'geo_haimbaugh2022', 'target_gene_hits_nominal_p05.csv')
if os.path.exists(haim_file):
    haim_genes = set()
    with open(haim_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) > 1:
                gene = parts[1].strip().upper() if len(parts) > 1 else ''
                if gene:
                    haim_genes.add(gene)
    haim_biv = haim_genes & bivalent
    print(f'Haimbaugh nominal DEGs: {len(haim_genes)} unique genes')
    print(f'Bivalent: {len(haim_biv)} ({100*len(haim_biv)/max(len(haim_genes),1):.1f}%)')
    if haim_biv:
        print(f'Bivalent hits: {sorted(haim_biv)}')

# =====================================================
# GSE288358 methylation targets
# =====================================================
print(f'\n{"=" * 60}')
print('GSE288358 METHYLATION TARGETS - BIVALENT CHECK')
print(f'{"=" * 60}')

# From our analysis: genes with significant PFAS-correlated methylation
methyl_targets = ['PHF8', 'KDM2B', 'HOXA3', 'HOXA5', 'HOXA9', 'ATP7B', 'SLC31A1',
                  'LAMP1', 'DNAJB1', 'DHCR24', 'FTL', 'CTSD', 'IREB2', 'NPC1',
                  'PRNP', 'EZH2', 'SUZ12']
for g in methyl_targets:
    status = 'BIVALENT' if g in bivalent else 'not bivalent'
    print(f'  {g:15s}: {status}')

# Liu 2022 HOX hits
print(f'\nLiu 2022 cord blood hypermethylated HOX genes:')
liu_hox = ['HOXA3', 'HOXB3', 'HOXB6', 'HOXC6', 'HOXC8']
for g in liu_hox:
    status = 'BIVALENT' if g in bivalent else 'not bivalent'
    print(f'  {g:15s}: {status}')
