"""
GSE79329: Genome-wide bivalent enrichment analysis
Dutch men (n=34), leukocytes, 450K array, 5 PFAS
Spearman correlation for each probe, then Fisher test for bivalent enrichment
"""
import gzip, pickle, re, sys, io, json
import numpy as np
from scipy import stats
from scipy.stats import fisher_exact

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 80)
print("GSE79329: GENOME-WIDE BIVALENT ENRICHMENT")
print("Dutch men consuming eel, leukocytes, 450K, n=34")
print("=" * 80)

# Load probe-to-gene mapping
print("\nLoading probe-to-gene mapping...")
with open(os.path.join('..', 'geo_GSE288358', 'probe_to_gene.pkl'), 'rb') as f:
    probe_to_gene = pickle.load(f)
print(f"  {len(probe_to_gene)} EPIC probes (450K is subset)")

# Load bivalent gene list
print("Loading bivalent genes...")
bivalent_genes = set()
with open(os.path.join('..', 'bivalent_domains_court_arnaud.txt'), 'r') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            bivalent_genes.add(line.upper())
print(f"  {len(bivalent_genes)} bivalent genes")

# Parse PFAS concentrations
print("\nParsing PFAS concentrations...")
with gzip.open('GSE79329_POP_levels.txt.gz', 'rt') as f:
    content = f.read()

parts_raw = content.split('Participant_')
header_text = parts_raw[0]
col_names_raw = re.findall(r'"([^"]+)"', header_text)
col_names = [c.split('\n')[0].strip() for c in col_names_raw]

pfas_cols = {}
for i, c in enumerate(col_names):
    if c in ['PFOA', 'PFNA', 'PFDA', 'PFHxS', 'PFOS']:
        pfas_cols[c] = i

pfas_data = {k: [] for k in pfas_cols}
for part in parts_raw[1:]:
    fields = part.strip().split('\n')[0].split('\t')
    for pfas_name, col_idx in pfas_cols.items():
        try:
            val = fields[col_idx + 1].strip()
            if val in ['n.d.', 'ND', '', 'NA', '<LOD']:
                pfas_data[pfas_name].append(np.nan)
            else:
                pfas_data[pfas_name].append(float(val))
        except (IndexError, ValueError):
            pfas_data[pfas_name].append(np.nan)

for k in pfas_data:
    pfas_data[k] = np.array(pfas_data[k])

n_participants = len(pfas_data['PFOS'])
print(f"  N = {n_participants}")
for pf in ['PFOS', 'PFOA', 'PFNA', 'PFDA', 'PFHxS']:
    v = pfas_data[pf][~np.isnan(pfas_data[pf])]
    print(f"  {pf}: median={np.median(v):.1f}, range={v.min():.1f}-{v.max():.1f}")

# Read ALL probes from series matrix (not just focus genes)
print("\nReading full 450K beta matrix...")
beta_matrix = {}
sample_order = []

with gzip.open('GSE79329_series_matrix.txt.gz', 'rt') as f:
    in_data = False
    for line in f:
        if line.startswith('"ID_REF"'):
            fields = line.strip().split('\t')
            sample_order = [s.strip('"') for s in fields[1:]]
            in_data = True
            continue
        if not in_data or line.startswith('!') or line.strip() == '':
            if in_data and (line.strip() == '' or line.startswith('!')):
                break
            continue
        fields = line.strip().split('\t')
        probe_id = fields[0].strip('"')
        if not probe_id.startswith('cg'):
            continue
        vals = []
        for v in fields[1:]:
            try:
                vals.append(float(v))
            except:
                vals.append(np.nan)
        beta_matrix[probe_id] = np.array(vals)

print(f"  {len(beta_matrix)} CpG probes loaded, {len(sample_order)} samples")

# Genome-wide Spearman correlation for each probe vs each PFAS
print("\nRunning genome-wide Spearman correlations...")

results_all = {}

for pfas_name in ['PFOS', 'PFOA', 'PFNA', 'PFDA', 'PFHxS']:
    pf_vals = pfas_data[pfas_name]
    pf_valid = ~np.isnan(pf_vals)

    probes_hyper = []
    probes_hypo = []
    probes_total = 0

    for probe_id, betas in beta_matrix.items():
        valid = pf_valid & ~np.isnan(betas)
        if valid.sum() < 20:
            continue
        probes_total += 1
        rho, p = stats.spearmanr(pf_vals[valid], betas[valid])
        if p < 0.05:
            if rho > 0:
                probes_hyper.append(probe_id)
            else:
                probes_hypo.append(probe_id)

    # Map probes to genes
    def probes_to_genes(probe_list):
        genes = set()
        for p in probe_list:
            g = probe_to_gene.get(p)
            if g:
                for gene in g:
                    genes.add(gene.upper())
        return genes

    genes_hyper = probes_to_genes(probes_hyper)
    genes_hypo = probes_to_genes(probes_hypo)

    # Background: all tested probes mapped to genes
    all_tested_probes = [p for p in beta_matrix if p in probe_to_gene]
    all_genes = set()
    for p in all_tested_probes:
        g = probe_to_gene.get(p)
        if g:
            for gene in g:
                all_genes.add(gene.upper())

    bivalent_on_array = all_genes & bivalent_genes
    nonbivalent_on_array = all_genes - bivalent_genes
    expected_frac = len(bivalent_on_array) / len(all_genes) if all_genes else 0

    print(f"\n--- {pfas_name} ---")
    print(f"  Tested: {probes_total} probes, {len(all_genes)} genes")
    print(f"  Background: {len(bivalent_on_array)} bivalent / {len(all_genes)} total = {100*expected_frac:.1f}%")
    print(f"  Sig p<0.05: HYPER={len(probes_hyper)} probes ({len(genes_hyper)} genes), "
          f"HYPO={len(probes_hypo)} probes ({len(genes_hypo)} genes)")

    pfas_results = {}

    for direction, gene_set in [('HYPER', genes_hyper), ('HYPO', genes_hypo)]:
        biv_hit = gene_set & bivalent_genes
        nonbiv_hit = gene_set - bivalent_genes
        biv_bg = bivalent_on_array - biv_hit
        nonbiv_bg = nonbivalent_on_array - nonbiv_hit

        table = [[len(biv_hit), len(nonbiv_hit)],
                 [len(biv_bg), len(nonbiv_bg)]]
        odds, pval = fisher_exact(table, alternative='two-sided')

        obs_frac = len(biv_hit) / len(gene_set) if gene_set else 0
        enrichment = obs_frac / expected_frac if expected_frac > 0 else 0

        print(f"  {direction}: bivalent {len(biv_hit)}/{len(gene_set)} ({100*obs_frac:.1f}%), "
              f"enrichment={enrichment:.2f}x, OR={odds:.2f}, Fisher p={pval:.2e}")

        pfas_results[direction] = {
            'n_probes': len(probes_hyper if direction == 'HYPER' else probes_hypo),
            'n_genes': len(gene_set),
            'n_bivalent': len(biv_hit),
            'enrichment': round(enrichment, 3),
            'odds_ratio': round(odds, 3),
            'fisher_p': pval,
        }

    results_all[pfas_name] = pfas_results

# Summary table
print("\n\n" + "=" * 80)
print("SUMMARY: GSE79329 Genome-wide Bivalent Enrichment")
print("=" * 80)
print(f"{'PFAS':<8s} {'Dir':<6s} {'Probes':>7s} {'Genes':>6s} {'Bival':>6s} "
      f"{'Enrich':>8s} {'OR':>6s} {'p':>12s}")
print("-" * 62)
for pfas in ['PFOS', 'PFOA', 'PFNA', 'PFDA', 'PFHxS']:
    for d in ['HYPER', 'HYPO']:
        r = results_all[pfas][d]
        sig = ""
        if r['fisher_p'] < 0.001: sig = " ***"
        elif r['fisher_p'] < 0.01: sig = " **"
        elif r['fisher_p'] < 0.05: sig = " *"
        print(f"{pfas:<8s} {d:<6s} {r['n_probes']:>7d} {r['n_genes']:>6d} "
              f"{r['n_bivalent']:>6d} {r['enrichment']:>7.2f}x "
              f"{r['odds_ratio']:>6.2f} {r['fisher_p']:>12.2e}{sig}")

# Save
output = {
    'dataset': 'GSE79329',
    'description': 'Dutch men consuming eel, leukocytes, 450K',
    'n': n_participants,
    'results': results_all,
}

def convert(obj):
    if isinstance(obj, (np.integer,)): return int(obj)
    if isinstance(obj, (np.floating, float)):
        if np.isnan(obj) or np.isinf(obj): return str(obj)
        return float(obj)
    if isinstance(obj, set): return sorted(obj)
    return obj

with open('bivalent_enrichment_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=convert)

print("\nResults saved to bivalent_enrichment_results.json")
print("DONE.")
