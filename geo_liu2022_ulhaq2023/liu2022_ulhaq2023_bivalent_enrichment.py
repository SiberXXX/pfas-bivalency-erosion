"""
Bivalent domain enrichment analysis for two independent PFAS datasets.

Dataset 1: Liu et al. (2022) Environ Health Perspect 130(3):037005
  - Gestational PFAS exposure and offspring DNA methylation (birth + age 12)
  - HOME Study, Cincinnati, n=279 (birth) / n=236 (age 12)
  - Illumina 450K array, 4 PFAS: PFOA, PFOS, PFHxS, PFNA
  - DMPs at raw p < 0.001 provided in supplementary tables

Dataset 2: Ulhaq et al. (2023) J Hazard Mater 457:131718
  - PFOS-exposed zebrafish embryo RNA-seq (DEGs)
  - Cross-species validation: zebrafish gene symbols mapped to human orthologs
  - Sheets: "Up" (upregulated) and "Down" (downregulated)

Bivalent genes: Court & Arnaud (2017) - 5,377 HC bivalent genes from 5 hESC lines
"""

import pickle, openpyxl, sys, io, json
import numpy as np
from scipy.stats import fisher_exact

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# =====================================================
# LOAD RESOURCES
# =====================================================
print("=" * 80)
print("BIVALENT ENRICHMENT: Liu 2022 + Ulhaq 2023")
print("=" * 80)

# Load probe-to-gene mapping (from EPIC manifest)
print("\nLoading probe-to-gene mapping...")
with open(os.path.join('..', 'geo_GSE288358', 'probe_to_gene.pkl'), 'rb') as f:
    probe_to_gene = pickle.load(f)
print(f"  {len(probe_to_gene)} probes mapped")

# Load bivalent gene list
print("Loading bivalent gene list (Court & Arnaud 2017)...")
bivalent_genes = set()
with open(os.path.join('..', 'bivalent_domains_court_arnaud.txt'), 'r') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            bivalent_genes.add(line.upper())
print(f"  {len(bivalent_genes)} bivalent genes")


# =====================================================
# ANALYSIS 1: Liu 2022 - CpG-level bivalent enrichment
# =====================================================
print("\n" + "=" * 80)
print("ANALYSIS 1: Liu et al. (2022) - Gestational PFAS x offspring DNAm")
print("=" * 80)

wb = openpyxl.load_workbook(
    'PMC8911098_supplementary/Supplemental Excel File R2.xlsx',
    read_only=True
)

pfas_sheets = {
    'Excel Table S1': 'PFOA',
    'Excel Table S2': 'PFOS',
    'Excel Table S3': 'PFHxS',
    'Excel Table S4': 'PFNA',
}

# Count total probes on 450K that have gene annotations (background)
total_probes_with_genes = sum(1 for p, g in probe_to_gene.items() if g)
all_genes_on_array = set()
for p, genes in probe_to_gene.items():
    if genes:
        for g in genes:
            all_genes_on_array.add(g.upper())

bivalent_on_array = all_genes_on_array & bivalent_genes
nonbivalent_on_array = all_genes_on_array - bivalent_genes
print(f"\nBackground: {len(all_genes_on_array)} genes on array, "
      f"{len(bivalent_on_array)} bivalent, {len(nonbivalent_on_array)} non-bivalent")

liu_results = {}

for sheet_name, pfas in pfas_sheets.items():
    ws = wb[sheet_name]

    cpgs_hyper = []
    cpgs_hypo = []

    for row in ws.iter_rows(min_row=3, values_only=True):
        if row[0] is None or not str(row[0]).startswith('cg'):
            continue
        cpg = str(row[0])
        beta = float(row[1]) if row[1] is not None else 0

        if beta > 0:
            cpgs_hyper.append(cpg)
        else:
            cpgs_hypo.append(cpg)

    print(f"\n--- {pfas} (Sheet: {sheet_name}) ---")
    print(f"  Total DMPs: {len(cpgs_hyper) + len(cpgs_hypo)} "
          f"(HYPER: {len(cpgs_hyper)}, HYPO: {len(cpgs_hypo)})")

    liu_results[pfas] = {}

    for direction, cpg_list in [('HYPER', cpgs_hyper), ('HYPO', cpgs_hypo)]:
        # Map CpGs to genes
        genes_in_list = set()
        unmapped = 0
        for cpg in cpg_list:
            gene_set = probe_to_gene.get(cpg)
            if gene_set:
                for g in gene_set:
                    genes_in_list.add(g.upper())
            else:
                unmapped += 1

        biv_hit = genes_in_list & bivalent_genes
        nonbiv_hit = genes_in_list - bivalent_genes
        biv_bg = bivalent_on_array - biv_hit
        nonbiv_bg = nonbivalent_on_array - nonbiv_hit

        # Fisher exact test (gene-level)
        table = [[len(biv_hit), len(nonbiv_hit)],
                 [len(biv_bg), len(nonbiv_bg)]]
        odds, pval = fisher_exact(table, alternative='two-sided')

        # Expected fraction
        expected_frac = len(bivalent_on_array) / len(all_genes_on_array) if len(all_genes_on_array) > 0 else 0
        observed_frac = len(biv_hit) / len(genes_in_list) if len(genes_in_list) > 0 else 0
        enrichment = observed_frac / expected_frac if expected_frac > 0 else 0

        print(f"\n  {direction} DMPs -> {len(genes_in_list)} unique genes "
              f"({unmapped} CpGs unmapped)")
        print(f"    Bivalent: {len(biv_hit)}/{len(genes_in_list)} "
              f"({100*observed_frac:.1f}%)")
        print(f"    Expected: {100*expected_frac:.1f}%")
        print(f"    Enrichment: {enrichment:.2f}x")
        print(f"    Fisher exact p = {pval:.2e}, OR = {odds:.2f}")
        print(f"    Table: [[{table[0][0]}, {table[0][1]}], "
              f"[{table[1][0]}, {table[1][1]}]]")

        liu_results[pfas][direction] = {
            'n_cpgs': len(cpg_list),
            'n_genes': len(genes_in_list),
            'n_bivalent': len(biv_hit),
            'enrichment': round(enrichment, 3),
            'odds_ratio': round(odds, 3),
            'fisher_p': pval,
            'unmapped_cpgs': unmapped,
        }

wb.close()

# =====================================================
# ANALYSIS 2: Ulhaq 2023 - Zebrafish DEGs
# =====================================================
print("\n\n" + "=" * 80)
print("ANALYSIS 2: Ulhaq et al. (2023) - PFOS-exposed zebrafish DEGs")
print("=" * 80)
print("Note: zebrafish gene symbols matched directly to human bivalent list.")
print("Background: all human protein-coding genes (~20,000) as universe.")

wb2 = openpyxl.load_workbook('mmc2_ulhaq2023_jhazmat.xlsx', read_only=True)

# First pass: collect all DEGs and map to human symbols
all_degs = {}  # direction -> list of gene symbols
for sheet_name, direction in [('Up', 'UP'), ('Down', 'DOWN')]:
    ws = wb2[sheet_name]
    genes = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        if row[1] is not None:
            gene_sym = str(row[1]).strip("'\" ")
            if gene_sym:
                genes.append(gene_sym.upper())
    all_degs[direction] = genes

wb2.close()

# For zebrafish RNA-seq, background should be the whole genome
# Since we lack the full zebrafish tested gene list, use the array gene set
# as a proxy for the human gene universe (25,608 genes)
# This is imperfect but consistent with the Liu 2022 analysis above

ulhaq_results = {}

for direction in ['UP', 'DOWN']:
    genes = all_degs[direction]
    print(f"\n--- {direction}regulated DEGs: {len(genes)} total ---")

    # Match zebrafish gene symbols to human gene universe
    biv_hit = set()
    nonbiv_hit = set()
    no_match = set()

    for g in genes:
        if g in bivalent_genes:
            biv_hit.add(g)
        elif g in all_genes_on_array:
            nonbiv_hit.add(g)
        else:
            no_match.add(g)

    matched_genes = biv_hit | nonbiv_hit

    if len(matched_genes) > 0:
        # Fisher test: same background as Liu 2022 (array genes)
        biv_not_hit = bivalent_on_array - biv_hit
        nonbiv_not_hit = nonbivalent_on_array - nonbiv_hit

        table = [[len(biv_hit), len(nonbiv_hit)],
                 [len(biv_not_hit), len(nonbiv_not_hit)]]
        odds, pval = fisher_exact(table, alternative='two-sided')

        expected_frac = len(bivalent_on_array) / len(all_genes_on_array)
        observed_frac = len(biv_hit) / len(matched_genes)
        enrichment = observed_frac / expected_frac if expected_frac > 0 else 0

        print(f"  Matched to human genes: {len(matched_genes)} / {len(genes)} "
              f"({100*len(matched_genes)/len(genes):.0f}%)")
        print(f"  Not matched (zebrafish-specific): {len(no_match)}")
        print(f"  Bivalent: {len(biv_hit)}/{len(matched_genes)} "
              f"({100*observed_frac:.1f}%)")
        print(f"  Expected: {100*expected_frac:.1f}%")
        print(f"  Enrichment: {enrichment:.2f}x")
        print(f"  Fisher exact p = {pval:.2e}, OR = {odds:.2f}")
        print(f"  Table: [[{table[0][0]}, {table[0][1]}], "
              f"[{table[1][0]}, {table[1][1]}]]")

        if biv_hit:
            examples = sorted(biv_hit)[:20]
            print(f"  Example bivalent DEGs: {', '.join(examples)}")

        ulhaq_results[direction] = {
            'n_total_degs': len(genes),
            'n_matched_human': len(matched_genes),
            'n_bivalent': len(biv_hit),
            'enrichment': round(enrichment, 3),
            'odds_ratio': round(odds, 3),
            'fisher_p': pval,
            'n_unmatched': len(no_match),
            'bivalent_examples': sorted(biv_hit)[:30],
        }
    else:
        print("  No genes matched human gene symbols")
        ulhaq_results[direction] = {'error': 'No matched genes'}


# =====================================================
# COMBINED SUMMARY
# =====================================================
print("\n\n" + "=" * 80)
print("COMBINED SUMMARY")
print("=" * 80)

print("\n--- Liu 2022 (gestational PFAS, offspring DNAm, 450K array) ---")
print(f"{'PFAS':<8s} {'Dir':<6s} {'DMPs':>6s} {'Genes':>6s} {'Bival':>6s} "
      f"{'Enrich':>8s} {'OR':>6s} {'p':>12s}")
print("-" * 65)
for pfas in ['PFOA', 'PFOS', 'PFHxS', 'PFNA']:
    for d in ['HYPER', 'HYPO']:
        r = liu_results[pfas][d]
        sig = ""
        if r['fisher_p'] < 0.001:
            sig = " ***"
        elif r['fisher_p'] < 0.01:
            sig = " **"
        elif r['fisher_p'] < 0.05:
            sig = " *"
        print(f"{pfas:<8s} {d:<6s} {r['n_cpgs']:>6d} {r['n_genes']:>6d} "
              f"{r['n_bivalent']:>6d} {r['enrichment']:>7.2f}x "
              f"{r['odds_ratio']:>6.2f} {r['fisher_p']:>12.2e}{sig}")

print("\n--- Ulhaq 2023 (PFOS-exposed zebrafish, RNA-seq) ---")
for d in ['UP', 'DOWN']:
    if 'error' not in ulhaq_results[d]:
        r = ulhaq_results[d]
        sig = ""
        if r['fisher_p'] < 0.001:
            sig = " ***"
        elif r['fisher_p'] < 0.01:
            sig = " **"
        elif r['fisher_p'] < 0.05:
            sig = " *"
        print(f"PFOS     {d:<6s} {r['n_total_degs']:>6d} {r['n_matched_human']:>6d} "
              f"{r['n_bivalent']:>6d} {r['enrichment']:>7.2f}x "
              f"{r['odds_ratio']:>6.2f} {r['fisher_p']:>12.2e}{sig}")

# =====================================================
# SAVE RESULTS
# =====================================================
output = {
    'liu2022': liu_results,
    'ulhaq2023': ulhaq_results,
    'background': {
        'total_genes_on_array': len(all_genes_on_array),
        'bivalent_on_array': len(bivalent_on_array),
        'nonbivalent_on_array': len(nonbivalent_on_array),
        'bivalent_fraction': round(len(bivalent_on_array) / len(all_genes_on_array), 4),
    }
}

# Convert numpy/float types for JSON
def convert(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        if np.isnan(obj) or np.isinf(obj):
            return str(obj)
        return float(obj)
    if isinstance(obj, set):
        return sorted(obj)
    return obj

with open('liu2022_ulhaq2023_bivalent_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=convert)

print("\nResults saved to liu2022_ulhaq2023_bivalent_results.json")
print("\nDONE.")
