"""
Build probe index files from EPIC gencode annotation (Zhou lab InfiniumAnnotationV1).

Generates two pickle files used by downstream analysis scripts:
  - probe_coords.pkl: {probe_id: (chr, position)} for coordinate-based mapping
  - probe_to_gene.pkl: {probe_id: {gene_symbols}} for gene-level enrichment

Source: EPIC.hg38.manifest.gencode.v36.tsv.gz
  Zhou W, Laird PW, Shen H. Nucleic Acids Research 2017.
  https://github.com/zhou-lab/InfiniumAnnotationV1

Columns used: probeID (4), CpG_chrm (0), CpG_beg (1), genesUniq (5)
"""

import gzip, pickle, os, sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))

MANIFEST = 'EPIC.hg38.manifest.gencode.v36.tsv.gz'

if not os.path.exists(MANIFEST):
    print(f"ERROR: {MANIFEST} not found.")
    print("Run: python download_data.py geo")
    sys.exit(1)

print("Building probe indices from EPIC gencode annotation...")
print(f"  Source: {MANIFEST}")

probe_coords = {}   # {probe_id: (chr_str, position_int)}
probe_to_gene = {}  # {probe_id: set(gene_symbols)}

with gzip.open(MANIFEST, 'rt') as f:
    header = f.readline().strip().split('\t')

    # Verify expected columns
    assert header[0] == 'CpG_chrm', f"Unexpected column 0: {header[0]}"
    assert header[4] == 'probeID', f"Unexpected column 4: {header[4]}"
    assert header[5] == 'genesUniq', f"Unexpected column 5: {header[5]}"

    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue

        probe_id = parts[4]
        if not probe_id.startswith('cg'):
            continue

        chrom = parts[0]  # e.g. "chr1"
        # Strip "chr" prefix for consistency with existing code
        chrom_short = chrom.replace('chr', '')

        try:
            pos = int(parts[1])
        except ValueError:
            continue

        probe_coords[probe_id] = (chrom_short, pos)

        gene = parts[5].strip()
        if gene and gene != 'NA' and gene != '.':
            # genesUniq may contain semicolons for multiple genes
            gene_set = set()
            for g in gene.split(';'):
                g = g.strip()
                if g and g != 'NA':
                    gene_set.add(g)
            probe_to_gene[probe_id] = gene_set
        else:
            probe_to_gene[probe_id] = set()

print(f"  {len(probe_coords)} probes with coordinates")
print(f"  {sum(1 for v in probe_to_gene.values() if v)} probes with gene annotations")

with open('probe_coords.pkl', 'wb') as f:
    pickle.dump(probe_coords, f)
print(f"  Saved probe_coords.pkl")

with open('probe_to_gene.pkl', 'wb') as f:
    pickle.dump(probe_to_gene, f)
print(f"  Saved probe_to_gene.pkl")

print("DONE.")
