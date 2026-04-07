# Analysis Scripts for "Bivalency Erosion as a Unified Mechanism of PFAS Transgenerational Neurotoxicity"

Andrey F. (2026). DOI: [10.5281/zenodo.19436133](https://doi.org/10.5281/zenodo.19436133)

Every number in the preprint can be reproduced by running the scripts below. No manual steps, no hidden parameters. If a script produces a different number, the preprint is wrong.

## Quick start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Download public datasets (~2 GB, mainly GEO beta matrices)
python download_data.py

# 3. Download Ensembl-to-symbol mapping (queries mygene.info API, ~2 min)
python qc/download_ensembl_mapping.py

# 4. Build probe-to-gene index from EPIC manifest (~1 min)
python geo_GSE288358/build_probe_index.py

# 5. Run any analysis script, e.g.:
python geo_GSE254408/analyze_GSE254408_bivalent_v3.py
```

Steps 2-4 are prerequisites. After that, all analysis scripts run independently.

## What each script does

### Molecular docking (Table 1)

| Script | What it computes |
|--------|-----------------|
| `docking/analyze_kdm6b_docking.py` | Binding affinities of PFOA/PFOS/PFHxA to KDM6B active site, Fe2+ distances from docked poses. Reads `.pdbqt` output files from AutoDock Vina runs. |

**Input:** `.pdbqt` docking output files (provided) + `5FJK.pdb` (downloaded by `download_data.py`)
**Key output:** PFOS -8.3 kcal/mol, PFOA -7.9 kcal/mol

### NHANES epidemiology (Section 3.1)

| Script | What it computes |
|--------|-----------------|
| `nhanes/analysis_cycles_EF.py` | PFAS vs serum iron and ferritin across NHANES 2003-2004 and 2005-2006. OLS regression with age, sex, ethnicity covariates. BH-FDR correction. |
| `nhanes/analysis_mega_copper.py` | PFAS vs serum copper. Only NHANES 2011-2012 (cycle G) has overlapping PFAS and copper subsamples. Cycles H and I are loaded but yield zero linkable participants (different NHANES subsamples). |

**Input:** `.XPT` files (downloaded by `download_data.py`)
**Key output:** All five PFAS negatively associated with ferritin in males (all q < 0.001)

### RNA-seq bivalent enrichment (Section 3.4)

| Script | What it computes |
|--------|-----------------|
| `geo_GSE254408/analyze_GSE254408_bivalent_v3.py` | Bivalent enrichment among DEGs in PFAS-exposed SH-SY5Y neurons (Running 2024). Tests all 6 PFAS compounds. Length-corrected by stratified permutation (50,000 iterations). Bonferroni correction across 6 compounds. |
| `geo_GSE301375/analyze_bivalent_enrichment.py` | Bivalent enrichment among DEGs in PFHxA-exposed mouse cerebellum (Plunk/Bader 2025). DESeq2-style analysis from raw counts. |

**Input:** GEO series matrix (GSE254408) or raw count files (GSE301375)
**Key output:**
- Running/PFOA: DOWN DEGs 1.46-fold enriched (Bonferroni p = 0.007)
- Plunk/PFHxA: DOWN DEGs 1.43-fold enriched (Fisher p = 0.002)
- UP DEGs show no enrichment in either dataset (internal negative control)

### Placental EWAS (Table 2)

| Script | What it computes |
|--------|-----------------|
| `geo_GSE288358/build_probe_index.py` | Builds probe-to-gene mapping from Illumina EPIC manifest. Generates `probe_coords.pkl` and `probe_to_gene.pkl`. **Must run before EWAS scripts.** |
| `geo_GSE288358/run_analysis_v2.py` | Full EWAS: regresses 365K CpG betas on PFAS concentration with 10 covariates (maternal age, gestational age, sex, BMI, batch, 5 cell types). Produces probe-level and gene-level results. |
| `geo_GSE288358/ewas_bivalent_enrichment.py` | Tests whether PFAS-associated DMPs are enriched at bivalent genes. Separate Fisher tests for hyper- and hypomethylated probes. |

**Input:** `GSE288358_Beta_Matrix.csv.gz` + series matrix (downloaded by `download_data.py`)
**Key output:** Hypermethylation enriched at bivalent genes (PFOS 1.17x probe-level, p = 1.8e-54). Hypomethylation depleted (0.87-0.90x). Directional asymmetry consistent with bivalency erosion model.

### Replication cohorts (Table 2)

| Script | What it computes |
|--------|-----------------|
| `geo_GSE79329/bivalent_enrichment_genomewide.py` | Genome-wide Spearman correlation of 464K CpG probes with 5 PFAS in Dutch eel consumers (van den Dungen 2017, n=34, 450K array). |
| `geo_liu2022_ulhaq2023/liu2022_ulhaq2023_bivalent_enrichment.py` | Bivalent enrichment in published DMR gene lists from Liu 2022 (cord blood) and Ulhaq 2023 (placenta). |

**Input:** GEO beta matrix (GSE79329) or supplementary Excel tables (Liu, Ulhaq)
**Key output:** PFHxS HYPER 2.07x (p = 5.2e-261), PFOA HYPER 2.00x (p < 1e-300)

### Zebrafish transgenerational (Section 4)

| Script | What it computes |
|--------|-----------------|
| `geo_haimbaugh2022/search_degs_v2.py` | Searches Haimbaugh 2022 zebrafish DEG tables for HOX genes and JmjC demethylases (phf8, kdm6b, etc.) across F0, F1, F2 generations. |

**Input:** Excel supplementary tables (Table S3-S5) from Haimbaugh 2022
**Key output:** phf8 log2FC = -0.42 (F0), -0.79 (F1) showing amplification across generations

### Quality control

| Script | What it computes |
|--------|-----------------|
| `qc/check_length_bias.py` | Verifies gene length bias: bivalent genes are 1.4x longer than non-bivalent. Shows that length correction is necessary and sufficient. |
| `qc/check_reliability.py` | Power analysis for n=3 RNA-seq: estimates pi0 (proportion of true nulls), permutation-based FDR, effect size distributions. |
| `qc/verify_enrichment_thorough.py` | Step-by-step verification of Fisher enrichment calculation with length-stratified permutation (50,000 iterations, 20 bins). Independent check of main result. |
| `qc/download_ensembl_mapping.py` | Downloads Ensembl ID to gene symbol mapping via mygene.info API. Generates `ensembl_to_symbol.tsv`. **Must run before GSE301375 analysis.** |

## Data files

### Downloaded automatically (`python download_data.py`)

- GEO series matrices and beta matrices (GSE254408, GSE288358, GSE79329)
- NHANES `.XPT` files (PFAS, iron, ferritin, copper, demographics)
- PDB structure `5FJK.pdb` (KDM6B crystal structure)
- Illumina EPIC manifest (for probe-to-gene mapping)

### Manual downloads (journal supplementary files)

These files come from journal supplementary materials and cannot be auto-downloaded:

| File | Where to get it | Put it in |
|------|----------------|-----------|
| Haimbaugh 2022 Table S3 (`.xlsx`) | [DOI: 10.1016/j.envpol.2022.119922](https://doi.org/10.1016/j.envpol.2022.119922) supplementary | `geo_haimbaugh2022/` |
| Haimbaugh 2022 Table S4 (`.xlsx`) | Same paper | `geo_haimbaugh2022/` |
| Haimbaugh 2022 Table S5 (`.xlsx`) | Same paper | `geo_haimbaugh2022/` |
| `ensembl_to_gene.tsv` | Haimbaugh supplementary or BioMart | `geo_haimbaugh2022/` |
| Liu 2022 supplementary tables | [DOI: 10.1016/j.envint.2022.107239](https://doi.org/10.1016/j.envint.2022.107239) | `geo_liu2022_ulhaq2023/` |
| Ulhaq 2023 supplementary Excel | [DOI: 10.1016/j.envres.2023.115733](https://doi.org/10.1016/j.envres.2023.115733) | `geo_liu2022_ulhaq2023/` |
| GSE301375 count files | [GEO: GSE301375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301375) supplementary | `geo_GSE301375/GSE301375_counts/` |

### Reference gene list

`bivalent_domains_court_arnaud.txt` contains 5,377 unique gene symbols from Court & Arnaud (2017) Oncotarget 8(3):4110-4124. These are high-confidence bivalent domain genes present in all 5 hESC lines analyzed. The file also includes a full region table (5,766 rows) for reference, but scripts use only the unique gene list section.

## Structure

```
pci_scripts/
  README.md
  requirements.txt
  download_data.py                         # Auto-downloader for public datasets
  bivalent_domains_court_arnaud.txt        # 5,377 bivalent genes (Court & Arnaud 2017)

  docking/
    analyze_kdm6b_docking.py               # Table 1: docking affinities

  nhanes/
    analysis_cycles_EF.py                  # Section 3.1: PFAS vs iron/ferritin
    analysis_mega_copper.py                # Section 3.1: PFAS vs copper (cycle G only)

  geo_GSE254408/
    analyze_GSE254408_bivalent_v3.py       # Section 3.4: Running 2024, human neurons

  geo_GSE301375/
    analyze_bivalent_enrichment.py         # Section 3.4: Plunk 2025, mouse cerebellum

  geo_GSE288358/
    build_probe_index.py                   # Utility: EPIC manifest -> probe-to-gene
    run_analysis_v2.py                     # Table 2: full EWAS with covariates
    ewas_bivalent_enrichment.py            # Table 2: bivalent enrichment in EWAS

  geo_GSE79329/
    bivalent_enrichment_genomewide.py      # Table 2: Dutch cohort replication

  geo_liu2022_ulhaq2023/
    liu2022_ulhaq2023_bivalent_enrichment.py  # Table 2: Liu + Ulhaq enrichment

  geo_haimbaugh2022/
    search_degs_v2.py                      # Section 4: zebrafish HOX/JmjC genes

  qc/
    check_length_bias.py                   # Gene length bias verification
    check_reliability.py                   # RNA-seq power analysis (n=3)
    verify_enrichment_thorough.py          # Independent enrichment verification
    download_ensembl_mapping.py            # Ensembl ID to symbol mapping
```

## Troubleshooting

**"FileNotFoundError: probe_to_gene.pkl"** — Run `python geo_GSE288358/build_probe_index.py` first. This generates `.pkl` files from the EPIC manifest.

**"FileNotFoundError: ensembl_to_symbol.tsv"** — Run `python qc/download_ensembl_mapping.py` first. This queries mygene.info to map Ensembl IDs to gene symbols.

**"FileNotFoundError: *.XPT"** — Run `python download_data.py` to download NHANES/GEO/PDB data files.

**NHANES copper script shows n=0 for cycles H and I** — This is expected. NHANES measured PFAS and copper (CUSEZN) on non-overlapping subsamples in 2013-2016. Only cycle G (2011-2012) has both analytes in the same participants.

**EWAS takes a long time** — The beta matrix is 151 samples x 365K probes. Reading and computing OLS for all probes takes several minutes. This is normal.

## License

CC BY 4.0
