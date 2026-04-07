# Analysis Scripts for "Bivalency Erosion as a Unified Mechanism of PFAS Transgenerational Neurotoxicity"

Andrey F. (2026). DOI: 10.5281/zenodo.19436133

## Structure

```
pci_scripts/
  bivalent_domains_court_arnaud.txt   # Reference gene list (Court & Arnaud 2017)
  download_data.py                    # Auto-downloader for public datasets
  requirements.txt

  docking/                # Molecular docking (Table 1)
  nhanes/                 # NHANES epidemiological analysis (Section 3.1)
  geo_GSE254408/          # Running 2024, SH-SY5Y neurons (Section 3.4)
  geo_GSE301375/          # Plunk 2025, mouse cerebellum (Section 3.4)
  geo_GSE288358/          # Everson 2025, placental EWAS (Table 2)
  geo_GSE79329/           # van den Dungen 2017, leukocytes (Table 2)
  geo_liu2022_ulhaq2023/  # Liu 2022 + Ulhaq 2023 enrichment (Table 2)
  geo_haimbaugh2022/      # Haimbaugh 2022, zebrafish DEGs (Section 4)
  qc/                     # Quality control and verification scripts
```

## Dependencies

```
pip install -r requirements.txt
```

## Data acquisition

Most public data (GEO, NHANES, PDB, EPIC manifest) can be downloaded automatically:

```
python download_data.py
```

This fetches 33 files from NCBI GEO, CDC NHANES, RCSB PDB, and GitHub (Zhou lab EPIC annotation). Run with `geo`, `nhanes`, or `pdb` argument to download selectively.

### Manual downloads (journal supplementary, not automatable)
- Haimbaugh 2022: Excel tables (Table S3-S5) and `ensembl_to_gene.tsv` -> `geo_haimbaugh2022/`
- Liu 2022: supplementary tables -> `geo_liu2022_ulhaq2023/`
- Ulhaq 2023: supplementary Excel -> `geo_liu2022_ulhaq2023/`
- GSE301375: count files -> `geo_GSE301375/GSE301375_counts/`

## Running order

1. `python download_data.py` (downloads public datasets)
2. `python qc/download_ensembl_mapping.py` (generates `ensembl_to_symbol.tsv` via mygene.info API)
3. `python geo_GSE288358/build_probe_index.py` (generates `probe_coords.pkl` and `probe_to_gene.pkl` from EPIC manifest)
4. `python geo_GSE288358/run_analysis_v2.py` (full EWAS with covariates)
5. All other scripts can run independently after steps 1-3

Note: `.pkl` files are generated artifacts and are not tracked in git. Run step 3 to produce them.

## Script descriptions

| Script | Produces | Preprint reference |
|--------|----------|-------------------|
| `geo_GSE288358/build_probe_index.py` | Probe-to-gene mapping from EPIC manifest | Utility |
| `docking/analyze_kdm6b_docking.py` | Docking affinities, Fe2+ distances | Table 1 |
| `nhanes/analysis_mega_copper.py` | PFAS vs serum copper (3 cycles) | Section 3.1 |
| `nhanes/analysis_cycles_EF.py` | PFAS vs iron/ferritin (2 cycles) | Section 3.1 |
| `geo_GSE254408/analyze_GSE254408_bivalent_v3.py` | Bivalent enrichment in Running DEGs | Section 3.4 |
| `geo_GSE301375/analyze_bivalent_enrichment.py` | Bivalent enrichment in Plunk DEGs | Section 3.4 |
| `geo_GSE288358/ewas_bivalent_enrichment.py` | Bivalent enrichment in placental EWAS | Table 2 |
| `geo_GSE288358/run_analysis_v2.py` | Full EWAS with covariates | Table 2 |
| `geo_GSE79329/bivalent_enrichment_genomewide.py` | Replication in Dutch cohort | Table 2 |
| `geo_liu2022_ulhaq2023/liu2022_ulhaq2023_bivalent_enrichment.py` | Liu + Ulhaq enrichment | Table 2 |
| `geo_haimbaugh2022/search_degs_v2.py` | HOX/JmjC genes across generations | Section 4 |
| `qc/check_length_bias.py` | Gene length bias verification | Section 3.4 |
| `qc/check_reliability.py` | RNA-seq power analysis (n=3) | Section 3.2 |
| `qc/verify_enrichment_thorough.py` | Step-by-step Fisher test verification | Section 3.4 |
| `qc/download_ensembl_mapping.py` | Ensembl ID to gene symbol mapping | Utility |

## License

CC BY 4.0
