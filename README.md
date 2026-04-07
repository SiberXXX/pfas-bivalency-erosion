# Analysis Scripts for "Bivalency Erosion as a Unified Mechanism of PFAS Transgenerational Neurotoxicity"

Andrey F. (2026). DOI: 10.5281/zenodo.19436133

## Structure

```
pci_scripts/
  bivalent_domains_court_arnaud.txt   # Reference gene list (Court & Arnaud 2017)
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

### GEO datasets (download before running)
- GSE254408: place `GSE254408_all_samples_raw_counts.txt.gz` in `geo_GSE254408/`
- GSE301375: place count files in `geo_GSE301375/GSE301375_counts/`
- GSE288358: place `GSE288358_Beta_Matrix.csv.gz` and `GSE288358_series_matrix.txt` in `geo_GSE288358/`
- GSE79329: place `GSE79329_POP_levels.txt.gz` and `GSE79329_series_matrix.txt.gz` in `geo_GSE79329/`

### NHANES XPT files (download from CDC)
Place XPT files in `nhanes/` with subdirectories per cycle:
- `nhanes/PFC_G.XPT`, `CUSEZN_G.XPT`, `DEMO_G.XPT`, `BIOPRO_G.XPT` (2011-2012)
- `nhanes/cycle_h_2013/PFAS_H.XPT`, etc. (2013-2014)
- `nhanes/cycle_i_2015/PFAS_I.XPT`, etc. (2015-2016)
- `nhanes/cycle_e_2007/PFC_E.XPT`, etc. (2007-2008)
- `nhanes/cycle_f_2009/PFC_F.XPT`, etc. (2009-2010)

### Docking output
AutoDock Vina PDBQT output files should be in `docking/`. PDB structures from RCSB: 2OX0, 5FUP, 3N9N, 2XUE, 1QLX, 1F41.

### Supplementary data
- Haimbaugh 2022: place Excel tables (Table S3-S5) and `ensembl_to_gene.tsv` in `geo_haimbaugh2022/`
- Liu 2022: supplementary tables referenced by `geo_liu2022_ulhaq2023/` script
- Ulhaq 2023: supplementary Excel in `geo_liu2022_ulhaq2023/`

## Running order

1. `qc/download_ensembl_mapping.py` (generates Ensembl-to-symbol mapping)
2. `geo_GSE288358/run_analysis_v2.py` (generates `probe_to_gene.pkl` used by other scripts)
3. All other scripts can run independently after steps 1-2

## Script descriptions

| Script | Produces | Preprint reference |
|--------|----------|--------------------|
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
