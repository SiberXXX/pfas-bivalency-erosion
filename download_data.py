#!/usr/bin/env python3
"""
Download all public datasets required to reproduce the analyses.

Sources: NCBI GEO, CDC NHANES, RCSB PDB.
All data is publicly available. No authentication required.

Usage:
    python download_data.py          # download everything
    python download_data.py geo      # GEO datasets only
    python download_data.py nhanes   # NHANES XPT files only
    python download_data.py pdb      # PDB structures only
"""

import os
import sys
import urllib.request
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# =====================================================================
# GEO datasets
# =====================================================================

GEO_FILES = [
    # GSE254408 - Running 2024, SH-SY5Y neurons
    {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254408/suppl/GSE254408_all_samples_raw_counts.txt.gz',
        'dest': 'geo_GSE254408/GSE254408_all_samples_raw_counts.txt.gz',
    },
    # GSE288358 - Everson 2025, placental EWAS
    {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288358/suppl/GSE288358_Beta_Matrix.csv.gz',
        'dest': 'geo_GSE288358/GSE288358_Beta_Matrix.csv.gz',
    },
    {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288358/matrix/GSE288358_series_matrix.txt.gz',
        'dest': 'geo_GSE288358/GSE288358_series_matrix.txt.gz',
    },
    # GSE79329 - van den Dungen 2017, Dutch men
    {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79329/suppl/GSE79329_POP_levels.txt.gz',
        'dest': 'geo_GSE79329/GSE79329_POP_levels.txt.gz',
    },
    {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79329/matrix/GSE79329_series_matrix.txt.gz',
        'dest': 'geo_GSE79329/GSE79329_series_matrix.txt.gz',
    },
]

# =====================================================================
# NHANES XPT files (CDC)
# =====================================================================

CDC_BASE = 'https://wwwn.cdc.gov/Nchs/Nhanes'

NHANES_FILES = [
    # Cycle E (2007-2008) - iron/ferritin analysis
    ('2007-2008', 'PFC_E.XPT',    'nhanes/cycle_e_2007/PFC_E.XPT'),
    ('2007-2008', 'DEMO_E.XPT',   'nhanes/cycle_e_2007/DEMO_E.XPT'),
    ('2007-2008', 'BIOPRO_E.XPT', 'nhanes/cycle_e_2007/BIOPRO_E.XPT'),
    ('2007-2008', 'FERTIN_E.XPT', 'nhanes/cycle_e_2007/FERTIN_E.XPT'),
    ('2007-2008', 'TFR_E.XPT',    'nhanes/cycle_e_2007/TFR_E.XPT'),
    # Cycle F (2009-2010) - iron/ferritin analysis
    ('2009-2010', 'PFC_F.XPT',    'nhanes/cycle_f_2009/PFC_F.XPT'),
    ('2009-2010', 'DEMO_F.XPT',   'nhanes/cycle_f_2009/DEMO_F.XPT'),
    ('2009-2010', 'BIOPRO_F.XPT', 'nhanes/cycle_f_2009/BIOPRO_F.XPT'),
    ('2009-2010', 'FERTIN_F.XPT', 'nhanes/cycle_f_2009/FERTIN_F.XPT'),
    ('2009-2010', 'TFR_F.XPT',    'nhanes/cycle_f_2009/TFR_F.XPT'),
    # Cycle G (2011-2012) - copper analysis
    ('2011-2012', 'PFC_G.XPT',    'nhanes/PFC_G.XPT'),
    ('2011-2012', 'CUSEZN_G.XPT', 'nhanes/CUSEZN_G.XPT'),
    ('2011-2012', 'DEMO_G.XPT',   'nhanes/DEMO_G.XPT'),
    ('2011-2012', 'BIOPRO_G.XPT', 'nhanes/BIOPRO_G.XPT'),
    # Cycle H (2013-2014) - copper analysis
    ('2013-2014', 'PFAS_H.XPT',   'nhanes/cycle_h_2013/PFAS_H.XPT'),
    ('2013-2014', 'CUSEZN_H.XPT', 'nhanes/cycle_h_2013/CUSEZN_H.XPT'),
    ('2013-2014', 'DEMO_H.XPT',   'nhanes/cycle_h_2013/DEMO_H.XPT'),
    ('2013-2014', 'BIOPRO_H.XPT', 'nhanes/cycle_h_2013/BIOPRO_H.XPT'),
    # Cycle I (2015-2016) - copper analysis
    ('2015-2016', 'PFAS_I.XPT',   'nhanes/cycle_i_2015/PFAS_I.XPT'),
    ('2015-2016', 'CUSEZN_I.XPT', 'nhanes/cycle_i_2015/CUSEZN_I.XPT'),
    ('2015-2016', 'DEMO_I.XPT',   'nhanes/cycle_i_2015/DEMO_I.XPT'),
    ('2015-2016', 'BIOPRO_I.XPT', 'nhanes/cycle_i_2015/BIOPRO_I.XPT'),
]

# =====================================================================
# PDB structures (RCSB)
# =====================================================================

PDB_IDS = ['2OX0', '5FUP', '3N9N', '2XUE']


def download_file(url, dest_path):
    """Download a file with progress indication. Skip if already exists."""
    full_path = os.path.join(SCRIPT_DIR, dest_path)

    if os.path.exists(full_path):
        size_mb = os.path.getsize(full_path) / (1024 * 1024)
        print(f"  SKIP  {dest_path} ({size_mb:.1f} MB, already exists)")
        return True

    os.makedirs(os.path.dirname(full_path), exist_ok=True)

    print(f"  GET   {dest_path} ... ", end='', flush=True)
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Python/analysis-scripts (academic research)'
        })
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
            with open(full_path, 'wb') as f:
                f.write(data)
        size_mb = len(data) / (1024 * 1024)
        print(f"OK ({size_mb:.1f} MB)")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False


def download_geo():
    print("\n" + "=" * 60)
    print("GEO DATASETS (NCBI)")
    print("=" * 60)
    ok = 0
    for item in GEO_FILES:
        if download_file(item['url'], item['dest']):
            ok += 1
        time.sleep(0.5)
    print(f"\nGEO: {ok}/{len(GEO_FILES)} files downloaded")

    # GSE301375 needs special handling (multiple count files)
    dest_dir = os.path.join(SCRIPT_DIR, 'geo_GSE301375', 'GSE301375_counts')
    if not os.path.exists(dest_dir) or not os.listdir(dest_dir):
        print(f"\n  NOTE: GSE301375 count files must be downloaded manually.")
        print(f"  Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301375")
        print(f"  Place count files in: geo_GSE301375/GSE301375_counts/")
    return ok


def download_nhanes():
    print("\n" + "=" * 60)
    print("NHANES XPT FILES (CDC)")
    print("=" * 60)
    ok = 0
    for cycle, fname, dest in NHANES_FILES:
        url = f"{CDC_BASE}/{cycle}/{fname}"
        if download_file(url, dest):
            ok += 1
        time.sleep(0.3)
    print(f"\nNHANES: {ok}/{len(NHANES_FILES)} files downloaded")
    return ok


def download_pdb():
    print("\n" + "=" * 60)
    print("PDB STRUCTURES (RCSB)")
    print("=" * 60)
    ok = 0
    for pdb_id in PDB_IDS:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        dest = f"docking/{pdb_id}.pdb"
        if download_file(url, dest):
            ok += 1
        time.sleep(0.3)
    print(f"\nPDB: {ok}/{len(PDB_IDS)} files downloaded")
    return ok


def main():
    print("Data downloader for PFAS Bivalency Erosion analysis")
    print("All sources are public. No authentication required.")

    targets = sys.argv[1:] if len(sys.argv) > 1 else ['geo', 'nhanes', 'pdb']

    total_ok = 0
    total_files = 0

    if 'geo' in targets:
        total_ok += download_geo()
        total_files += len(GEO_FILES)

    if 'nhanes' in targets:
        total_ok += download_nhanes()
        total_files += len(NHANES_FILES)

    if 'pdb' in targets:
        total_ok += download_pdb()
        total_files += len(PDB_IDS)

    print("\n" + "=" * 60)
    print(f"TOTAL: {total_ok}/{total_files} files downloaded successfully")

    # Remind about manual steps
    print("\nManual steps still required:")
    print("  1. GSE301375 count files (see note above)")
    print("  2. Haimbaugh 2022 supplementary Excel tables (Table S3-S5)")
    print("     -> place in geo_haimbaugh2022/")
    print("  3. Liu 2022 supplementary tables")
    print("     -> place in geo_liu2022_ulhaq2023/")
    print("  4. Ulhaq 2023 supplementary Excel")
    print("     -> place in geo_liu2022_ulhaq2023/")
    print("\nThen run:  python qc/download_ensembl_mapping.py")
    print("           python geo_GSE288358/run_analysis_v2.py")
    print("=" * 60)


if __name__ == '__main__':
    main()
