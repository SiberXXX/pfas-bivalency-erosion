"""
NHANES 2007-2008 (E) + 2009-2010 (F): PFAS vs albumin, iron, ferritin
EXPLORATORY / HYPOTHESIS-GENERATING — no serum copper in these cycles

Replicates the analytical approach from cycles G/H/I (2011-2016)
but limited to variables available: albumin (LBXSAL), iron (LBXSIR),
ferritin (LBXFER), transferrin receptor (LBXTFR).
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import warnings
warnings.filterwarnings('ignore')
import os
BASE = os.path.dirname(os.path.abspath(__file__))

# ── Load and merge ──────────────────────────────────────────────────

def load_cycle(suffix, path):
    pfc = pd.read_sas(f'{BASE}/{path}/PFC_{suffix}.XPT', format='xport')
    demo = pd.read_sas(f'{BASE}/{path}/DEMO_{suffix}.XPT', format='xport')
    bio = pd.read_sas(f'{BASE}/{path}/BIOPRO_{suffix}.XPT', format='xport')
    fer = pd.read_sas(f'{BASE}/{path}/FERTIN_{suffix}.XPT', format='xport')
    tfr_ = pd.read_sas(f'{BASE}/{path}/TFR_{suffix}.XPT', format='xport')

    # PFAS variables
    pfas_vars = ['SEQN', 'LBXPFOA', 'LBXPFOS', 'LBXPFHS', 'LBXPFDE', 'LBXPFNA']
    pfas_vars = [v for v in pfas_vars if v in pfc.columns]

    # Demographics
    demo_vars = ['SEQN', 'RIAGENDR', 'RIDAGEYR', 'RIDRETH1', 'INDFMPIR']
    demo_vars = [v for v in demo_vars if v in demo.columns]

    # Merge
    df = pfc[pfas_vars].merge(demo[demo_vars], on='SEQN')
    df = df.merge(bio[['SEQN', 'LBXSAL', 'LBXSIR']], on='SEQN')
    df = df.merge(fer[['SEQN', 'LBXFER']], on='SEQN', how='left')
    df = df.merge(tfr_[['SEQN', 'LBXTFR']], on='SEQN', how='left')
    df['CYCLE'] = suffix
    return df

df_e = load_cycle('E', 'cycle_e_2007')
df_f = load_cycle('F', 'cycle_f_2009')
df = pd.concat([df_e, df_f], ignore_index=True)

print(f"Pooled cycles E+F: n={len(df)}")
print(f"  With PFOA: {df['LBXPFOA'].notna().sum()}")
print(f"  With albumin: {df['LBXSAL'].notna().sum()}")
print(f"  With iron: {df['LBXSIR'].notna().sum()}")
print(f"  With ferritin: {df['LBXFER'].notna().sum()}")
print(f"  With TFR: {df['LBXTFR'].notna().sum()}")

# Sex labels
df['SEX'] = df['RIAGENDR'].map({1: 'Male', 2: 'Female'})

# Age filter: adults 20+
df = df[df['RIDAGEYR'] >= 20].copy()
print(f"\nAdults 20+: n={len(df)}")
print(f"  Males: {(df['SEX']=='Male').sum()}")
print(f"  Females: {(df['SEX']=='Female').sum()}")

# ── Log-transform PFAS ─────────────────────────────────────────────

pfas_list = ['LBXPFOA', 'LBXPFOS', 'LBXPFHS', 'LBXPFDE', 'LBXPFNA']
pfas_list = [v for v in pfas_list if v in df.columns]

for v in pfas_list:
    df[f'ln_{v}'] = np.log(df[v].clip(lower=0.01))

# ── Analysis function ───────────────────────────────────────────────

def run_regression(df_sub, pfas_var, outcome_var, label='', covariates=None):
    """Weighted linear regression: outcome ~ ln(PFAS) + covariates"""
    ln_pfas = f'ln_{pfas_var}'

    cols = [outcome_var, ln_pfas, 'RIDAGEYR', 'RIDRETH1', 'INDFMPIR']
    if covariates:
        cols += covariates

    tmp = df_sub[cols].dropna()
    if len(tmp) < 30:
        return None

    # Build design matrix
    X_vars = [ln_pfas, 'RIDAGEYR']

    # Race dummies
    for r in [2, 3, 4, 5]:
        tmp[f'race_{r}'] = (tmp['RIDRETH1'] == r).astype(int)
        X_vars.append(f'race_{r}')

    X_vars.append('INDFMPIR')

    if covariates:
        X_vars += covariates

    X = sm.add_constant(tmp[X_vars])
    y = tmp[outcome_var]

    model = sm.OLS(y, X).fit()

    coef = model.params[ln_pfas]
    pval = model.pvalues[ln_pfas]
    ci = model.conf_int().loc[ln_pfas]
    n = len(tmp)

    return {
        'label': label,
        'pfas': pfas_var.replace('LBXPF', '').replace('LBX', ''),
        'outcome': outcome_var,
        'n': n,
        'beta': coef,
        'ci_lo': ci[0],
        'ci_hi': ci[1],
        'p': pval,
        'r2': model.rsquared
    }


# ── Main analysis ──────────────────────────────────────────────────

outcomes = {
    'LBXSAL': 'Albumin (g/dL)',
    'LBXSIR': 'Iron (ug/dL)',
    'LBXFER': 'Ferritin (ng/mL)',
    'LBXTFR': 'Transferrin Receptor (mg/L)'
}

print("\n" + "="*80)
print("ANALYSIS 1: PFAS vs Biomarkers — ALL ADULTS, sex-stratified")
print("Covariates: age, race/ethnicity, poverty-income ratio")
print("="*80)

results = []

for sex in ['Male', 'Female', 'All']:
    df_sub = df if sex == 'All' else df[df['SEX'] == sex]

    for pfas_var in pfas_list:
        for out_var, out_label in outcomes.items():
            r = run_regression(df_sub, pfas_var, out_var, label=f'{sex}')
            if r:
                results.append(r)

# Display results
res_df = pd.DataFrame(results)

# FDR correction within each sex group
from statsmodels.stats.multitest import multipletests

for sex in ['Male', 'Female', 'All']:
    mask = res_df['label'] == sex
    pvals = res_df.loc[mask, 'p'].values
    if len(pvals) > 0:
        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        res_df.loc[mask, 'q'] = qvals

print(f"\n{'Sex':<8} {'PFAS':<6} {'Outcome':<12} {'n':>5} {'beta':>8} {'p':>10} {'q':>10} {'sig'}")
print("-"*75)
for _, r in res_df.sort_values(['label', 'outcome', 'pfas']).iterrows():
    sig = ''
    if r['q'] < 0.05:
        sig = '***'
    elif r['q'] < 0.10:
        sig = '**'
    elif r['p'] < 0.05:
        sig = '*'
    print(f"{r['label']:<8} {r['pfas']:<6} {r['outcome']:<12} {r['n']:>5} {r['beta']:>8.4f} {r['p']:>10.4f} {r['q']:>10.4f} {sig}")


# ── Analysis 2: with albumin adjustment ─────────────────────────────

print("\n" + "="*80)
print("ANALYSIS 2: PFAS vs Iron/Ferritin/TFR — ADJUSTED FOR ALBUMIN")
print("Covariates: age, race, PIR, albumin")
print("="*80)

results2 = []
for sex in ['Male', 'Female']:
    df_sub = df[df['SEX'] == sex]
    for pfas_var in pfas_list:
        for out_var in ['LBXSIR', 'LBXFER', 'LBXTFR']:
            r = run_regression(df_sub, pfas_var, out_var,
                             label=f'{sex}_adj', covariates=['LBXSAL'])
            if r:
                results2.append(r)

res2_df = pd.DataFrame(results2)

for sex_adj in res2_df['label'].unique():
    mask = res2_df['label'] == sex_adj
    pvals = res2_df.loc[mask, 'p'].values
    if len(pvals) > 0:
        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        res2_df.loc[mask, 'q'] = qvals

print(f"\n{'Sex':<12} {'PFAS':<6} {'Outcome':<12} {'n':>5} {'beta':>8} {'p':>10} {'q':>10} {'sig'}")
print("-"*75)
for _, r in res2_df.sort_values(['label', 'outcome', 'pfas']).iterrows():
    sig = ''
    if r['q'] < 0.05:
        sig = '***'
    elif r['q'] < 0.10:
        sig = '**'
    elif r['p'] < 0.05:
        sig = '*'
    print(f"{r['label']:<12} {r['pfas']:<6} {r['outcome']:<12} {r['n']:>5} {r['beta']:>8.4f} {r['p']:>10.4f} {r['q']:>10.4f} {sig}")


# ── Analysis 3: Iron paradox check ──────────────────────────────────

print("\n" + "="*80)
print("ANALYSIS 3: IRON PARADOX — Ferritin and Iron vs PFAS quartiles")
print("="*80)

for sex in ['Male', 'Female']:
    df_sub = df[(df['SEX'] == sex) & df['LBXPFOS'].notna() & df['LBXSIR'].notna()].copy()
    df_sub['PFOS_Q'] = pd.qcut(df_sub['LBXPFOS'], 4, labels=['Q1','Q2','Q3','Q4'])

    print(f"\n{sex} — Iron by PFOS quartile:")
    for q in ['Q1','Q2','Q3','Q4']:
        sub = df_sub[df_sub['PFOS_Q']==q]
        iron_mean = sub['LBXSIR'].mean()
        iron_se = sub['LBXSIR'].std() / np.sqrt(len(sub))
        fer_sub = sub['LBXFER'].dropna()
        fer_mean = fer_sub.mean() if len(fer_sub) > 0 else np.nan
        pfos_range = f"{sub['LBXPFOS'].min():.1f}-{sub['LBXPFOS'].max():.1f}"
        print(f"  {q} (PFOS {pfos_range:>15} ng/mL, n={len(sub):>4}): "
              f"Iron={iron_mean:.1f}+/-{iron_se:.1f}, Ferritin={fer_mean:.1f}")

    # Trend test
    from scipy.stats import spearmanr
    r_iron, p_iron = spearmanr(df_sub['LBXPFOS'], df_sub['LBXSIR'])
    print(f"  Spearman PFOS-Iron: r={r_iron:.3f}, p={p_iron:.4f}")

    fer_sub = df_sub.dropna(subset=['LBXFER'])
    if len(fer_sub) > 30:
        r_fer, p_fer = spearmanr(fer_sub['LBXPFOS'], fer_sub['LBXFER'])
        print(f"  Spearman PFOS-Ferritin: r={r_fer:.3f}, p={p_fer:.4f}")


# ── Analysis 4: Comparison with 2011-2016 albumin effect ──────────

print("\n" + "="*80)
print("ANALYSIS 4: PFOA-ALBUMIN ASSOCIATION — Replication of 2011-2016 finding")
print("="*80)

for sex in ['Male', 'Female']:
    df_sub = df[(df['SEX'] == sex) & df['LBXPFOA'].notna() & df['LBXSAL'].notna()].copy()
    r, p = spearmanr(df_sub['LBXPFOA'], df_sub['LBXSAL'])
    print(f"\n{sex} (n={len(df_sub)}): Spearman PFOA-Albumin r={r:.3f}, p={p:.1e}")

    # Also PFOS
    r2, p2 = spearmanr(df_sub['LBXPFOS'], df_sub['LBXSAL'])
    print(f"{sex} (n={len(df_sub)}): Spearman PFOS-Albumin r={r2:.3f}, p={p2:.1e}")


# ── Summary ─────────────────────────────────────────────────────────

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"""
Cycles E+F (2007-2010), adults 20+, pooled n={len(df)}
Males: {(df['SEX']=='Male').sum()}, Females: {(df['SEX']=='Female').sum()}

NOTE: This analysis is EXPLORATORY. Serum copper was not measured
in NHANES 2007-2010. The primary PFAS-copper hypothesis cannot be
directly tested in these cycles. Results for albumin, iron, and
ferritin serve as hypothesis-generating observations only.

Available biomarkers: albumin, serum iron, ferritin, transferrin receptor
Missing: serum copper, selenium, zinc
""")
