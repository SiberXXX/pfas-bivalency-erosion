import sys, io, gzip, os
import numpy as np
from scipy import stats
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# =====================================================
# ПРОВЕРКА: обогащение бивалентных — реальное или
# артефакт gene length bias?
# =====================================================

# Загрузка
ensembl_to_sym = {}
with open('../geo_GSE254408/ensembl_to_symbol.tsv', 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            ensembl_to_sym[parts[0]] = parts[1]

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

gene_ids = []
gene_lengths = []
count_rows = []

with gzip.open('../geo_GSE254408/GSE254408_all_samples_raw_counts.txt.gz', 'rt') as f:
    f.readline()
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        gid = parts[0]
        sym = ensembl_to_sym.get(gid, '')
        if not sym:
            continue
        gene_ids.append(sym)
        gene_lengths.append(int(parts[5]))
        count_rows.append([int(parts[i]) for i in range(6, 27)])

count_matrix = np.array(count_rows, dtype=np.float64)
gene_lengths = np.array(gene_lengths)
n_genes = len(gene_ids)

# CPM
lib_sizes = count_matrix.sum(axis=0)
cpm = count_matrix / lib_sizes * 1e6

# Groups
ctrl_idx = [2, 9, 16]
pfoa_idx = [5, 12, 19]

ctrl_cpm = cpm[:, ctrl_idx]
pfoa_cpm = cpm[:, pfoa_idx]
ctrl_mean = ctrl_cpm.mean(axis=1)
pfoa_mean = pfoa_cpm.mean(axis=1)
expressed = (ctrl_mean > 1) | (pfoa_mean > 1)

is_biv = np.array([g.upper() in bivalent for g in gene_ids])

print(f'Всего генов с символами: {n_genes}')
print(f'Экспрессированные (CPM>1): {expressed.sum()}')
print(f'Бивалентные среди экспрессированных: {(is_biv & expressed).sum()}')

# =====================================================
# ТЕСТ 1: Длина генов у DEGs vs не-DEGs
# =====================================================
print(f'\n{"=" * 60}')
print('ТЕСТ 1: Длина DEGs vs не-DEGs (PFOA)')
print(f'{"=" * 60}')

from statsmodels.stats.multitest import multipletests

pvals_all = []
log2fcs = []
gene_idx_tested = []

for i in range(n_genes):
    if not expressed[i]:
        continue
    c = ctrl_cpm[i, :]
    h = pfoa_cpm[i, :]
    if c.std() == 0 and h.std() == 0:
        continue
    t, p = stats.ttest_ind(h, c, equal_var=False)
    if np.isnan(p):
        continue
    fc = (pfoa_mean[i] + 0.1) / (ctrl_mean[i] + 0.1)
    log2fc = np.log2(fc) if fc > 0 else 0
    pvals_all.append(p)
    log2fcs.append(log2fc)
    gene_idx_tested.append(i)

pvals_all = np.array(pvals_all)
log2fcs = np.array(log2fcs)
gene_idx_tested = np.array(gene_idx_tested)

is_deg = pvals_all < 0.05
is_deg_down = is_deg & (log2fcs < 0)

deg_lengths = gene_lengths[gene_idx_tested[is_deg]]
nondeg_lengths = gene_lengths[gene_idx_tested[~is_deg]]
deg_down_lengths = gene_lengths[gene_idx_tested[is_deg_down]]

print(f'DEG (nom p<0.05) median length: {np.median(deg_lengths):.0f}')
print(f'Non-DEG median length: {np.median(nondeg_lengths):.0f}')
print(f'DEG DOWN median length: {np.median(deg_down_lengths):.0f}')
print(f'Отношение DEG/nonDEG: {np.median(deg_lengths)/np.median(nondeg_lengths):.2f}x')

u, up = stats.mannwhitneyu(deg_lengths, nondeg_lengths, alternative='greater')
print(f'Mann-Whitney DEG > nonDEG: p={up:.4e}')

# =====================================================
# ТЕСТ 2: Length-matched permutation test
# Для каждого DEG находим nonDEG с похожей длиной
# и проверяем бивалентное обогащение
# =====================================================
print(f'\n{"=" * 60}')
print('ТЕСТ 2: Обогащение с контролем на длину (length-stratified)')
print(f'{"=" * 60}')

# Стратификация по квинтилям длины
tested_lengths = gene_lengths[gene_idx_tested]
tested_is_biv = is_biv[gene_idx_tested]

# Для DOWN DEGs
quintiles = np.percentile(tested_lengths, [20, 40, 60, 80])

def stratified_enrichment(is_target, tested_lengths, tested_is_biv, label):
    """Вычисляет обогащение с контролем на длину через стратификацию"""
    bins = np.digitize(tested_lengths, quintiles)

    observed_biv = 0
    expected_biv = 0
    total_target = 0

    for b in range(len(quintiles) + 1):
        in_bin = bins == b
        target_in_bin = is_target & in_bin
        n_target = target_in_bin.sum()
        n_bin = in_bin.sum()

        if n_target == 0 or n_bin == 0:
            continue

        total_target += n_target
        observed_biv += (target_in_bin & tested_is_biv).sum()

        # Ожидаемая доля бивалентных в этом бине (среди ВСЕХ генов бина)
        bin_biv_rate = tested_is_biv[in_bin].sum() / n_bin
        expected_biv += n_target * bin_biv_rate

    obs_rate = observed_biv / total_target if total_target > 0 else 0
    exp_rate = expected_biv / total_target if total_target > 0 else 0
    enrichment = obs_rate / exp_rate if exp_rate > 0 else 0

    print(f'{label}:')
    print(f'  Target genes: {total_target}')
    print(f'  Observed bivalent: {observed_biv} ({100*obs_rate:.1f}%)')
    print(f'  Expected (length-matched): {expected_biv:.1f} ({100*exp_rate:.1f}%)')
    print(f'  Enrichment (length-corrected): {enrichment:.2f}x')

    # Биномиальный тест
    if total_target > 0:
        binom_p = stats.binomtest(observed_biv, total_target, exp_rate, alternative='greater').pvalue
        print(f'  Binomial test p={binom_p:.4e}')

    return enrichment, observed_biv, expected_biv

# Все DEGs
stratified_enrichment(is_deg, tested_lengths, tested_is_biv, 'Все DEGs (nom p<0.05)')

# DOWN DEGs
stratified_enrichment(is_deg_down, tested_lengths, tested_is_biv, 'DOWN DEGs (nom p<0.05)')

# UP DEGs
is_deg_up = is_deg & (log2fcs > 0)
stratified_enrichment(is_deg_up, tested_lengths, tested_is_biv, 'UP DEGs (nom p<0.05)')

# =====================================================
# ТЕСТ 3: Permutation test (рандомизация)
# Берём N случайных генов с такой же длиной, считаем
# долю бивалентных, повторяем 10000 раз
# =====================================================
print(f'\n{"=" * 60}')
print('ТЕСТ 3: Permutation test с контролем на длину')
print(f'{"=" * 60}')

np.random.seed(42)

def length_matched_permutation(target_mask, tested_lengths, tested_is_biv, n_perm=10000, label=''):
    target_idx = np.where(target_mask)[0]
    n_target = len(target_idx)
    observed_biv = tested_is_biv[target_idx].sum()

    # Для каждого target гена, найдём генов с похожей длиной (±30%)
    all_idx = np.arange(len(tested_lengths))

    # Стратифицированный permutation:
    # разбиваем по децилям длины, сэмплируем в каждом
    decile_bins = np.digitize(tested_lengths, np.percentile(tested_lengths, np.arange(10, 100, 10)))
    target_bins = decile_bins[target_idx]

    perm_biv_counts = np.zeros(n_perm)
    for p in range(n_perm):
        sampled = []
        for b in range(10):
            n_in_bin = (target_bins == b).sum()
            pool = all_idx[decile_bins == b]
            if len(pool) > 0 and n_in_bin > 0:
                chosen = np.random.choice(pool, size=n_in_bin, replace=True)
                sampled.extend(chosen)
        perm_biv_counts[p] = tested_is_biv[sampled].sum()

    perm_mean = perm_biv_counts.mean()
    perm_p = (perm_biv_counts >= observed_biv).sum() / n_perm

    print(f'{label}:')
    print(f'  n={n_target}, observed bivalent={observed_biv}')
    print(f'  Permutation mean={perm_mean:.1f}, sd={perm_biv_counts.std():.1f}')
    print(f'  Permutation p-value: {perm_p:.4f} ({n_perm} permutations)')
    print(f'  Z-score: {(observed_biv - perm_mean) / max(perm_biv_counts.std(), 0.01):.2f}')

    return perm_p

# PFOA DOWN
length_matched_permutation(is_deg_down, tested_lengths, tested_is_biv, 10000, 'PFOA DOWN (nom p<0.05)')

# PFOA ALL
length_matched_permutation(is_deg, tested_lengths, tested_is_biv, 10000, 'PFOA ALL (nom p<0.05)')

# PFOA UP
length_matched_permutation(is_deg_up, tested_lengths, tested_is_biv, 10000, 'PFOA UP (nom p<0.05)')

# =====================================================
# ТЕСТ 4: То же для GSE301375 (мышь)
# =====================================================
print(f'\n{"=" * 60}')
print('ТЕСТ 4: GSE301375 (мышь) — контроль на длину')
print(f'{"=" * 60}')

# У нас нет длин для мышиных генов напрямую
# Но мы можем проверить: в GSE301375, DEGs при nominal p<0.05
# имеют медианную длину больше чем non-DEGs?
# К сожалению, формат GSE301375 — простые count файлы без длины
# Мы можем только проверить по expression level

# Загрузим DEG-лист и проверим среднюю экспрессию
degs_301375 = []
with open(os.path.join('..', 'geo_GSE301375', 'GSE301375_degs_nominal.txt'), 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            degs_301375.append(parts[0].upper())

degs_301375_set = set(degs_301375)
mouse_degs_biv = degs_301375_set & bivalent
print(f'GSE301375 DEGs (nom): {len(degs_301375_set)}')
print(f'Бивалентные (by name match): {len(mouse_degs_biv)}')
print(f'ПРЕДУПРЕЖДЕНИЕ: без контроля на длину для мышиных данных')
print(f'Длина генов НЕ доступна в GSE301375 формате')
print(f'Bias от длины ВЕРОЯТЕН и для этого датасета')

print(f'\n{"=" * 60}')
print('ИТОГ')
print(f'{"=" * 60}')
print(f'Бивалентные гены в 3.3x длиннее чем не-бивалентные')
print(f'DEGs в RNA-seq имеют bias к длинным генам')
print(f'→ Обогащение МОЖЕТ быть артефактом length bias')
print(f'→ Нужно смотреть на length-corrected p-value')
