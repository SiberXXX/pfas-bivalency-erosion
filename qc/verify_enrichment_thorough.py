import sys, io, gzip, os
import numpy as np
from scipy import stats
from collections import defaultdict, Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# =====================================================
# ТЩАТЕЛЬНАЯ ПРОВЕРКА КАЖДОГО ЗВЕНА
# =====================================================

print('=' * 70)
print('ПРОВЕРКА 1: Дубликаты в Ensembl→Symbol маппинге')
print('=' * 70)

ensembl_to_sym = {}
with open('../geo_GSE254408/ensembl_to_symbol.tsv', 'r') as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            ensembl_to_sym[parts[0]] = parts[1]

# Какие символы имеют >1 Ensembl ID в нашем файле каунтов?
gene_ids_raw = []
with gzip.open('../geo_GSE254408/GSE254408_all_samples_raw_counts.txt.gz', 'rt') as f:
    f.readline()
    f.readline()
    for line in f:
        gene_ids_raw.append(line.split('\t')[0])

sym_to_ensgs_in_data = defaultdict(list)
for gid in gene_ids_raw:
    if gid in ensembl_to_sym:
        sym_to_ensgs_in_data[ensembl_to_sym[gid]].append(gid)

dups = {s: ids for s, ids in sym_to_ensgs_in_data.items() if len(ids) > 1}
print(f'Символов с >1 Ensembl ID в данных: {len(dups)}')
print(f'Из них с >5 ID: {sum(1 for s, ids in dups.items() if len(ids) > 5)}')

# ПРОБЛЕМА: если один символ встречается от нескольких ENSG,
# то при t-тесте каждый ENSG тестируется ОТДЕЛЬНО.
# Один и тот же ген может попасть в DEGs несколько раз
# (или попасть И в DEG И в non-DEG одновременно!)

# Сколько дубликатных символов среди бивалентных?
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

dup_biv = {s for s in dups if s.upper() in bivalent}
print(f'Дубликатных символов среди бивалентных: {len(dup_biv)}')
print(f'Примеры: {sorted(list(dup_biv))[:10]}')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 2: Правильный DEG анализ БЕЗ дубликатов')
print('Стратегия: для дублей берём ENSG с максимальной экспрессией')
print('=' * 70)

# Загрузка данных
all_data = []  # (ensembl_id, symbol, length, counts[21])
with gzip.open('../geo_GSE254408/GSE254408_all_samples_raw_counts.txt.gz', 'rt') as f:
    f.readline()
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        gid = parts[0]
        sym = ensembl_to_sym.get(gid, '')
        if not sym:
            continue
        length = int(parts[5])
        counts = [int(parts[i]) for i in range(6, 27)]
        all_data.append((gid, sym, length, counts))

print(f'Всего записей с символами: {len(all_data)}')

# Дедупликация: для каждого символа берём ENSG с максимальным суммарным каунтом
best_per_sym = {}
for gid, sym, length, counts in all_data:
    total = sum(counts)
    key = sym.upper()
    if key not in best_per_sym or total > best_per_sym[key][3]:
        best_per_sym[key] = (gid, sym, length, total, counts)

print(f'Уникальных символов после дедупликации: {len(best_per_sym)}')
deduped_before = len(all_data)
deduped_after = len(best_per_sym)
print(f'Удалено дубликатов: {deduped_before - deduped_after}')

# Собираем дедуплицированный массив
gene_symbols_dedup = []
gene_lengths_dedup = []
count_matrix_dedup = []

for key in sorted(best_per_sym.keys()):
    gid, sym, length, total, counts = best_per_sym[key]
    gene_symbols_dedup.append(sym)
    gene_lengths_dedup.append(length)
    count_matrix_dedup.append(counts)

count_matrix = np.array(count_matrix_dedup, dtype=np.float64)
gene_lengths = np.array(gene_lengths_dedup)
n_genes = len(gene_symbols_dedup)

# CPM
lib_sizes = count_matrix.sum(axis=0)
cpm = count_matrix / lib_sizes * 1e6

# Groups (0-based in the 21-sample block)
ctrl_idx = [2, 9, 16]
pfoa_idx = [5, 12, 19]

ctrl_cpm = cpm[:, ctrl_idx]
pfoa_cpm = cpm[:, pfoa_idx]
ctrl_mean = ctrl_cpm.mean(axis=1)
pfoa_mean = pfoa_cpm.mean(axis=1)

# Фильтр экспрессии
expressed = (ctrl_mean > 1) | (pfoa_mean > 1)
is_biv = np.array([g.upper() in bivalent for g in gene_symbols_dedup])

print(f'Экспрессированные: {expressed.sum()}')
print(f'Бивалентные среди экспрессированных: {(is_biv & expressed).sum()}')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 3: DEG анализ (дедуплицированный)')
print('=' * 70)

from statsmodels.stats.multitest import multipletests

pvals = []
log2fcs = []
tested_idx = []

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
    pvals.append(p)
    log2fcs.append(log2fc)
    tested_idx.append(i)

pvals = np.array(pvals)
log2fcs = np.array(log2fcs)
tested_idx = np.array(tested_idx)

n_tested = len(pvals)
_, padj, _, _ = multipletests(pvals, method='fdr_bh')

is_deg = pvals < 0.05
is_deg_down = is_deg & (log2fcs < 0)
is_deg_up = is_deg & (log2fcs > 0)

print(f'Протестировано (дедупл.): {n_tested}')
print(f'DEGs nom p<0.05: {is_deg.sum()} (UP:{is_deg_up.sum()}, DOWN:{is_deg_down.sum()})')
print(f'DEGs FDR<0.05: {(padj < 0.05).sum()}')

# Базовые проценты
tested_biv = is_biv[tested_idx]
tested_lengths_arr = gene_lengths[tested_idx]
bg_biv_rate = tested_biv.sum() / n_tested
print(f'Background bivalent rate: {tested_biv.sum()}/{n_tested} = {100*bg_biv_rate:.1f}%')

# Наивное обогащение (без коррекции на длину)
for label, mask in [('ALL DEGs', is_deg), ('DOWN DEGs', is_deg_down), ('UP DEGs', is_deg_up)]:
    n = mask.sum()
    biv_n = (mask & tested_biv).sum()
    rate = biv_n / n if n > 0 else 0
    naive_enrich = rate / bg_biv_rate if bg_biv_rate > 0 else 0
    # Fisher
    a = biv_n
    b = n - a
    c = tested_biv.sum() - a
    d = n_tested - tested_biv.sum() - b
    ort, fp = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    print(f'  {label:12s}: {n:5d} genes, {biv_n:3d} biv ({100*rate:.1f}%) naive={naive_enrich:.2f}x Fisher p={fp:.4e}')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 4: Длина генов — подробная диагностика')
print('=' * 70)

biv_lens = tested_lengths_arr[tested_biv]
nonbiv_lens = tested_lengths_arr[~tested_biv]
deg_lens = tested_lengths_arr[is_deg]
nondeg_lens = tested_lengths_arr[~is_deg]
deg_down_lens = tested_lengths_arr[is_deg_down]
deg_up_lens = tested_lengths_arr[is_deg_up]

print(f'Длина (медиана):')
print(f'  Бивалентные:     {np.median(biv_lens):8.0f}')
print(f'  Не-бивалентные:  {np.median(nonbiv_lens):8.0f}')
print(f'  DEG all:         {np.median(deg_lens):8.0f}')
print(f'  DEG DOWN:        {np.median(deg_down_lens):8.0f}')
print(f'  DEG UP:          {np.median(deg_up_lens):8.0f}')
print(f'  NOT DEG:         {np.median(nondeg_lens):8.0f}')

# Ключевой вопрос: DEG DOWN короче чем DEG UP?
print(f'\n  DEG DOWN vs DEG UP длина: {np.median(deg_down_lens):.0f} vs {np.median(deg_up_lens):.0f}')
print(f'  DEG DOWN vs NOT DEG длина: {np.median(deg_down_lens):.0f} vs {np.median(nondeg_lens):.0f}')

u1, p1 = stats.mannwhitneyu(deg_down_lens, nondeg_lens)
print(f'  Mann-Whitney DOWN vs nonDEG: p={p1:.4e} (two-sided)')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 5: Length-stratified enrichment (20 бинов)')
print('=' * 70)

n_bins = 20
bin_edges = np.percentile(tested_lengths_arr, np.linspace(0, 100, n_bins + 1)[1:-1])
length_bins = np.digitize(tested_lengths_arr, bin_edges)

def length_corrected_enrichment(target_mask, label):
    obs_biv = 0
    exp_biv = 0.0
    total = 0

    for b in range(n_bins):
        in_bin = length_bins == b
        target_in_bin = target_mask & in_bin
        n_t = target_in_bin.sum()
        n_b = in_bin.sum()

        if n_t == 0 or n_b == 0:
            continue

        total += n_t
        obs_biv += (target_in_bin & tested_biv).sum()
        bin_biv_rate = tested_biv[in_bin].sum() / n_b
        exp_biv += n_t * bin_biv_rate

    if total == 0:
        print(f'  {label}: нет данных')
        return

    obs_rate = obs_biv / total
    exp_rate = exp_biv / total
    enrich = obs_rate / exp_rate if exp_rate > 0 else 0

    # Binomial test
    bp = stats.binomtest(obs_biv, total, exp_rate, alternative='greater').pvalue

    print(f'  {label}:')
    print(f'    n={total}, obs_biv={obs_biv} ({100*obs_rate:.1f}%), '
          f'exp_biv={exp_biv:.1f} ({100*exp_rate:.1f}%)')
    print(f'    Enrichment (length-corrected): {enrich:.3f}x')
    print(f'    Binomial p={bp:.6e}')

    return enrich, bp

length_corrected_enrichment(is_deg, 'ALL DEGs')
length_corrected_enrichment(is_deg_down, 'DOWN DEGs')
length_corrected_enrichment(is_deg_up, 'UP DEGs')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 6: Permutation test (10 бинов, 50000 итераций)')
print('=' * 70)

np.random.seed(42)
N_PERM = 50000

decile_edges = np.percentile(tested_lengths_arr, np.arange(10, 100, 10))
decile_bins = np.digitize(tested_lengths_arr, decile_edges)

def permutation_test(target_mask, label):
    target_indices = np.where(target_mask)[0]
    n_target = len(target_indices)
    obs_biv = tested_biv[target_indices].sum()

    # Считаем сколько target генов в каждом бине
    target_bin_counts = np.bincount(decile_bins[target_indices], minlength=10)

    # Для каждого бина, получаем пул индексов
    bin_pools = []
    bin_biv_flags = []
    for b in range(10):
        pool = np.where(decile_bins == b)[0]
        bin_pools.append(pool)
        bin_biv_flags.append(tested_biv[pool])

    perm_counts = np.zeros(N_PERM)
    for p in range(N_PERM):
        total_biv = 0
        for b in range(10):
            n_need = target_bin_counts[b]
            if n_need == 0 or len(bin_pools[b]) == 0:
                continue
            chosen = np.random.choice(len(bin_pools[b]), size=n_need, replace=True)
            total_biv += bin_biv_flags[b][chosen].sum()
        perm_counts[p] = total_biv

    perm_mean = perm_counts.mean()
    perm_std = perm_counts.std()
    p_value = (perm_counts >= obs_biv).sum() / N_PERM

    print(f'  {label}:')
    print(f'    n={n_target}, observed bivalent={obs_biv}')
    print(f'    Permutation: mean={perm_mean:.1f}, std={perm_std:.1f}')
    print(f'    Z-score: {(obs_biv - perm_mean) / max(perm_std, 0.01):.2f}')
    print(f'    p-value: {p_value:.6f} ({N_PERM} permutations)')

    return p_value

permutation_test(is_deg_down, 'PFOA DOWN')
permutation_test(is_deg_up, 'PFOA UP')
permutation_test(is_deg, 'PFOA ALL')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 7: А может DOWN-DEGs просто КОРОЧЕ и поэтому')
print('бивалентные среди них не выбиваются по длине?')
print('=' * 70)

# Посмотрим: среди DOWN DEGs, различается ли длина
# между бивалентными и не-бивалентными?
down_idx = tested_idx[is_deg_down]
down_biv_mask = is_biv[down_idx]
down_biv_lens = gene_lengths[down_idx[down_biv_mask]]
down_nonbiv_lens = gene_lengths[down_idx[~down_biv_mask]]

print(f'DOWN DEGs: {len(down_idx)} genes')
print(f'  Бивалентные: n={down_biv_mask.sum()}, median length={np.median(down_biv_lens):.0f}')
print(f'  Не-бивалентные: n={(~down_biv_mask).sum()}, median length={np.median(down_nonbiv_lens):.0f}')
print(f'  Отношение длин: {np.median(down_biv_lens)/np.median(down_nonbiv_lens):.2f}x')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 8: Множественное сравнение по PFAS')
print('Мы тестируем 6 PFAS → Bonferroni = 6')
print('=' * 70)

all_pfas = {
    'PFOA': [5, 12, 19],
    'PFOS': [6, 13, 20],
    'PFDA': [3, 10, 17],
    'PFDS': [4, 11, 18],
    'FTOH': [0, 7, 14],
    'FTS': [1, 8, 15],
}

down_pvals = []
for pfas_name, pfas_idx in all_pfas.items():
    p_cpm = cpm[:, pfas_idx]
    p_mean = p_cpm.mean(axis=1)
    expr = (ctrl_mean > 1) | (p_mean > 1)

    pv = []
    lf = []
    ti = []
    for i in range(n_genes):
        if not expr[i]:
            continue
        c = ctrl_cpm[i, :]
        h = p_cpm[i, :]
        if c.std() == 0 and h.std() == 0:
            continue
        t, p = stats.ttest_ind(h, c, equal_var=False)
        if np.isnan(p):
            continue
        fc = (p_mean[i] + 0.1) / (ctrl_mean[i] + 0.1)
        log2fc = np.log2(fc) if fc > 0 else 0
        pv.append(p)
        lf.append(log2fc)
        ti.append(i)

    pv = np.array(pv)
    lf = np.array(lf)
    ti = np.array(ti)
    is_d_down = (pv < 0.05) & (lf < 0)
    ti_biv = is_biv[ti]

    # Length-corrected permutation (quick, 10000)
    if is_d_down.sum() > 0:
        target_i = np.where(is_d_down)[0]
        obs = ti_biv[target_i].sum()
        dec_bins_local = np.digitize(gene_lengths[ti], decile_edges)
        t_bin_counts = np.bincount(dec_bins_local[target_i], minlength=10)
        perm_c = np.zeros(10000)
        for pp in range(10000):
            tb = 0
            for b in range(10):
                nn = t_bin_counts[b]
                pool = np.where(dec_bins_local == b)[0]
                if nn == 0 or len(pool) == 0:
                    continue
                ch = np.random.choice(pool, size=nn, replace=True)
                tb += ti_biv[ch].sum()
            perm_c[pp] = tb
        perm_p = (perm_c >= obs).sum() / 10000
    else:
        obs = 0
        perm_p = 1.0

    n_down = is_d_down.sum()
    print(f'  {pfas_name:6s}: DOWN={n_down:4d}, biv={obs}, perm_p={perm_p:.4f}')
    down_pvals.append(perm_p)

# Bonferroni correction
print(f'\nBonferroni коррекция (×6):')
for i, (pfas_name, _) in enumerate(all_pfas.items()):
    corrected = min(down_pvals[i] * 6, 1.0)
    sig = '***' if corrected < 0.001 else '**' if corrected < 0.01 else '*' if corrected < 0.05 else 'NS'
    print(f'  {pfas_name:6s}: raw p={down_pvals[i]:.4f}, Bonferroni p={corrected:.4f} {sig}')

# =====================================================
print('\n' + '=' * 70)
print('ПРОВЕРКА 9: Комбинированный анализ (все PFAS) —')
print('те же тесты (length-corrected)')
print('=' * 70)

all_pfas_idx = []
for v in all_pfas.values():
    all_pfas_idx.extend(v)

all_pfas_cpm = cpm[:, all_pfas_idx]
all_pfas_mean = all_pfas_cpm.mean(axis=1)
expr_comb = (ctrl_mean > 1) | (all_pfas_mean > 1)

pv_c = []
lf_c = []
ti_c = []
for i in range(n_genes):
    if not expr_comb[i]:
        continue
    c = ctrl_cpm[i, :]
    h = all_pfas_cpm[i, :]
    if c.std() == 0 and h.std() == 0:
        continue
    t, p = stats.ttest_ind(h, c, equal_var=False)
    if np.isnan(p):
        continue
    fc = (all_pfas_mean[i] + 0.1) / (ctrl_mean[i] + 0.1)
    log2fc = np.log2(fc) if fc > 0 else 0
    pv_c.append(p)
    lf_c.append(log2fc)
    ti_c.append(i)

pv_c = np.array(pv_c)
lf_c = np.array(lf_c)
ti_c = np.array(ti_c)
_, padj_c, _, _ = multipletests(pv_c, method='fdr_bh')

n_fdr05 = (padj_c < 0.05).sum()
is_down_fdr = (padj_c < 0.05) & (lf_c < 0)
is_up_fdr = (padj_c < 0.05) & (lf_c > 0)
is_down_nom = (pv_c < 0.05) & (lf_c < 0)

print(f'Combined: {len(pv_c)} tested, {n_fdr05} FDR<0.05')
print(f'DOWN FDR<0.05: {is_down_fdr.sum()}, UP FDR<0.05: {is_up_fdr.sum()}')

# Length-corrected для combined DOWN FDR<0.05
ti_c_biv = is_biv[ti_c]
ti_c_lens = gene_lengths[ti_c]
dec_bins_c = np.digitize(ti_c_lens, decile_edges)

for label, mask in [('DOWN FDR<0.05', is_down_fdr), ('UP FDR<0.05', is_up_fdr),
                     ('DOWN nom', is_down_nom)]:
    if mask.sum() == 0:
        print(f'  {label}: 0 genes')
        continue
    target_i = np.where(mask)[0]
    obs = ti_c_biv[target_i].sum()
    n_t = len(target_i)
    t_bin_counts = np.bincount(dec_bins_c[target_i], minlength=10)
    perm_c2 = np.zeros(10000)
    for pp in range(10000):
        tb = 0
        for b in range(10):
            nn = t_bin_counts[b]
            pool = np.where(dec_bins_c == b)[0]
            if nn == 0 or len(pool) == 0:
                continue
            ch = np.random.choice(pool, size=nn, replace=True)
            tb += ti_c_biv[ch].sum()
        perm_c2[pp] = tb
    perm_p2 = (perm_c2 >= obs).sum() / 10000
    obs_rate = obs / n_t
    exp_rate = perm_c2.mean() / n_t
    enrich = obs_rate / exp_rate if exp_rate > 0 else 0
    print(f'  {label}: n={n_t}, biv={obs} ({100*obs_rate:.1f}%), '
          f'expected={perm_c2.mean():.1f} ({100*exp_rate:.1f}%), '
          f'enrichment={enrich:.2f}x, perm_p={perm_p2:.4f}')

print('\n' + '=' * 70)
print('ФИНАЛЬНЫЙ ВЕРДИКТ')
print('=' * 70)
