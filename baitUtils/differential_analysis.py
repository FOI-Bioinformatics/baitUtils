#!/usr/bin/env python3

"""
differential_analysis.py

Statistical comparison of oligo sets evaluated against the same reference.

All tests run on observed data:

- Coverage distributions are compared on mean depth per non-overlapping
  window (default 1000 bp) taken from the per-base coverage arrays. Windows
  are used rather than bases because neighbouring bases under one bait are
  not independent observations. Both sets cover the same windows, so the
  location test is paired (Wilcoxon signed-rank on the per-window
  differences); Kolmogorov-Smirnov and Levene compare the distributions'
  shape and spread and treat the windows as independent.
- Per-reference metrics (breadth, mean depth, gap count) are compared with
  paired tests across references (Wilcoxon signed-rank for two sets,
  Friedman for more). With fewer than two references the test is reported
  as not applicable rather than given a p-value.
- Gap sizes are compared with a Mann-Whitney U test.

P-values within each family of tests are adjusted for multiple comparisons
with the configured method and stored in p_adjusted.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import scipy.stats as stats

from baitUtils.comparative_analyzer import OligoSetResult


@dataclass
class StatisticalTest:
    """Result of one statistical test."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: float
    significance_level: str        # '***', '**', '*', 'ns' or 'n/a'
    interpretation: str
    p_adjusted: Optional[float] = None
    applicable: bool = True
    n: int = 0
    what: str = ""      # quantity compared, used to word the interpretation
    detail: str = ""    # statistic summary shown in the interpretation

    @property
    def p_reported(self) -> float:
        """Adjusted p-value when available, otherwise the raw p-value."""
        return self.p_value if self.p_adjusted is None else self.p_adjusted


@dataclass
class CoverageDistributionComparison:
    """Comparison of window-level coverage depth between two oligo sets."""
    set1_name: str
    set2_name: str
    ks_test: StatisticalTest
    wilcoxon_test: StatisticalTest      # paired, per window
    levene_test: StatisticalTest
    summary_stats: Dict[str, Dict[str, float]]
    window_size: int = 1000


def _not_applicable(test_name: str, reason: str, n: int = 0) -> StatisticalTest:
    return StatisticalTest(test_name, float('nan'), float('nan'), float('nan'),
                           'n/a', f"Not applicable: {reason}", None, False, n)


def _matched_rank_biserial(diff: np.ndarray) -> float:
    """
    Matched-pairs rank-biserial correlation of paired differences, in -1..1:
    (sum of ranks of positive differences - sum of ranks of negative ones) /
    total rank sum, zero differences dropped as in the Wilcoxon test.
    Positive when the first set is larger. Defined also when every
    difference is the same, unlike a standardised mean difference.
    """
    diff = np.asarray(diff, dtype=float)
    diff = diff[diff != 0]
    if diff.size == 0:
        return float('nan')
    ranks = stats.rankdata(np.abs(diff))
    return float((ranks[diff > 0].sum() - ranks[diff < 0].sum()) / ranks.sum())


def _all_tied(*samples) -> bool:
    """True when every value in every sample is the same (rank tests are undefined)."""
    values = np.concatenate([np.asarray(s, dtype=float) for s in samples])
    return values.size == 0 or bool(np.allclose(values, values[0]))


class DifferentialAnalyzer:
    """Statistical comparison of oligo set performance."""

    METRIC_LABELS = {
        'coverage_breadth': 'coverage breadth',
        'mean_depth': 'mean coverage depth',
        'gap_count': 'gap count',
    }

    def __init__(self, significance_level: float = 0.05, correction_method: str = 'fdr',
                 window_size: int = 1000):
        """
        Args:
            significance_level: Alpha used to call a difference significant
            correction_method: 'bonferroni', 'holm', 'fdr' or 'none'
            window_size: Window length in bp for coverage distribution tests
        """
        if correction_method not in ('bonferroni', 'holm', 'fdr', 'none'):
            raise ValueError(f"Unknown correction method: {correction_method}")
        self.significance_level = significance_level
        self.correction_method = correction_method
        self.window_size = window_size
        # Star labels scale with alpha: '*' at alpha, '**' at alpha/5, '***' at alpha/50
        # (0.05, 0.01 and 0.001 for the default alpha of 0.05)
        self.alpha_levels = {significance_level / 50: '***', significance_level / 5: '**',
                             significance_level: '*'}

    # ------------------------------------------------------------------ data
    def _window_means(self, oligo_set: OligoSetResult, references: Optional[List[str]] = None) -> np.ndarray:
        """
        Mean depth per non-overlapping window, concatenated over the given
        references in that order (default: all, sorted by name), so that two
        sets give aligned windows.
        """
        arrays = oligo_set.coverage_arrays or {}
        if not arrays:
            raise ValueError(
                f"Oligo set '{oligo_set.name}' has no coverage arrays; "
                "coverage distribution tests need per-base coverage")
        means = []
        for ref in (sorted(arrays) if references is None else references):
            arr = np.asarray(arrays[ref], dtype=float)
            n_full = len(arr) // self.window_size
            if n_full:
                means.append(arr[:n_full * self.window_size].reshape(n_full, self.window_size).mean(axis=1))
            if len(arr) % self.window_size:
                means.append(np.array([arr[n_full * self.window_size:].mean()]))
        return np.concatenate(means) if means else np.array([])

    # --------------------------------------------------- coverage distribution
    def compare_coverage_distributions(self, set1: OligoSetResult,
                                       set2: OligoSetResult) -> CoverageDistributionComparison:
        """
        Compare window-level coverage depth between two oligo sets. Windows
        are paired by reference and position (the references both sets
        share), so the location test is a Wilcoxon signed-rank test on the
        per-window differences; KS and Levene treat windows as independent.
        """
        shared = sorted(set(set1.coverage_arrays or {}) & set(set2.coverage_arrays or {}))
        if not shared:
            # raises the usual error for a set without coverage arrays
            self._window_means(set1), self._window_means(set2)
        w1 = self._window_means(set1, shared)
        w2 = self._window_means(set2, shared)
        n1, n2 = len(w1), len(w2)

        if n1 < 2 or n2 < 2 or n1 != n2:
            reason = (f"only {min(n1, n2)} window(s) of {self.window_size} bp" if min(n1, n2) < 2
                      else "the sets' references differ in length")
            ks = _not_applicable("Kolmogorov-Smirnov", reason, n1 + n2)
            wx = _not_applicable("Wilcoxon signed-rank (paired windows)", reason, n1)
            lev = _not_applicable("Levene", reason, n1 + n2)
        else:
            ks_stat, ks_p = stats.ks_2samp(w1, w2)
            ks = self._test("Kolmogorov-Smirnov", ks_stat, ks_p, ks_stat,
                            "window depth distributions (windows treated as independent)",
                            f"D={ks_stat:.3f}", n1 + n2)

            diff = w1 - w2
            nonzero = int(np.count_nonzero(diff))
            if nonzero == 0:
                wx = _not_applicable("Wilcoxon signed-rank (paired windows)",
                                     "identical depth in every window", n1)
            else:
                wx_stat, wx_p = stats.wilcoxon(w1, w2)
                effect = _matched_rank_biserial(diff)     # positive when set1 is deeper
                wx = self._test("Wilcoxon signed-rank (paired windows)", wx_stat, wx_p, effect,
                                "window depth (paired by window)",
                                f"W={wx_stat:g}, {nonzero} of {n1} windows differ", nonzero)

            if np.allclose(w1, w1[0]) and np.allclose(w2, w2[0]):
                lev = _not_applicable("Levene", "no variance in either set", n1 + n2)
            else:
                lev_stat, lev_p = stats.levene(w1, w2)
                var_ratio = (np.var(w1, ddof=1) / np.var(w2, ddof=1)) if np.var(w2, ddof=1) > 0 else float('inf')
                lev = self._test("Levene", lev_stat, lev_p, var_ratio,
                                 "window depth variability (windows treated as independent)",
                                 f"W={lev_stat:.3f}", n1 + n2)

        self._adjust([ks, wx, lev])

        def summary(w: np.ndarray) -> Dict[str, float]:
            return {
                'mean': float(np.mean(w)), 'median': float(np.median(w)), 'std': float(np.std(w, ddof=1)) if len(w) > 1 else 0.0,
                'q25': float(np.percentile(w, 25)), 'q75': float(np.percentile(w, 75)),
                'windows': int(len(w)),
            }

        return CoverageDistributionComparison(
            set1_name=set1.name, set2_name=set2.name, ks_test=ks, wilcoxon_test=wx,
            levene_test=lev, summary_stats={set1.name: summary(w1), set2.name: summary(w2)},
            window_size=self.window_size)

    # ------------------------------------------------------ per-reference tests
    def compare_quality_metrics(self, oligo_sets: List[OligoSetResult]) -> Dict[str, StatisticalTest]:
        """
        Paired comparison of per-reference metrics across oligo sets.

        Each reference sequence is one paired observation. Two sets are
        compared with the Wilcoxon signed-rank test, more with the Friedman
        test. Metrics with fewer than two references, or with no difference
        between sets, are reported as not applicable.
        """
        if len(oligo_sets) < 2:
            raise ValueError("Need at least 2 oligo sets for comparison")

        per_ref = [s.coverage_stats.get('per_reference', {}) for s in oligo_sets]
        refs = sorted(set.intersection(*(set(p) for p in per_ref))) if per_ref else []
        # Gap counts per reference use the gaps of the gap analysis (at or above
        # min_gap_size), the same definition as total_gaps and the gap-size test.
        # Gap lists without coordinates (synthetic data) fall back to the
        # per-reference run counts of the coverage statistics.
        located = all('chromosome' in g for s in oligo_sets for g in s.gap_analysis.get('gaps', []))
        gap_counts = [self._gaps_per_reference(s) for s in oligo_sets] if located else None
        source_keys = {'coverage_breadth': 'coverage_breadth', 'mean_depth': 'mean_depth',
                       'gap_count': None if located else 'gaps'}

        results: Dict[str, StatisticalTest] = {}
        for metric, key in source_keys.items():
            label = self.METRIC_LABELS[metric]
            if len(refs) < 2:
                results[metric] = _not_applicable(
                    "Wilcoxon signed-rank" if len(oligo_sets) == 2 else "Friedman",
                    f"{len(refs)} shared reference(s); paired tests need at least 2", len(refs))
                continue

            if key is None:
                matrix = np.array([[float(counts.get(r, 0)) for r in refs] for counts in gap_counts])
            else:
                matrix = np.array([[float(p[r].get(key, 0.0)) for r in refs] for p in per_ref])
            if len(oligo_sets) == 2:
                diff = matrix[0] - matrix[1]
                if np.all(diff == 0):
                    results[metric] = _not_applicable("Wilcoxon signed-rank",
                                                      f"identical {label} on every reference", len(refs))
                    continue
                stat, p = stats.wilcoxon(matrix[0], matrix[1])
                effect = _matched_rank_biserial(diff)
                # scipy drops zero differences, so n is the references that differ
                nonzero = int(np.count_nonzero(diff))
                results[metric] = self._test("Wilcoxon signed-rank", stat, p, effect, label,
                                             f"W={stat:g}, {nonzero} of {len(refs)} references differ",
                                             nonzero)
            else:
                if np.all(matrix == matrix[0]):
                    results[metric] = _not_applicable("Friedman", f"identical {label} on every reference", len(refs))
                    continue
                stat, p = stats.friedmanchisquare(*matrix)
                k, n = matrix.shape
                kendall_w = float(stat / (n * (k - 1))) if n * (k - 1) > 0 else float('nan')
                results[metric] = self._test("Friedman", stat, p, kendall_w,
                                             label, f"chi2={stat:.2f}, n={n}", n)

        self._adjust(list(results.values()))
        return results

    # -------------------------------------------------------------- gap sizes
    def analyze_gap_patterns(self, set1: OligoSetResult, set2: OligoSetResult) -> Dict[str, StatisticalTest]:
        """Compare gap size distributions between two oligo sets."""
        sizes1 = [g['length'] for g in set1.gap_analysis.get('gaps', [])]
        sizes2 = [g['length'] for g in set2.gap_analysis.get('gaps', [])]
        if len(sizes1) < 2 or len(sizes2) < 2:
            return {'gap_size_distribution': _not_applicable(
                "Mann-Whitney U (gap sizes)", "fewer than two gaps in a set", len(sizes1) + len(sizes2))}
        if _all_tied(sizes1, sizes2):
            return {'gap_size_distribution': _not_applicable(
                "Mann-Whitney U (gap sizes)", "every gap has the same length", len(sizes1) + len(sizes2))}
        stat, p = stats.mannwhitneyu(sizes1, sizes2, alternative='two-sided')
        effect = 2.0 * float(stat) / (len(sizes1) * len(sizes2)) - 1.0
        result = self._test("Mann-Whitney U (gap sizes)", stat, p, effect,
                            "gap size distributions", f"U={stat:g}", len(sizes1) + len(sizes2))
        self._adjust([result])
        return {'gap_size_distribution': result}

    # ------------------------------------------------------- per-oligo identity
    def compare_oligo_identity(self, set1: OligoSetResult, set2: OligoSetResult) -> StatisticalTest:
        """
        Mann-Whitney U test on best-hit identity per mapped bait. Baits are
        independent observations, so this is the most direct comparison of
        design quality between two sets.
        """
        def identities(oligo_set):
            table = oligo_set.per_oligo
            if table is None or len(table) == 0:
                return np.array([])
            # baits that --max-hits turned into unmapped baits are not mapped baits
            if 'multi_mapping_excluded' in table.columns:
                table = table[~table['multi_mapping_excluded'].astype(bool)]
            return table['best_identity'].to_numpy(dtype=float)

        a, b = identities(set1), identities(set2)
        if len(a) < 2 or len(b) < 2:
            return _not_applicable("Mann-Whitney U (bait identity)", "fewer than two mapped baits in a set",
                                   len(a) + len(b))
        if _all_tied(a, b):
            return _not_applicable("Mann-Whitney U (bait identity)", "identical identity for every bait",
                                   len(a) + len(b))
        stat, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        effect = 2.0 * float(stat) / (len(a) * len(b)) - 1.0
        result = self._test("Mann-Whitney U (bait identity)", stat, p, effect,
                            "best-hit identity per bait", f"U={stat:g}", len(a) + len(b))
        self._adjust([result])
        return result

    # ----------------------------------------------------- multiple comparison
    def multiple_comparison_correction(self, p_values: List[float], method: Optional[str] = None) -> List[float]:
        """Adjust p-values with 'bonferroni', 'holm', 'fdr' (Benjamini-Hochberg) or 'none'."""
        method = method or self.correction_method
        p = np.asarray(p_values, dtype=float)
        if p.size == 0:
            return []
        if method == 'none':
            return p.tolist()
        if method == 'bonferroni':
            return np.minimum(p * p.size, 1.0).tolist()
        if method == 'holm':
            order = np.argsort(p)
            adjusted = np.empty_like(p)
            running = 0.0
            for rank, idx in enumerate(order):
                running = max(running, min(1.0, p[idx] * (p.size - rank)))
                adjusted[idx] = running
            return adjusted.tolist()
        if method == 'fdr':
            order = np.argsort(p)
            adjusted = np.empty_like(p)
            running = 1.0
            for rank in range(p.size - 1, -1, -1):
                idx = order[rank]
                running = min(running, p[idx] * p.size / (rank + 1))
                adjusted[idx] = min(1.0, running)
            return adjusted.tolist()
        raise ValueError(f"Unknown correction method: {method}")

    def _adjust(self, tests: List[StatisticalTest]) -> None:
        """
        Fill p_adjusted for the applicable tests in one family. Stars and the
        interpretation are rewritten from the adjusted p-value so that the
        three never disagree.
        """
        live = [t for t in tests if t.applicable]
        if not live:
            return
        for test, p_adj in zip(live, self.multiple_comparison_correction([t.p_value for t in live])):
            test.p_adjusted = p_adj
            test.significance_level = self._get_significance_level(p_adj)
            test.interpretation = self._interpret(p_adj, test.what, test.detail)

    @staticmethod
    def _gaps_per_reference(oligo_set: OligoSetResult) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for gap in oligo_set.gap_analysis.get('gaps', []):
            chrom = gap.get('chromosome')
            if chrom is not None:
                counts[chrom] = counts.get(chrom, 0) + 1
        return counts

    def _test(self, test_name: str, statistic, p_value, effect_size, what: str, detail: str,
              n: int) -> StatisticalTest:
        """Build an applicable test result; p_adjusted is filled by _adjust."""
        return StatisticalTest(test_name, float(statistic), float(p_value), float(effect_size),
                               self._get_significance_level(p_value),
                               self._interpret(p_value, what, detail), None, True, n, what, detail)

    # ----------------------------------------------------------------- helpers
    def _get_significance_level(self, p_value: float) -> str:
        """Star label for a p-value relative to the configured alpha; 'ns' above alpha."""
        if p_value is None or np.isnan(p_value):
            return 'n/a'
        for threshold, level in sorted(self.alpha_levels.items()):
            if p_value <= threshold:
                return level
        return 'ns'

    def significance_legend(self) -> str:
        """Text explaining the star labels for the configured alpha."""
        a = self.significance_level
        return (f"*** p<={a / 50:g}, ** p<={a / 5:g}, * p<={a:g}, ns p>{a:g} "
                f"(adjusted p-values, {self.correction_method})")

    def _interpret(self, p_value: float, what: str, detail: str) -> str:
        if p_value <= self.significance_level:
            return f"Difference in {what} at alpha {self.significance_level} ({detail})"
        return f"No difference in {what} at alpha {self.significance_level} ({detail})"
