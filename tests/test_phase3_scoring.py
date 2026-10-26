"""
Tests for phase 3 of the 2026-10 code audit: the quality scorer's
uniformity target, weight normalisation and reference difficulty credit,
the consolidated recommendations, the benchmark without a Gini
coefficient, the generated methodology text and the seaborn calls.
"""

import warnings

import pandas as pd
import pytest

from baitUtils.quality_scorer import QualityScorer
from baitUtils.benchmark import BenchmarkAnalyzer
from baitUtils.evaluate import generate_recommendations
from baitUtils.report_generator import InteractiveReportGenerator
from baitUtils.plotting_utils import StatisticalPlotter


STATS = {'coverage_breadth': 80.0, 'mean_depth': 10.0, 'coverage_gini': 0.3, 'coverage_cv': 0.9,
         'mapping_efficiency': 90.0}
GAPS = {'gap_percentage': 5.0, 'max_gap_size': 800, 'total_gaps': 12}


def _summary(total, extreme_gc=0, repeats=0, challenging=0):
    return {'analysis_summary': {'total_sequences': total, 'sequences_with_extreme_gc': extreme_gc,
                                 'sequences_with_high_repeats': repeats, 'challenging_sequences': challenging}}


class TestUniformity:
    def test_uniformity_follows_the_gini_target(self):
        """F-19: the depth_uniformity target was never used."""
        even = QualityScorer({**STATS, 'coverage_gini': 0.1}, GAPS)._score_coverage_depth()
        uneven = QualityScorer({**STATS, 'coverage_gini': 0.7}, GAPS)._score_coverage_depth()
        assert even > uneven
        # 1 - 0.2 = 0.8 reaches the target exactly; the depth half is 10/20 -> 0.5
        at_target = QualityScorer({**STATS, 'coverage_gini': 0.2}, GAPS)._score_coverage_depth()
        assert at_target == pytest.approx((0.5 + 1.0) / 2)
        # the target can be overridden
        strict = QualityScorer({**STATS, 'coverage_gini': 0.2}, GAPS, benchmarks={'depth_uniformity': 0.9})
        assert strict._score_coverage_depth() == pytest.approx((0.5 + 0.8 / 0.9) / 2)

    def test_cv_is_only_a_fallback(self):
        without_gini = {k: v for k, v in STATS.items() if k != 'coverage_gini'}
        score = QualityScorer({**without_gini, 'coverage_cv': 0.5}, GAPS)._score_coverage_depth()
        assert score == pytest.approx((0.5 + 1.0) / 2)


class TestWeights:
    def test_overrides_are_normalised(self):
        """F-21: weights summing to more than 1 produced scores above 1."""
        scorer = QualityScorer(STATS, GAPS, weights={'coverage_breadth': 1.3})
        assert sum(scorer.weights.values()) == pytest.approx(1.0)
        assert 0.0 <= scorer.calculate_quality_score().overall_score <= 1.0

    def test_defaults_are_untouched(self):
        assert QualityScorer(STATS, GAPS).weights == QualityScorer.DEFAULT_WEIGHTS

    def test_non_positive_weights_are_rejected(self):
        with pytest.raises(ValueError):
            QualityScorer(STATS, GAPS, weights={k: 0.0 for k in QualityScorer.DEFAULT_WEIGHTS})


class TestReferenceDifficulty:
    def test_same_breadth_scores_higher_on_a_harder_reference(self):
        """F-20: the docstring promised credit; the code penalised."""
        easy = QualityScorer(STATS, GAPS, _summary(10))
        hard = QualityScorer(STATS, GAPS, _summary(10, extreme_gc=10, repeats=10, challenging=10))
        assert easy.reference_difficulty_index() == 0.0
        assert hard.reference_difficulty_index() == 1.0
        assert hard.expected_breadth() == pytest.approx(QualityScorer.DEFAULT_BENCHMARKS['coverage_breadth'] / 2)
        assert hard._score_reference_difficulty() > easy._score_reference_difficulty()
        assert hard._score_reference_difficulty() == 1.0        # 80 % on the hardest reference
        assert easy._score_reference_difficulty() == pytest.approx(80.0 / 95.0)

    def test_no_reference_analysis_counts_as_easy(self):
        """F-45: no more undocumented 0.7 placeholder."""
        none = QualityScorer(STATS, GAPS)._score_reference_difficulty()
        easy = QualityScorer(STATS, GAPS, _summary(5))._score_reference_difficulty()
        assert none == easy == pytest.approx(80.0 / 95.0)

    def test_component_stays_in_range(self):
        scorer = QualityScorer({**STATS, 'coverage_breadth': 100.0}, GAPS, _summary(2, challenging=2))
        assert scorer._score_reference_difficulty() == 1.0
        scorer = QualityScorer({**STATS, 'coverage_breadth': 0.0}, GAPS, _summary(2, challenging=2))
        assert scorer._score_reference_difficulty() == 0.0


class TestRecommendations:
    def test_text_report_uses_the_scorer_and_gap_suggestions(self):
        """F-43: one source of recommendations for the text and HTML reports."""
        gaps = {**GAPS, 'gap_percentage': 40.0, 'suggestions': ['Use baitUtils fill to close the gaps.']}
        quality = QualityScorer({**STATS, 'coverage_breadth': 40.0}, gaps).calculate_quality_score()
        recs = generate_recommendations({**STATS, 'coverage_breadth': 40.0}, gaps, quality)
        assert recs[0] == 'Use baitUtils fill to close the gaps.'
        assert recs[1:] == [r for r in quality.recommendations if r not in recs[:1]]
        # computed on the fly when no quality score is passed
        assert generate_recommendations({**STATS, 'coverage_breadth': 40.0}, gaps) == recs


class TestBenchmarkWithoutGini:
    def test_missing_gini_is_not_applicable(self):
        """F-47: a missing Gini used to become a fabricated 0.5."""
        result = BenchmarkAnalyzer({'coverage_breadth': 90.0}, {'gap_percentage': 1.0}).benchmark_depth_uniformity()
        assert result.category == "Not applicable"
        assert result.actual_score != result.actual_score
        analyzer = BenchmarkAnalyzer({'coverage_breadth': 90.0}, {'gap_percentage': 1.0})
        text = analyzer.generate_benchmark_report(analyzer.run_full_benchmark())
        assert "Not applicable" in text and "0.50" not in text.split("Depth Uniformity")[1].split("Gap Reduction")[0]


class TestMethodology:
    def test_methodology_names_the_mapper_and_the_constants(self, tmp_path):
        """F-40"""
        gen = InteractiveReportGenerator(
            coverage_stats={'parameters': {'min_identity': 92.0, 'min_length': 80, 'min_coverage': 2.0}},
            gap_analysis={}, quality_scores={}, reference_analysis={}, output_dir=tmp_path, mapper="minimap2")
        gen._create_methodology()
        html = gen.html_components[-1]
        assert "mapped to the reference with minimap2" in html
        assert "at least 92.0%" in html and "at least 80 bp" in html
        assert "Excellent &gt;= 0.85" in html and "Coverage Breadth</strong>: 30%" in html
        assert "pblat" not in html


class TestSeabornCalls:
    def test_category_plots_do_not_use_deprecated_palette_without_hue(self, tmp_path):
        """F-38: seaborn 0.14 removes palette without hue."""
        df = pd.DataFrame({'gc_content': [40.0, 45.0, 50.0, 55.0, 60.0, 62.0],
                           'kept': pd.Categorical([True, False, True, False, True, False])})
        plotter = StatisticalPlotter()
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            plotter.plot_statistics(df=df, columns=['gc_content'], plot_types=['boxplot', 'violinplot', 'boxenplot'],
                                    color_col='kept', outdir=str(tmp_path), file_format='png')
        assert (tmp_path / "boxplot_gc_content.png").exists()
