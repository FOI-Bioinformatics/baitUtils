"""
Content checks for the HTML reports and interactive plots produced by
evaluate and compare on the fixture dataset.
"""

import re

import pytest


@pytest.fixture(scope="module")
def evaluate_output(dataset, tmp_path_factory):
    """Run evaluate once per module with the fake pblat."""
    import os
    import stat
    import sys
    from baitUtils.__main__ import main

    root = tmp_path_factory.mktemp("reports")
    bin_dir = root / "bin"
    bin_dir.mkdir()
    script = bin_dir / "pblat"
    script.write_text(f"#!/bin/sh\nif [ $# -eq 0 ]; then exit 0; fi\nfor last; do :; done\ncp '{dataset['psl']}' \"$last\"\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    old_path, old_argv = os.environ["PATH"], sys.argv
    os.environ["PATH"] = f"{bin_dir}{os.pathsep}{old_path}"
    try:
        out = root / "eval"
        sys.argv = ["baitUtils", "evaluate", "-i", str(dataset["baits"]), "-r", str(dataset["reference"]),
                    "-o", str(out), "--offline-plots"]
        main()
        cmp_out = root / "cmp"
        sys.argv = ["baitUtils", "compare", "-r", str(dataset["reference"]), "-o", str(cmp_out),
                    "--sets", f"A:{dataset['baits']}", f"B:{dataset['baits']}"]
        main()
    finally:
        os.environ["PATH"], sys.argv = old_path, old_argv
    return {"eval": out, "cmp": cmp_out}


def _text(html: str) -> str:
    """HTML with tags and scripts removed, whitespace collapsed."""
    html = re.sub(r"<script.*?</script>", " ", html, flags=re.S)
    html = re.sub(r"<style.*?</style>", " ", html, flags=re.S)
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))


class TestEvaluateReport:
    def test_report_states_the_computed_numbers(self, evaluate_output, dataset):
        html = (evaluate_output["eval"] / "coverage_evaluation_report.html").read_text()
        text = _text(html)
        exp = dataset["expected"]
        efficiency = 100.0 * exp["n_mapped"] / exp["n_baits"]
        assert f"{efficiency:.1f}" in text
        breadth = 100.0 * sum(exp["covered"].values()) / 5000
        assert f"{breadth:.1f}" in text
        assert "Total Gaps" in text or "gaps" in text.lower()
        assert re.search(r"\b(Excellent|Good|Fair|Poor)\b", text)

    def test_report_is_self_contained_when_offline(self, evaluate_output):
        html = (evaluate_output["eval"] / "coverage_evaluation_report.html").read_text()
        assert '<script src="https://cdn.plot.ly' not in html
        assert "Plotly" in html
        assert "<style>" in html and "</style>" in html
        assert "font-family" in html

    def test_report_sections_and_plots(self, evaluate_output):
        html = (evaluate_output["eval"] / "coverage_evaluation_report.html").read_text()
        for div_id in ("coverage-pie", "depth-dist", "gap-dist", "quality-radar"):
            assert f'id="{div_id}"' in html or f"'{div_id}'" in html or div_id in html
        text = _text(html)
        assert "Reference" in text
        assert "Recommendations" in text or "recommend" in text.lower()

    def test_interactive_plots_are_plotly_documents(self, evaluate_output):
        plots = sorted((evaluate_output["eval"] / "interactive_plots").glob("*.html"))
        assert [p.name for p in plots] == ["coverage_dashboard.html", "detailed_gap_analysis.html",
                                           "quality_assessment.html", "reference_analysis.html"]
        for path in plots:
            content = path.read_text()
            assert "plotly" in content.lower()
            assert "<div" in content

    def test_static_plots_exist(self, evaluate_output):
        names = {p.name for p in (evaluate_output["eval"] / "plots").glob("*.png")}
        assert {"coverage_overview.png", "depth_distribution.png", "gap_analysis.png",
                "per_reference_coverage.png", "coverage_heatmap.png"} <= names


class TestCompareReport:
    def test_report_lists_sets_ranking_and_statistics(self, evaluate_output):
        html = (evaluate_output["cmp"] / "comparative_analysis_report.html").read_text()
        text = _text(html)
        assert " A " in text and " B " in text
        assert "Statistical Analysis" in text
        assert "adjusted for multiple comparisons (fdr)" in text
        assert "Not applicable" in text       # identical sets, paired tests cannot run
        assert "P-value: nan" not in text
        assert "Kolmogorov-Smirnov" in text and "Levene" in text
        assert "Best-hit identity per bait" in text

    def test_comparison_outputs(self, evaluate_output):
        cmp_out = evaluate_output["cmp"]
        assert (cmp_out / "comparison_matrix.csv").read_text().count("\n") >= 3
        plots = {p.name for p in (cmp_out / "comparative_plots").glob("*.png")}
        assert {"coverage_distributions.png", "gap_analysis_comparison.png",
                "performance_heatmap.png", "statistical_comparison.png"} <= plots
        interactive = {p.name for p in (cmp_out / "interactive_comparative_plots").glob("*.html")}
        assert {"comparison_dashboard.html", "performance_radar.html"} <= interactive
