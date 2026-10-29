#!/usr/bin/env python3

"""
coverage_viz.py

Static figures of the evaluate command (matplotlib), one concern per figure:

- coverage_profile: depth along each reference with gaps shaded
- coverage_overview: base composition per reference and the cumulative depth curve
- depth_distribution: depth histogram and depth spread per reference
- gap_analysis: gap length distribution and the largest gaps
- per_reference_coverage: breadth per reference and breadth against length
- coverage_heatmap: a coverage map, references by relative position
- quality_score: the five component scores and the overall score

Every figure uses the design system in baitUtils.design (one palette, one
typography, recessive axes) and says when it has nothing to show instead
of drawing empty axes.
"""

import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from baitUtils import design
from baitUtils.design import CHROME, SEQUENTIAL, style_axes, empty_panel
from baitUtils.quality_scorer import QualityScorer

warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

ACCENT = CHROME["accent"]
GAP = CHROME["gap_fill"]
MAX_PROFILE_REFERENCES = 8
MAX_BAR_REFERENCES = 25


def bp_formatter():
    return FuncFormatter(lambda value, _pos: design.fmt_bp(value))


class CoverageVisualizer:
    """Static figures for the coverage evaluation."""

    def __init__(
        self,
        coverage_stats: Dict[str, Any],
        gap_analysis: Dict[str, Any],
        output_dir: Path,
        plot_format: str = 'png',
        dpi: int = 300,
        coverage_arrays: Optional[Dict[str, np.ndarray]] = None,
        quality_scores: Optional[Dict[str, Any]] = None,
    ):
        """
        Args:
            coverage_stats: Output of CoverageAnalyzer.analyze
            gap_analysis: Output of GapAnalyzer.analyze
            output_dir: Directory that receives the plots/ folder
            plot_format: png, pdf or svg
            dpi: Resolution of raster output
            coverage_arrays: Per-reference depth arrays (for the profile, the
                depth spread and the coverage map); figures that need them are
                skipped when absent
            quality_scores: QualityScore.to_dict() for the quality figure
        """
        self.coverage_stats = coverage_stats
        self.gap_analysis = gap_analysis
        self.output_dir = Path(output_dir)
        self.plot_format = plot_format
        self.dpi = dpi
        self.coverage_arrays = coverage_arrays or {}
        self.quality_scores = quality_scores or {}

        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        design.apply_matplotlib_style(dpi=dpi)

    # ------------------------------------------------------------- driver
    def generate_all_plots(self) -> Dict[str, str]:
        """Generate every figure; returns {name: path} of the files written."""
        written: Dict[str, str] = {}
        steps = [
            ("coverage_profile", self.plot_coverage_profile),
            ("coverage_overview", self.plot_coverage_overview),
            ("depth_distribution", self.plot_depth_distribution),
            ("gap_analysis", self.plot_gap_analysis),
            ("per_reference_coverage", self.plot_per_reference_coverage),
            ("coverage_heatmap", self.plot_coverage_heatmap),
            ("quality_score", self.plot_quality_score),
        ]
        for name, method in steps:
            logging.info(f"Generating {name} figure...")
            path = method()
            if path:
                written[name] = str(path)
        logging.info(f"{len(written)} figures written to {self.plots_dir}")
        return written

    # ----------------------------------------------------------- helpers
    @property
    def _per_reference(self) -> Dict[str, Dict[str, Any]]:
        return self.coverage_stats.get('per_reference', {}) or {}

    @property
    def _min_coverage(self) -> float:
        return float(self.coverage_stats.get('parameters', {}).get('min_coverage', 1.0))

    @property
    def _target_coverage(self) -> float:
        return float(self.coverage_stats.get('parameters', {}).get('target_coverage', 10.0))

    def _gaps_by_reference(self) -> Dict[str, List[Dict[str, Any]]]:
        gaps: Dict[str, List[Dict[str, Any]]] = {}
        for gap in self.gap_analysis.get('gaps', []) or []:
            gaps.setdefault(gap.get('chromosome', ''), []).append(gap)
        return gaps

    def _save(self, fig, name: str) -> Optional[Path]:
        path = self.plots_dir / f"{name}.{self.plot_format}"
        try:
            fig.savefig(path, dpi=self.dpi, bbox_inches='tight', facecolor=CHROME["surface"], edgecolor='none')
            return path
        except Exception as e:
            logging.error(f"Error saving figure {name}: {e}")
            return None
        finally:
            plt.close(fig)

    # ----------------------------------------------------------- figures
    def plot_coverage_profile(self) -> Optional[Path]:
        """Depth along each of the longest references, gaps shaded."""
        if not self.coverage_arrays:
            logging.info("No coverage arrays; the coverage profile is skipped")
            return None
        refs = design.select_references(
            {r: {'length': len(a)} for r, a in self.coverage_arrays.items()}, MAX_PROFILE_REFERENCES)
        gaps = self._gaps_by_reference()
        n = len(refs)
        fig, axes = plt.subplots(n, 1, figsize=(11, 1.6 * n + 1.2), sharex=False, squeeze=False)
        max_depth = max(float(np.max(self.coverage_arrays[r])) for r in refs) if refs else 1.0
        ymax = max(max_depth, self._min_coverage) * 1.15 or 1.0
        for ax, ref in zip(axes[:, 0], refs):
            array = self.coverage_arrays[ref]
            x, y = design.bin_coverage(array)
            ax.fill_between(x, 0, y, color=ACCENT, alpha=0.12, linewidth=0)
            ax.plot(x, y, color=ACCENT, linewidth=1.4)
            for gap in gaps.get(ref, []):
                ax.axvspan(gap['start'], gap['end'], color=GAP, alpha=0.18, linewidth=0)
            ax.axhline(self._min_coverage, color=CHROME["ink_muted"], linewidth=0.8, linestyle=(0, (3, 3)))
            ax.set_xlim(0, len(array))
            ax.set_ylim(0, ymax)
            ax.xaxis.set_major_formatter(bp_formatter())
            ax.yaxis.set_major_locator(MaxNLocator(3, integer=True))
            style_axes(ax, grid="y")
            stats = self._per_reference.get(ref, {})
            label = (f"{design.truncate(ref, 32)}   breadth {stats.get('coverage_breadth', 0):.1f}%, "
                     f"mean depth {design.fmt_depth(stats.get('mean_depth', float(np.mean(array))))}, "
                     f"{len(gaps.get(ref, []))} gap(s)")
            ax.set_title(label, fontsize=9.5, loc="left")
        axes[-1, 0].set_xlabel("Position along reference")
        fig.supylabel("Depth (aligned baits per base)", fontsize=9.5, color=CHROME["ink_secondary"])
        total = len(self.coverage_arrays)
        title = "Coverage profile"
        if total > n:
            title += f" of the {n} longest of {total} references"
        fig.suptitle(title, x=0.01, ha="left", fontsize=12)
        fig.text(0.01, 0.005,
                 f"Mean depth per bin; shaded: gaps at or above the minimum gap size; "
                 f"dotted line: minimum depth {self._min_coverage:g}",
                 fontsize=8, color=CHROME["ink_muted"])
        fig.tight_layout(rect=(0.02, 0.02, 1, 0.97))
        return self._save(fig, "coverage_profile")

    def plot_coverage_overview(self) -> Optional[Path]:
        """Base composition (covered, uncovered, N) per reference and the cumulative depth curve."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.2, 1]})

        per_ref = self._per_reference
        rows = []
        total = self.coverage_stats
        if total.get('total_bases'):
            rows.append(("All references", total.get('covered_bases', 0), total.get('n_bases', 0), total['total_bases']))
        for ref in design.select_references(per_ref, 12):
            s = per_ref[ref]
            rows.append((design.truncate(ref, 22), s.get('covered_bases', 0), s.get('n_bases', 0), s.get('length', 0)))
        if rows:
            labels = [r[0] for r in rows][::-1]
            covered = np.array([r[1] / max(r[3], 1) * 100 for r in rows])[::-1]
            n_part = np.array([r[2] / max(r[3], 1) * 100 for r in rows])[::-1]
            uncovered = np.clip(100 - covered - n_part, 0, 100)
            y = np.arange(len(rows))
            height = 0.62
            ax1.barh(y, covered, height=height, color=ACCENT, label="Covered")
            ax1.barh(y, uncovered, left=covered, height=height, color=GAP, alpha=0.75, label="Uncovered")
            ax1.barh(y, n_part, left=covered + uncovered, height=height, color=CHROME["axis"], label="N bases")
            ax1.set_yticks(y)
            ax1.set_yticklabels(labels)
            ax1.set_xlim(0, 100)
            ax1.set_xlabel("Share of bases (%)")
            ax1.set_title("Coverage breadth by reference")
            ax1.legend(loc="lower right", ncol=3, bbox_to_anchor=(1, 1.0), frameon=False, fontsize=8.5)
            for yi, value in zip(y, covered):
                if value >= 12:
                    ax1.text(value - 1, yi, f"{value:.1f}%", ha="right", va="center", fontsize=8, color="#ffffff")
                else:
                    ax1.text(value + 1, yi, f"{value:.1f}%", ha="left", va="center", fontsize=8, color=CHROME["ink_secondary"])
            if len(per_ref) > 12:
                ax1.text(0, -0.9, f"{len(per_ref) - 12} shorter references not shown", fontsize=8, color=CHROME["ink_muted"])
            style_axes(ax1, grid="x")
            ax1.spines["left"].set_visible(False)
        else:
            empty_panel(ax1, "No per-reference statistics")

        self._draw_cumulative_depth(ax2)
        fig.tight_layout()
        return self._save(fig, "coverage_overview")

    def _draw_cumulative_depth(self, ax) -> None:
        """Fraction of bases at or above each depth, from the arrays when available."""
        if self.coverage_arrays:
            values = np.concatenate([np.asarray(a) for a in self.coverage_arrays.values()])
            if values.size:
                upper = int(max(np.percentile(values, 99.5), self._target_coverage, 1))
                depths = np.arange(0, upper + 2)
                fraction = [100.0 * np.mean(values >= d) for d in depths]
                ax.step(depths, fraction, where="post", color=ACCENT, linewidth=2)
                ax.fill_between(depths, fraction, step="post", color=ACCENT, alpha=0.1, linewidth=0)
                ax.set_xlim(0, upper + 1)
        else:
            dist = self.coverage_stats.get('depth_distribution', {}) or {}
            if not dist:
                empty_panel(ax, "No depth distribution")
                return
            thresholds = sorted(dist)
            ax.plot(thresholds, [dist[t] for t in thresholds], color=ACCENT, marker="o", markersize=5)
        ax.axvline(self._min_coverage, color=CHROME["ink_muted"], linewidth=0.8, linestyle=(0, (3, 3)))
        ax.text(self._min_coverage, 101, f" min {self._min_coverage:g}x", fontsize=8, color=CHROME["ink_muted"], va="bottom")
        if self._target_coverage > self._min_coverage:
            ax.axvline(self._target_coverage, color=CHROME["ink_secondary"], linewidth=0.8, linestyle=(0, (3, 3)))
            ax.text(self._target_coverage, 101, f" target {self._target_coverage:g}x", fontsize=8,
                    color=CHROME["ink_secondary"], va="bottom")
        ax.set_ylim(0, 112)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xlabel("Depth d")
        ax.set_ylabel("Bases with depth ≥ d (%)")
        ax.set_title("Cumulative depth")
        style_axes(ax, grid="y")

    def plot_depth_distribution(self) -> Optional[Path]:
        """Depth histogram and the depth spread per reference."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 1.1]})

        hist = self.coverage_stats.get('depth_histogram', {}) or {}
        if self.coverage_arrays:
            values = np.concatenate([np.asarray(a) for a in self.coverage_arrays.values()])
            upper = int(max(np.percentile(values, 99.5), 1)) + 1
            bins = np.arange(0, upper + 1) - 0.5
            counts, edges = np.histogram(values, bins=bins)
            centres = (edges[:-1] + edges[1:]) / 2
            ax1.bar(centres, counts / max(values.size, 1) * 100, width=0.8, color=ACCENT)
            ax1.set_xlabel("Depth (aligned baits per base)")
            ax1.set_ylabel("Bases (%)")
            ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        elif hist.get('counts'):
            counts = np.asarray(hist['counts'], dtype=float)
            edges = np.asarray(hist['bin_edges'], dtype=float)
            ax1.bar((edges[:-1] + edges[1:]) / 2, counts / counts.sum() * 100, width=np.diff(edges) * 0.9, color=ACCENT)
            ax1.set_xlabel("Depth")
            ax1.set_ylabel("Bases (%)")
        else:
            empty_panel(ax1, "No depth histogram")
        ax1.set_title("Depth distribution")
        style_axes(ax1, grid="y")

        per_ref = self._per_reference
        refs = design.select_references(per_ref, MAX_BAR_REFERENCES)
        if refs and self.coverage_arrays:
            rows = []
            for ref in refs:
                array = np.asarray(self.coverage_arrays.get(ref, []))
                if array.size:
                    rows.append((ref, np.percentile(array, 10), np.median(array), np.percentile(array, 90), float(np.mean(array))))
            if rows:
                y = np.arange(len(rows))[::-1]
                for yi, (ref, p10, med, p90, mean) in zip(y, rows):
                    ax2.plot([p10, p90], [yi, yi], color=SEQUENTIAL[3], linewidth=3, solid_capstyle="round")
                    ax2.plot(med, yi, "o", color=ACCENT, markersize=6, markeredgecolor=CHROME["surface"], markeredgewidth=1.5)
                ax2.set_yticks(y)
                ax2.set_yticklabels([design.truncate(r[0], 22) for r in rows])
                ax2.set_xlabel("Depth: median (dot) and 10th to 90th percentile (bar)")
                ax2.set_title("Depth per reference")
                ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
                style_axes(ax2, grid="x")
                ax2.spines["left"].set_visible(False)
        elif refs:
            depths = [per_ref[r].get('mean_depth', 0) for r in refs]
            y = np.arange(len(refs))[::-1]
            ax2.barh(y, depths, height=0.6, color=ACCENT)
            ax2.set_yticks(y)
            ax2.set_yticklabels([design.truncate(r, 22) for r in refs])
            ax2.set_xlabel("Mean depth")
            ax2.set_title("Mean depth per reference")
            style_axes(ax2, grid="x")
            ax2.spines["left"].set_visible(False)
        else:
            empty_panel(ax2, "No per-reference statistics")
        fig.tight_layout()
        return self._save(fig, "depth_distribution")

    def plot_gap_analysis(self) -> Optional[Path]:
        """Gap length distribution and the largest gaps."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 1.2]})
        gaps = self.gap_analysis.get('gaps', []) or []
        lengths = np.array([g['length'] for g in gaps], dtype=float)

        if lengths.size:
            lo, hi = max(lengths.min(), 1), max(lengths.max(), 2)
            if hi / lo > 10 and lengths.size >= 5:
                bins = np.logspace(np.log10(lo), np.log10(hi * 1.001), min(20, max(5, int(np.sqrt(lengths.size)) + 3)))
                ax1.set_xscale("log")
            else:
                spec = design.histogram_bins(lengths, min(15, max(3, int(np.sqrt(lengths.size)) + 2)))
                bins = np.arange(spec["start"], spec["end"] + spec["size"] / 2, spec["size"])
            ax1.hist(lengths, bins=bins, color=ACCENT, rwidth=0.92)
            ax1.set_xlabel("Gap length (bp)")
            ax1.set_ylabel("Gaps")
            ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax1.set_title(f"Gap lengths ({lengths.size:,} gaps, {design.fmt_bp(lengths.sum())} in total)")
            style_axes(ax1, grid="y")
        else:
            empty_panel(ax1, "No gaps at or above the minimum gap size")
            ax1.set_title("Gap lengths")

        largest = sorted(gaps, key=lambda g: g['length'], reverse=True)[:12]
        if largest:
            y = np.arange(len(largest))[::-1]
            sizes = [g['length'] for g in largest]
            ax2.barh(y, sizes, height=0.6, color=GAP, alpha=0.85)
            ax2.set_yticks(y)
            ax2.set_yticklabels([f"{design.truncate(g['chromosome'], 18)}:{g['start']:,}-{g['end']:,}" for g in largest],
                                fontsize=8)
            for yi, size in zip(y, sizes):
                ax2.text(size + max(sizes) * 0.01, yi, design.fmt_bp(size), va="center", fontsize=8,
                         color=CHROME["ink_secondary"])
            ax2.set_xlim(0, max(sizes) * 1.18)
            ax2.set_ylim(-0.7, max(len(largest), 4) - 0.3)
            ax2.set_xlabel("Gap length")
            ax2.xaxis.set_major_formatter(bp_formatter())
            ax2.set_title("Largest gaps")
            style_axes(ax2, grid="x")
            ax2.spines["left"].set_visible(False)
        else:
            empty_panel(ax2, "No gaps to rank")
            ax2.set_title("Largest gaps")
        fig.tight_layout()
        return self._save(fig, "gap_analysis")

    def plot_per_reference_coverage(self) -> Optional[Path]:
        """Breadth per reference and breadth against reference length."""
        per_ref = self._per_reference
        if not per_ref:
            logging.warning("No per-reference statistics available")
            return None
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.2, 1]})

        refs = design.select_references(per_ref, MAX_BAR_REFERENCES)
        breadth = [per_ref[r].get('coverage_breadth', 0) for r in refs]
        y = np.arange(len(refs))[::-1]
        ax1.barh(y, breadth, height=0.6, color=ACCENT)
        overall = self.coverage_stats.get('coverage_breadth')
        if overall is not None:
            ax1.axvline(overall, color=CHROME["ink_secondary"], linewidth=0.8, linestyle=(0, (3, 3)))
            ax1.text(overall, len(refs) - 0.35, f" overall {overall:.1f}%", fontsize=8, va="bottom",
                     ha="right" if overall > 75 else "left", color=CHROME["ink_secondary"])
            ax1.set_ylim(-0.7, len(refs) + 0.1)
        ax1.set_yticks(y)
        ax1.set_yticklabels([design.truncate(r, 22) for r in refs])
        ax1.set_xlim(0, 100)
        ax1.set_xlabel("Coverage breadth (%)")
        title = "Breadth per reference"
        if len(per_ref) > len(refs):
            title += f" ({len(refs)} longest of {len(per_ref)})"
        ax1.set_title(title)
        style_axes(ax1, grid="x")
        ax1.spines["left"].set_visible(False)

        lengths = np.array([s.get('length', 0) for s in per_ref.values()], dtype=float)
        breadths = np.array([s.get('coverage_breadth', 0) for s in per_ref.values()], dtype=float)
        depths = np.array([s.get('mean_depth', 0) for s in per_ref.values()], dtype=float)
        sizes = 20 + 60 * (depths / depths.max() if depths.max() > 0 else 0)
        ax2.scatter(lengths, breadths, s=sizes, color=ACCENT, alpha=0.75, edgecolors=CHROME["surface"], linewidths=1)
        if len(per_ref) <= 12:
            for ref, x, yv in zip(per_ref, lengths, breadths):
                ax2.annotate(design.truncate(ref, 16), (x, yv), xytext=(5, 4), textcoords="offset points",
                             fontsize=8, color=CHROME["ink_secondary"])
        if lengths.size and lengths.max() / max(lengths.min(), 1) > 20:
            ax2.set_xscale("log")
        ax2.xaxis.set_major_formatter(bp_formatter())
        ax2.set_ylim(-3, 105)
        ax2.set_xlabel("Reference length")
        ax2.set_ylabel("Coverage breadth (%)")
        ax2.set_title("Breadth against length (marker size: mean depth)")
        style_axes(ax2, grid="both")
        fig.tight_layout()
        return self._save(fig, "per_reference_coverage")

    def plot_coverage_heatmap(self) -> Optional[Path]:
        """Coverage map: mean depth in 100 relative-position bins for each reference."""
        if not self.coverage_arrays:
            fig, ax = plt.subplots(figsize=(11, 2.2))
            empty_panel(ax, "Coverage map needs per-position coverage (none was supplied)")
            ax.set_title("Coverage map")
            fig.tight_layout()
            return self._save(fig, "coverage_heatmap")
        refs = design.select_references({r: {'length': len(a)} for r, a in self.coverage_arrays.items()}, 40)
        n_bins = 100
        matrix = np.zeros((len(refs), n_bins))
        for i, ref in enumerate(refs):
            _x, y = design.bin_coverage(self.coverage_arrays[ref], n_bins)
            if y.size == n_bins:
                matrix[i] = y
            elif y.size:
                matrix[i] = np.interp(np.linspace(0, y.size - 1, n_bins), np.arange(y.size), y)
        fig_h = 0.28 * len(refs) + 1.8
        fig, ax = plt.subplots(figsize=(11, fig_h))
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("baitutils_blue", ["#ffffff"] + SEQUENTIAL)
        vmax = max(float(matrix.max()), self._min_coverage, 1.0)
        im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_yticks(range(len(refs)))
        ax.set_yticklabels([design.truncate(r, 26) for r in refs], fontsize=8)
        ax.set_xticks([0, 24.5, 49.5, 74.5, 99])
        ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
        ax.set_xlabel("Relative position along the reference")
        ax.tick_params(length=0)
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.grid(False)
        title = "Coverage map"
        if len(self.coverage_arrays) > len(refs):
            title += f" ({len(refs)} longest of {len(self.coverage_arrays)} references)"
        ax.set_title(title)
        cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("Mean depth per bin", fontsize=9, color=CHROME["ink_secondary"])
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(labelsize=8, colors=CHROME["ink_muted"])
        fig.tight_layout()
        return self._save(fig, "coverage_heatmap")

    def plot_quality_score(self) -> Optional[Path]:
        """Component scores with their weights, and the overall score against the category thresholds."""
        scores = self.quality_scores.get('component_scores', {}) or {}
        if not scores:
            logging.info("No quality scores; the quality figure is skipped")
            return None
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={'width_ratios': [2.2, 1]})
        weights = QualityScorer.DEFAULT_WEIGHTS
        names = list(scores)
        values = [scores[n] for n in names]
        y = np.arange(len(names))[::-1]
        ax1.barh(y, values, height=0.58, color=ACCENT)
        ax1.set_yticks(y)
        ax1.set_yticklabels([f"{n.replace('_', ' ').capitalize()}  ({weights.get(n, 0):.0%})" for n in names])
        for yi, value in zip(y, values):
            ax1.text(value + 0.015, yi, f"{value:.2f}", va="center", fontsize=8.5, color=CHROME["ink_secondary"])
        ax1.set_xlim(0, 1.12)
        ax1.set_xlabel("Component score (0-1); the weight of each component in brackets")
        ax1.set_title("Quality score components")
        style_axes(ax1, grid="x")
        ax1.spines["left"].set_visible(False)

        overall = float(self.quality_scores.get('overall_score', 0.0))
        category = str(self.quality_scores.get('category', ''))
        ax2.barh([0], [overall], height=0.5, color=design.category_color(category))
        thresholds = QualityScorer.DEFAULT_THRESHOLDS
        for name in ("fair", "good", "excellent"):
            t = thresholds[name]
            ax2.axvline(t, color=CHROME["ink_muted"], linewidth=0.8, linestyle=(0, (3, 3)))
            ax2.text(t, 0.42, name.capitalize(), fontsize=8, color=CHROME["ink_muted"], ha="center", va="bottom")
        ax2.text(min(overall + 0.02, 0.9), 0, f"{overall:.2f} {category}", va="center", fontsize=9.5,
                 color=CHROME["ink"], fontweight="semibold")
        ax2.set_xlim(0, 1)
        ax2.set_ylim(-0.6, 0.7)
        ax2.set_yticks([])
        ax2.set_xlabel("Overall score")
        ax2.set_title("Overall quality")
        style_axes(ax2, grid="x", baseline=False)
        fig.tight_layout()
        return self._save(fig, "quality_score")

    # Kept for callers of the previous API
    def _save_plot(self, fig, filename: str) -> None:
        self._save(fig, filename)
