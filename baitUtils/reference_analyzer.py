#!/usr/bin/env python3

"""
reference_analyzer.py

Reference sequence features and their relation to observed coverage.

Sequence-level features are computed per reference with vectorized numpy
routines (sequence_features.py). Window-level features are computed on
non-overlapping windows and, when per-base coverage arrays are supplied,
each window also carries its mean depth and breadth. Correlations between
features and coverage are Spearman rank correlations across all windows,
with the number of windows and a p-value reported. Without coverage arrays
the correlations are reported as not applicable.

Challenging references are flagged from sequence-level features with stated
thresholds; challenging windows are windows with low breadth that also show
an extreme feature value.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats
from Bio import SeqIO

from baitUtils.sequence_features import (
    encode, sequence_summary, window_features, window_coverage,
)


class ReferenceAnalyzer:
    """Analyze reference sequence features and their relation to coverage."""

    # Thresholds used to flag challenging sequences and windows
    EXTREME_GC_LOW = 20.0
    EXTREME_GC_HIGH = 80.0
    LOW_ENTROPY = 1.0
    HIGH_REPEAT = 30.0
    LONG_HOMOPOLYMER = 10
    HIGH_N = 5.0
    LOW_BREADTH_WINDOW = 50.0
    WINDOW_FEATURES = ("gc_content", "entropy", "homopolymer_content", "repeat_density", "n_content")

    def __init__(
        self,
        reference_file: Path,
        coverage_data: Dict[str, Any],
        window_size: int = 1000,
        overlap: int = 0,
        coverage_arrays: Optional[Dict[str, np.ndarray]] = None,
        min_coverage: float = 1.0,
        max_challenging_windows: int = 200,
    ):
        """
        Args:
            reference_file: Reference FASTA
            coverage_data: Output of CoverageAnalyzer.analyze (per-reference stats)
            window_size: Window length in bp
            overlap: Overlap between windows in bp (0 gives non-overlapping windows)
            coverage_arrays: Per-base depth arrays keyed by reference id
            min_coverage: Depth at which a base counts as covered
            max_challenging_windows: Cap on reported challenging windows per reference
        """
        self.reference_file = Path(reference_file)
        self.coverage_data = coverage_data or {}
        self.window_size = int(window_size)
        self.overlap = int(overlap)
        self.step_size = max(1, self.window_size - self.overlap)
        self.coverage_arrays = coverage_arrays or {}
        self.min_coverage = min_coverage
        self.max_challenging_windows = max_challenging_windows

        self.reference_sequences: Dict[str, str] = {}
        self.sequence_features: Dict[str, Dict[str, float]] = {}
        self.window_features: Dict[str, Dict[str, np.ndarray]] = {}
        self.coverage_correlations: Dict[str, Any] = {}
        self.challenging_regions: Dict[str, Dict[str, Any]] = {}
        self.challenging_windows: Dict[str, List[Dict[str, Any]]] = {}

    # ------------------------------------------------------------------ main
    def analyze(self) -> Dict[str, Any]:
        """Run the analysis and return a dictionary of results."""
        logging.info("Analyzing reference sequence features...")
        self._load_reference_sequences()
        self._analyze_sequence_features()
        self._perform_window_analysis()
        self._correlate_with_coverage()
        self._identify_challenging_regions()
        return {
            "sequence_features": self.sequence_features,
            "window_features": {ref: self._windows_as_records(ref) for ref in self.window_features},
            "coverage_correlations": self.coverage_correlations,
            "challenging_regions": self.challenging_regions,
            "challenging_windows": self.challenging_windows,
            "analysis_summary": self._generate_summary(),
        }

    def _load_reference_sequences(self) -> None:
        with open(self.reference_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.reference_sequences[record.id] = str(record.seq)
        logging.info(f"Loaded {len(self.reference_sequences)} reference sequences")

    # -------------------------------------------------------------- features
    def _analyze_sequence_features(self) -> None:
        for ref_id, sequence in self.reference_sequences.items():
            features = sequence_summary(sequence)
            self.sequence_features[ref_id] = features
        # problematic_score needs the window features; filled in below
        logging.info("Sequence-level features analyzed")

    def _perform_window_analysis(self) -> None:
        for ref_id, sequence in self.reference_sequences.items():
            codes = encode(sequence)
            windows = window_features(codes, self.window_size, self.step_size)
            if ref_id in self.coverage_arrays:
                windows.update(window_coverage(self.coverage_arrays[ref_id], self.window_size,
                                               self.step_size, self.min_coverage))
            self.window_features[ref_id] = windows

            flagged = self._flag_windows(windows)
            n_windows = int(windows["start"].size)
            self.sequence_features[ref_id]["problematic_score"] = (
                float(flagged.any(axis=0).mean() * 100.0) if n_windows else 0.0)
            self.sequence_features[ref_id]["windows"] = n_windows
        logging.info("Window-based analysis completed")

    def _flag_windows(self, windows: Dict[str, np.ndarray]) -> np.ndarray:
        """Boolean matrix (5 flags x windows) of extreme feature values."""
        n = windows["start"].size
        if n == 0:
            return np.zeros((5, 0), dtype=bool)
        return np.stack([
            (windows["gc_content"] < self.EXTREME_GC_LOW) | (windows["gc_content"] > self.EXTREME_GC_HIGH),
            windows["entropy"] < self.LOW_ENTROPY,
            windows["repeat_density"] > self.HIGH_REPEAT,
            windows["homopolymer_content"] > 10.0,
            windows["n_content"] > self.HIGH_N,
        ])

    FLAG_NAMES = ("extreme GC", "low complexity", "repetitive", "homopolymer-rich", "N-rich")

    # ---------------------------------------------------------- correlations
    def _correlate_with_coverage(self) -> None:
        """Spearman correlation between window features and window coverage, pooled over references."""
        refs = [r for r in self.window_features if "breadth" in self.window_features[r]]
        if not refs:
            self.coverage_correlations = {
                "applicable": False,
                "reason": "no per-base coverage arrays supplied",
                "windows": 0, "window_size": self.window_size, "features": {},
            }
            logging.info("Coverage correlations not applicable: no coverage arrays")
            return

        pooled = {key: np.concatenate([self.window_features[r][key] for r in refs])
                  for key in self.WINDOW_FEATURES + ("breadth", "mean_depth")}
        n = int(pooled["breadth"].size)
        features: Dict[str, Dict[str, float]] = {}
        for name in self.WINDOW_FEATURES:
            x = pooled[name]
            entry: Dict[str, float] = {}
            for target in ("breadth", "mean_depth"):
                y = pooled[target]
                if n < 3 or np.all(x == x[0]) or np.all(y == y[0]):
                    entry[f"rho_{target}"] = float("nan")
                    entry[f"p_{target}"] = float("nan")
                else:
                    rho, p = stats.spearmanr(x, y)
                    entry[f"rho_{target}"] = float(rho)
                    entry[f"p_{target}"] = float(p)
            features[name] = entry

        # Breadth by GC decile, a direct view of composition effects
        gc = pooled["gc_content"]
        edges = np.arange(0, 101, 10)
        gc_bins = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (gc >= lo) & ((gc < hi) if hi < 100 else (gc <= hi))
            if sel.any():
                gc_bins.append({
                    "gc_min": int(lo), "gc_max": int(hi), "windows": int(sel.sum()),
                    "mean_breadth": float(pooled["breadth"][sel].mean()),
                    "mean_depth": float(pooled["mean_depth"][sel].mean()),
                })

        self.coverage_correlations = {
            "applicable": n >= 3,
            "method": "Spearman rank correlation across windows",
            "windows": n, "window_size": self.window_size,
            "features": features, "breadth_by_gc": gc_bins,
        }
        logging.info(f"Coverage-feature correlations computed over {n} windows")

    # ------------------------------------------------------------ challenging
    def _identify_challenging_regions(self) -> None:
        for ref_id, features in self.sequence_features.items():
            score = 0
            reasons = []
            gc = features["gc_content"]
            if gc < self.EXTREME_GC_LOW or gc > self.EXTREME_GC_HIGH:
                score += 2
                reasons.append(f"Extreme GC content ({gc:.1f}%)")
            if features["shannon_entropy"] < self.LOW_ENTROPY:
                score += 2
                reasons.append("Low sequence complexity")
            if features["repeat_content"] > self.HIGH_REPEAT:
                score += 1
                reasons.append(f"High repetitive content ({features['repeat_content']:.1f}% duplicated 12-mers)")
            if features["max_homopolymer"] > self.LONG_HOMOPOLYMER:
                score += 1
                reasons.append(f"Long homopolymer runs (max {features['max_homopolymer']} bp)")
            if features["n_content"] > self.HIGH_N:
                score += 1
                reasons.append(f"High N content ({features['n_content']:.1f}%)")
            if features.get("problematic_score", 0.0) > 20.0:
                score += 1
                reasons.append(f"{features['problematic_score']:.0f}% of windows show an extreme feature")
            self.challenging_regions[ref_id] = {
                "challenging_score": score,
                "reasons": reasons,
                "difficulty_level": self._get_difficulty_level(score),
            }

            windows = self.window_features.get(ref_id, {})
            if "breadth" in windows and windows["start"].size:
                flagged = self._flag_windows(windows)
                low = windows["breadth"] < self.LOW_BREADTH_WINDOW
                hits = np.flatnonzero(low & flagged.any(axis=0))
                self.challenging_windows[ref_id] = [{
                    "start": int(windows["start"][i]), "end": int(windows["end"][i]),
                    "breadth": float(windows["breadth"][i]),
                    "reasons": [name for name, f in zip(self.FLAG_NAMES, flagged[:, i]) if f],
                } for i in hits[:self.max_challenging_windows]]

    @staticmethod
    def _get_difficulty_level(score: int) -> str:
        if score >= 5:
            return "Very Difficult"
        if score >= 3:
            return "Difficult"
        if score >= 1:
            return "Moderate"
        return "Easy"

    # ---------------------------------------------------------------- output
    def _windows_as_records(self, ref_id: str) -> List[Dict[str, float]]:
        windows = self.window_features[ref_id]
        keys = list(windows)
        return [{k: (int(windows[k][i]) if k in ("start", "end") else float(windows[k][i])) for k in keys}
                for i in range(int(windows["start"].size))]

    def _generate_summary(self) -> Dict[str, Any]:
        if not self.sequence_features:
            return {}
        feats = list(self.sequence_features.values())
        return {
            "total_sequences": len(feats),
            "total_length": int(sum(f["length"] for f in feats)),
            "mean_gc_content": float(np.mean([f["gc_content"] for f in feats])),
            "mean_complexity": float(np.mean([f["shannon_entropy"] for f in feats])),
            "sequences_with_extreme_gc": sum(1 for f in feats
                                             if f["gc_content"] < self.EXTREME_GC_LOW or f["gc_content"] > self.EXTREME_GC_HIGH),
            "sequences_with_high_repeats": sum(1 for f in feats if f["repeat_content"] > self.HIGH_REPEAT),
            "challenging_sequences": sum(1 for r in self.challenging_regions.values() if r["challenging_score"] >= 3),
            "challenging_windows": sum(len(w) for w in self.challenging_windows.values()),
            "window_size": self.window_size,
            "correlations_applicable": bool(self.coverage_correlations.get("applicable", False)),
        }
