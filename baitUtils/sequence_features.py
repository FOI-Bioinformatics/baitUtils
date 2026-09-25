"""
sequence_features.py

Vectorized sequence feature calculations shared by the reference analyzer
and the gap analyzer. All functions run in time linear in sequence length
and work on numpy arrays, so they scale to genome-sized references.

Sequences are encoded as uint8 arrays with A, C, G, T as 0..3, N (or any
other symbol) as 4. Lowercase input is treated as masked and reported
separately by the caller.
"""

from typing import Dict, Tuple

import numpy as np

_LOOKUP = np.full(256, 4, dtype=np.uint8)
for _base, _code in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    _LOOKUP[ord(_base)] = _code
    _LOOKUP[ord(_base.lower())] = _code

BASE_CODES = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}


def encode(sequence: str) -> np.ndarray:
    """Encode a sequence as uint8 codes (A=0, C=1, G=2, T=3, other=4)."""
    return _LOOKUP[np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)]


def masked_mask(sequence: str) -> np.ndarray:
    """Boolean array marking lowercase (soft-masked) positions."""
    raw = np.frombuffer(sequence.encode("ascii", "replace"), dtype=np.uint8)
    return (raw >= ord("a")) & (raw <= ord("z"))


def base_counts(codes: np.ndarray) -> np.ndarray:
    """Counts of A, C, G, T, other in order."""
    return np.bincount(codes, minlength=5)[:5]


def gc_content(codes: np.ndarray) -> float:
    """GC percentage of the A/C/G/T bases (N excluded from the denominator)."""
    counts = base_counts(codes)
    acgt = counts[:4].sum()
    return float((counts[1] + counts[2]) / acgt * 100.0) if acgt else 0.0


def n_content(codes: np.ndarray) -> float:
    """Percentage of positions that are not A, C, G or T."""
    return float(np.mean(codes == 4) * 100.0) if codes.size else 0.0


def shannon_entropy(codes: np.ndarray) -> float:
    """Shannon entropy (bits) of the A/C/G/T composition."""
    counts = base_counts(codes)[:4].astype(float)
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log2(p)).sum())


def run_lengths(codes: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Start, length and code of every run of identical codes."""
    if codes.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=np.uint8)
    change = np.flatnonzero(np.diff(codes)) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [codes.size]))
    return starts, ends - starts, codes[starts]


def max_homopolymer(codes: np.ndarray) -> int:
    """Longest run of one base (N runs excluded)."""
    starts, lengths, run_codes = run_lengths(codes)
    lengths = lengths[run_codes != 4]
    return int(lengths.max()) if lengths.size else 0


def homopolymer_mask(codes: np.ndarray, min_run: int = 5) -> np.ndarray:
    """Boolean array marking positions inside homopolymer runs of at least min_run bases."""
    mask = np.zeros(codes.size, dtype=bool)
    starts, lengths, run_codes = run_lengths(codes)
    keep = (lengths >= min_run) & (run_codes != 4)
    for start, length in zip(starts[keep], lengths[keep]):
        mask[start:start + length] = True
    return mask


def kmer_codes(codes: np.ndarray, k: int) -> np.ndarray:
    """Integer code of every k-mer without N; positions with N give -1."""
    n = codes.size - k + 1
    if n <= 0:
        return np.array([], dtype=np.int64)
    values = np.zeros(n, dtype=np.int64)
    invalid = np.zeros(n, dtype=bool)
    for offset in range(k):
        window = codes[offset:offset + n]
        values = values * 4 + np.minimum(window, 3)
        invalid |= window == 4
    values[invalid] = -1
    return values


def duplicated_kmer_mask(codes: np.ndarray, k: int = 12) -> np.ndarray:
    """
    Boolean array over k-mer start positions marking k-mers that occur more
    than once in the sequence. Used as a fast proxy for repetitive content.
    """
    kmers = kmer_codes(codes, k)
    if kmers.size == 0:
        return np.array([], dtype=bool)
    valid = kmers >= 0
    unique, inverse, counts = np.unique(kmers[valid], return_inverse=True, return_counts=True)
    mask = np.zeros(kmers.size, dtype=bool)
    mask[valid] = counts[inverse] > 1
    return mask


def repeat_content(codes: np.ndarray, k: int = 12) -> float:
    """Percentage of k-mer positions whose k-mer occurs more than once."""
    mask = duplicated_kmer_mask(codes, k)
    return float(mask.mean() * 100.0) if mask.size else 0.0


def linguistic_complexity(codes: np.ndarray, max_k: int = 4) -> float:
    """
    Product over k = 1..max_k of observed distinct k-mers divided by the
    maximum possible (min(4^k, L - k + 1)); 1.0 for maximally diverse input.
    """
    if codes.size == 0:
        return 0.0
    value = 1.0
    for k in range(1, max_k + 1):
        kmers = kmer_codes(codes, k)
        kmers = kmers[kmers >= 0]
        if kmers.size == 0:
            return 0.0
        possible = min(4 ** k, kmers.size)
        value *= np.unique(kmers).size / possible
    return float(value)


def _window_sums(indicator: np.ndarray, window: int, step: int) -> np.ndarray:
    """Sum of an indicator array over windows [start, start + window) for each start."""
    cumulative = np.concatenate(([0], np.cumsum(indicator, dtype=np.int64)))
    starts = window_starts(indicator.size, window, step)
    ends = np.minimum(starts + window, indicator.size)
    return cumulative[ends] - cumulative[starts]


def window_starts(length: int, window: int, step: int) -> np.ndarray:
    """Start positions of windows; the last window may be shorter than window."""
    if length <= 0:
        return np.array([], dtype=np.int64)
    starts = np.arange(0, max(length - window, 0) + 1, step, dtype=np.int64)
    if starts.size == 0 or starts[-1] + window < length:
        tail = starts[-1] + step if starts.size else 0
        if tail < length:
            starts = np.append(starts, tail)
    return starts


def window_features(codes: np.ndarray, window: int, step: int,
                    homopolymer_min_run: int = 5, repeat_k: int = 12) -> Dict[str, np.ndarray]:
    """
    Per-window features as arrays aligned with window_starts(len, window, step):
    start, end, gc_content (%), n_content (%), entropy (bits),
    homopolymer_content (% bases in runs >= homopolymer_min_run),
    repeat_density (% k-mer positions that are duplicated in the sequence).
    """
    starts = window_starts(codes.size, window, step)
    ends = np.minimum(starts + window, codes.size)
    lengths = (ends - starts).astype(float)
    if starts.size == 0:
        return {key: np.array([]) for key in
                ("start", "end", "gc_content", "n_content", "entropy", "homopolymer_content", "repeat_density")}

    counts = np.stack([_window_sums(codes == code, window, step) for code in range(5)], axis=1).astype(float)
    acgt = counts[:, :4].sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        gc = np.where(acgt > 0, (counts[:, 1] + counts[:, 2]) / acgt * 100.0, 0.0)
        p = np.where(acgt[:, None] > 0, counts[:, :4] / np.maximum(acgt[:, None], 1), 0.0)
        entropy = -np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0.0).sum(axis=1)
    n_frac = counts[:, 4] / lengths * 100.0
    homopolymer = _window_sums(homopolymer_mask(codes, homopolymer_min_run), window, step) / lengths * 100.0

    dup = duplicated_kmer_mask(codes, repeat_k)
    if dup.size:
        padded = np.concatenate((dup, np.zeros(codes.size - dup.size, dtype=bool)))
        kmer_positions = np.maximum(_window_sums(np.ones(codes.size, dtype=np.int8), window, step) - repeat_k + 1, 1)
        repeat = _window_sums(padded, window, step) / kmer_positions * 100.0
    else:
        repeat = np.zeros(starts.size)

    return {
        "start": starts, "end": ends, "gc_content": gc, "n_content": n_frac, "entropy": entropy,
        "homopolymer_content": homopolymer, "repeat_density": np.minimum(repeat, 100.0),
    }


def window_coverage(coverage: np.ndarray, window: int, step: int, min_coverage: float = 1.0) -> Dict[str, np.ndarray]:
    """Mean depth and breadth (%) per window of a per-base coverage array."""
    coverage = np.asarray(coverage)
    starts = window_starts(coverage.size, window, step)
    ends = np.minimum(starts + window, coverage.size)
    lengths = (ends - starts).astype(float)
    if starts.size == 0:
        return {"mean_depth": np.array([]), "breadth": np.array([])}
    depth_sum = _window_sums(coverage.astype(np.int64), window, step)
    covered = _window_sums(coverage >= min_coverage, window, step)
    return {"mean_depth": depth_sum / lengths, "breadth": covered / lengths * 100.0}


def sequence_summary(sequence: str, repeat_k: int = 12) -> Dict[str, float]:
    """Sequence-level features for one reference sequence."""
    codes = encode(sequence)
    masked = masked_mask(sequence)
    counts = base_counts(codes)
    acgt = int(counts[:4].sum())
    length = int(codes.size)
    return {
        "length": length,
        "gc_content": gc_content(codes),
        "at_content": float((counts[0] + counts[3]) / acgt * 100.0) if acgt else 0.0,
        "n_content": n_content(codes),
        "masked_content": float(masked.mean() * 100.0) if length else 0.0,
        "shannon_entropy": shannon_entropy(codes),
        "linguistic_complexity": linguistic_complexity(codes),
        "max_homopolymer": max_homopolymer(codes),
        "homopolymer_content": float(homopolymer_mask(codes).mean() * 100.0) if length else 0.0,
        "repeat_content": repeat_content(codes, repeat_k),
    }
