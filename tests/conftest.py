"""
Shared test fixtures for baitUtils.

The fixtures build a small, deterministic dataset whose expected values can be
derived by hand:

Reference
    chrA: 3000 bp
    chrB: 2000 bp

Baits (all 120 bp unless noted)
    A_00..A_08   chrA, starts 0, 120, ..., 960      (9 baits, cover 0-1080)
    hole         chrA 1080-1380 has no bait          (300 bp gap)
    A_09..A_21   chrA, starts 1380, 1500, ..., 2820 (13 baits, cover 1380-2940)
    tail         chrA 2940-3000 uncovered           (60 bp)
    B_00..B_15   chrB, starts 0, 120, ..., 1800     (16 baits, cover 0-1920)
    tail         chrB 1920-2000 uncovered           (80 bp)
    G_00         chrB 500-630, two blocks of 60 bp with a 10 bp target insert
    L_00         chrA 500-620 with 20 mismatches    (identity 83.3 %)
    U_00..U_02   random, no PSL row                 (unmapped)

Totals
    43 baits, 39 with identity >= 90 % and length >= 100
    chrA covered 2640 / 3000, chrB covered 1920 / 2000
"""

import os
import random
import stat
import sys
from pathlib import Path

import pytest

BAIT_LEN = 120
CHR_LENGTHS = {"chrA": 3000, "chrB": 2000}
A_STARTS_1 = list(range(0, 1080, 120))            # 9 baits
A_STARTS_2 = list(range(1380, 2940, 120))         # 13 baits
B_STARTS = list(range(0, 1920, 120))              # 16 baits
HOLE = ("chrA", 1080, 1380)

EXPECTED = {
    "n_baits": 43,
    "n_mapped": 39,
    "n_unmapped": 4,
    "covered": {"chrA": 2640, "chrB": 1920},
    "hole": HOLE,
}


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng: random.Random, seq: str, n: int) -> str:
    seq = list(seq)
    for pos in rng.sample(range(len(seq)), n):
        seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
    return "".join(seq)


def _write_fasta(path: Path, records) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def _psl_row(matches, mismatches, strand, qname, qsize, qstart, qend,
             tname, tsize, tstart, tend, blocks):
    """Build one 21-column PSL row. blocks is a list of (size, qstart, tstart)."""
    q_num_insert = 0
    q_base_insert = 0
    t_num_insert = max(0, len(blocks) - 1)
    t_base_insert = 0
    for (s1, _, t1), (_, _, t2) in zip(blocks, blocks[1:]):
        t_base_insert += t2 - (t1 + s1)
    cols = [
        matches, mismatches, 0, 0,
        q_num_insert, q_base_insert, t_num_insert, t_base_insert,
        strand, qname, qsize, qstart, qend, tname, tsize, tstart, tend,
        len(blocks),
        ",".join(str(b[0]) for b in blocks) + ",",
        ",".join(str(b[1]) for b in blocks) + ",",
        ",".join(str(b[2]) for b in blocks) + ",",
    ]
    return "\t".join(str(c) for c in cols)


@pytest.fixture(scope="session")
def dataset(tmp_path_factory) -> dict:
    """Reference FASTA, bait FASTA and matching PSL with known coordinates."""
    rng = random.Random(20260925)
    root = tmp_path_factory.mktemp("dataset")

    ref = {name: _random_seq(rng, n) for name, n in CHR_LENGTHS.items()}
    baits = []
    psl_rows = []

    def add_tile(prefix, idx, chrom, start):
        name = f"{prefix}_{idx:02d}"
        seq = ref[chrom][start:start + BAIT_LEN]
        baits.append((name, seq))
        psl_rows.append(_psl_row(
            BAIT_LEN, 0, "+", name, BAIT_LEN, 0, BAIT_LEN,
            chrom, CHR_LENGTHS[chrom], start, start + BAIT_LEN,
            [(BAIT_LEN, 0, start)]))

    for i, s in enumerate(A_STARTS_1 + A_STARTS_2):
        add_tile("A", i, "chrA", s)
    for i, s in enumerate(B_STARTS):
        add_tile("B", i, "chrB", s)

    # Gapped alignment: two 60 bp blocks separated by a 10 bp target insert.
    g_seq = ref["chrB"][500:560] + ref["chrB"][570:630]
    baits.append(("G_00", g_seq))
    psl_rows.append(_psl_row(
        BAIT_LEN, 0, "+", "G_00", BAIT_LEN, 0, BAIT_LEN,
        "chrB", CHR_LENGTHS["chrB"], 500, 630,
        [(60, 0, 500), (60, 60, 570)]))

    # Low identity bait: 20 mismatches over 120 bp.
    l_seq = _mutate(rng, ref["chrA"][500:620], 20)
    baits.append(("L_00", l_seq))
    psl_rows.append(_psl_row(
        100, 20, "+", "L_00", BAIT_LEN, 0, BAIT_LEN,
        "chrA", CHR_LENGTHS["chrA"], 500, 620,
        [(BAIT_LEN, 0, 500)]))

    # Unmapped baits.
    for i in range(3):
        baits.append((f"U_{i:02d}", _random_seq(rng, BAIT_LEN)))

    ref_fa = root / "reference.fa"
    bait_fa = root / "baits.fa"
    psl = root / "hits.psl"
    _write_fasta(ref_fa, ref.items())
    _write_fasta(bait_fa, baits)
    psl.write_text("\n".join(psl_rows) + "\n")

    return {
        "root": root,
        "reference": ref_fa,
        "baits": bait_fa,
        "psl": psl,
        "reference_seqs": ref,
        "bait_seqs": dict(baits),
        "expected": EXPECTED,
    }


@pytest.fixture
def fake_pblat(dataset, tmp_path, monkeypatch) -> Path:
    """Put a fake pblat on PATH that writes the fixture PSL to its output argument."""
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "pblat"
    script.write_text(
        "#!/bin/sh\n"
        "if [ $# -eq 0 ]; then echo 'pblat - fake' >&2; exit 0; fi\n"
        "for last; do :; done\n"
        f"cp '{dataset['psl']}' \"$last\"\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return script


@pytest.fixture
def run_cli(monkeypatch):
    """Run the baitUtils CLI in-process with the given argument list."""
    from baitUtils.__main__ import main

    def _run(argv):
        monkeypatch.setattr(sys, "argv", ["baitUtils"] + [str(a) for a in argv])
        return main()

    return _run
