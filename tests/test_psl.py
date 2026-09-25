"""
Tests for the shared PSL parser in mapping_utils.
"""

import pytest

from baitUtils.mapping_utils import PSLHit, parse_psl, parse_psl_line, filter_hits, PSLParser


def row(*cols):
    return "\t".join(str(c) for c in cols)


PERFECT = row(120, 0, 0, 0, 0, 0, 0, 0, "+", "b0", 120, 0, 120, "chrA", 3000, 0, 120, 1, "120,", "0,", "0,")
MISMATCH = row(100, 20, 0, 0, 0, 0, 0, 0, "+", "L", 120, 0, 120, "chrA", 3000, 500, 620, 1, "120,", "0,", "500,")
T_INSERT = row(120, 0, 0, 0, 0, 0, 1, 10, "+", "G", 120, 0, 120, "chrB", 2000, 500, 630, 2, "60,60,", "0,60,", "500,570,")
Q_INSERT = row(110, 0, 0, 0, 1, 10, 0, 0, "-", "Q", 120, 0, 120, "chrB", 2000, 100, 210, 2, "50,60,", "0,60,", "100,150,")
PARTIAL = row(40, 0, 0, 0, 0, 0, 0, 0, "+", "P", 120, 80, 120, "chrA", 3000, 900, 940, 1, "40,", "80,", "900,")


class TestParseLine:
    def test_all_fields_parsed(self):
        hit = parse_psl_line(T_INSERT)
        assert isinstance(hit, PSLHit)
        assert (hit.matches, hit.t_num_insert, hit.t_base_insert) == (120, 1, 10)
        assert hit.strand == "+"
        assert hit.block_sizes == [60, 60]
        assert hit.t_starts == [500, 570]
        assert hit.target_blocks == [(500, 560), (570, 630)]
        assert hit.aligned_length == 120
        assert hit.target_span == 130

    def test_short_line_is_rejected(self):
        assert parse_psl_line("120\t0\t0") is None

    def test_non_numeric_is_rejected(self):
        assert parse_psl_line(PERFECT.replace("120\t0\t0\t0", "x\t0\t0\t0", 1)) is None

    def test_space_separated_line_is_accepted(self):
        assert parse_psl_line(PERFECT.replace("\t", " ")).q_name == "b0"


class TestIdentity:
    def test_perfect_hit(self):
        assert parse_psl_line(PERFECT).identity == pytest.approx(100.0)

    def test_mismatches_reduce_identity(self):
        assert parse_psl_line(MISMATCH).identity == pytest.approx(100 - 1000 * 20 / 120 / 10)

    def test_target_insert_not_penalised_for_mrna_mode(self):
        hit = parse_psl_line(T_INSERT)
        assert hit.identity == pytest.approx(100.0)
        assert hit.milli_bad(is_mrna=False) > 0

    def test_query_insert_is_penalised(self):
        hit = parse_psl_line(Q_INSERT)
        # milliBad = 1000 * (0 mismatches + 1 insert + round(3 * ln(1 + 10))) / 110
        expected = 1000.0 * (1 + round(3 * __import__("math").log(11))) / 110
        assert hit.milli_bad() == pytest.approx(expected)
        assert hit.identity == pytest.approx(100 - expected / 10)


class TestParseFile:
    def test_headers_blank_and_malformed_lines_skipped(self, tmp_path, caplog):
        psl = tmp_path / "hits.psl"
        psl.write_text("psLayout version 3\n\nmatch\tmis-\n-----\n" + PERFECT + "\nbad line\n" + MISMATCH + "\n")
        hits = list(parse_psl(psl))
        assert [h.q_name for h in hits] == ["b0", "L"]
        assert "Skipped 1 malformed" in caplog.text

    def test_filter_hits(self, tmp_path):
        psl = tmp_path / "hits.psl"
        psl.write_text("\n".join([PERFECT, MISMATCH, T_INSERT, PARTIAL]) + "\n")
        names = [h.q_name for h in filter_hits(parse_psl(psl), min_identity=90, min_length=100)]
        assert names == ["b0", "G"]
        names = [h.q_name for h in filter_hits(parse_psl(psl), min_identity=0, min_length=0, min_matches=50)]
        assert names == ["b0", "L", "G"]


class TestPSLParserMappedIds:
    def test_mapped_ids_and_filtered_output(self, tmp_path):
        psl = tmp_path / "hits.psl"
        psl.write_text("psLayout version 3\n\n" + "\n".join([PERFECT, MISMATCH, PARTIAL]) + "\n")
        out = tmp_path / "filtered.psl"
        mapped = PSLParser.parse_psl_file(str(psl), min_identity=90, min_match_count=0,
                                          filtered_output=str(out), min_length=100)
        assert mapped == {"b0"}
        text = out.read_text()
        assert text.startswith("psLayout version 3")
        assert "\tb0\t" in text and "\tL\t" not in text


class TestHitTable:
    def test_hit_table_ranks_hits_per_bait(self, tmp_path):
        from baitUtils.mapping_utils import build_hit_table
        second = row(110, 10, 0, 0, 0, 0, 0, 0, "-", "b0", 120, 0, 120, "chrB", 2000, 40, 160, 1, "120,", "0,", "40,")
        psl = tmp_path / "hits.psl"
        psl.write_text("\n".join([PERFECT, second, MISMATCH, T_INSERT]) + "\n")
        table = build_hit_table(parse_psl(psl))
        assert list(table["oligo_id"]) == ["G", "L", "b0"]
        b0 = table.set_index("oligo_id").loc["b0"]
        assert b0["n_hits"] == 2 and b0["n_targets"] == 2
        assert b0["best_identity"] == 100.0
        assert b0["second_best_identity"] == pytest.approx(100 - 1000 * 10 / 120 / 10, abs=0.01)
        assert (b0["best_target"], b0["best_start"], b0["best_end"], b0["best_strand"]) == ("chrA", 0, 120, "+")
        g = table.set_index("oligo_id").loc["G"]
        assert g["n_hits"] == 1 and g["second_best_identity"] != g["second_best_identity"]  # NaN

    def test_empty_hit_table_has_columns(self):
        from baitUtils.mapping_utils import build_hit_table
        table = build_hit_table([])
        assert len(table) == 0 and "best_identity" in table.columns
