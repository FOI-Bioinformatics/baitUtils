"""
Value-level tests for SequenceAnalyzer thermodynamics and composition metrics.
"""

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from baitUtils.sequence_analysis import SequenceAnalyzer, HAS_VIENNA_RNA

SEQ20 = "ACGTTGCAACGTAGCTAGGC"
SEQ120 = ("ACGTTGCAACGTAGCTAGGCTTAGCCGATCGGATCCAGTAGCTAGCTTGACCGATAGGCAT"
          "TAGCCAGTAGGCTAGCTAGCTAGGCTAGCTAGCAGTCGATCGATCGATCGATTCGAGCAGT")


class TestMeltingTemperature:
    def test_matches_biopython_nearest_neighbour(self):
        analyzer = SequenceAnalyzer(na_equivalent=50.0, dnac1_equivalent=250.0, dnac2_equivalent=250.0)
        expected = round(mt.Tm_NN(SEQ20, Na=50.0, dnac1=250.0, dnac2=250.0), 2)
        assert analyzer.calculate_melting_temperature(SEQ20) == pytest.approx(expected)

    def test_returns_number_for_default_parameters(self):
        tm = SequenceAnalyzer().calculate_melting_temperature(SEQ120)
        assert isinstance(tm, float)
        assert 60 < tm < 95

    def test_salt_changes_tm(self):
        low = SequenceAnalyzer(na_equivalent=10.0).calculate_melting_temperature(SEQ120)
        high = SequenceAnalyzer(na_equivalent=200.0).calculate_melting_temperature(SEQ120)
        assert high > low

    def test_lowercase_input_is_accepted(self):
        upper = SequenceAnalyzer().calculate_melting_temperature(SEQ20)
        lower = SequenceAnalyzer().calculate_melting_temperature(SEQ20.lower())
        assert upper == pytest.approx(lower)


@pytest.mark.skipif(not HAS_VIENNA_RNA, reason="ViennaRNA not installed")
class TestSecondaryStructure:
    def test_hairpin_uses_dna_parameters(self):
        import RNA
        analyzer = SequenceAnalyzer(hybridization_temp=37.0)
        md = RNA.md()
        md.temperature = 37.0
        expected = RNA.fold_compound(SEQ120, md).mfe()[1]
        assert analyzer.calculate_hairpin_dg(SEQ120) == pytest.approx(expected, abs=1e-6)
        assert analyzer.calculate_mfe(SEQ120) == analyzer.calculate_hairpin_dg(SEQ120)

    def test_higher_temperature_gives_weaker_structure(self):
        cold = SequenceAnalyzer(hybridization_temp=37.0).calculate_hairpin_dg(SEQ120)
        hot = SequenceAnalyzer(hybridization_temp=65.0).calculate_hairpin_dg(SEQ120)
        assert hot >= cold

    def test_unstructured_sequence_is_near_zero(self):
        dg = SequenceAnalyzer(hybridization_temp=65.0).calculate_hairpin_dg("A" * 60 + "C" * 60)
        assert dg == pytest.approx(0.0, abs=0.5)

    def test_self_dimer_of_palindrome_is_strongly_negative(self):
        analyzer = SequenceAnalyzer(hybridization_temp=37.0)
        palindrome = "GAATTCGGATCCGAATTC" * 3   # self-complementary
        dimer = analyzer.calculate_self_dimer_dg(palindrome)
        hairpin = analyzer.calculate_hairpin_dg(palindrome)
        assert dimer < -20
        assert dimer <= hairpin

    def test_self_dimer_of_homopolymer_is_near_zero(self):
        dg = SequenceAnalyzer(hybridization_temp=65.0).calculate_self_dimer_dg("A" * 60)
        assert dg == pytest.approx(0.0, abs=0.5)

    def test_analyze_sequence_reports_both_energies(self):
        stats = SequenceAnalyzer().analyze_sequence(SEQ120, "x")
        assert isinstance(stats["hairpin_dg"], float)
        assert isinstance(stats["self_dimer_dg"], float)
        assert "mfe" not in stats


class TestCaseHandling:
    def test_homopolymer_runs_ignore_case(self):
        analyzer = SequenceAnalyzer()
        assert analyzer.count_homopolymer_runs("AAaaAA") == [("A", 0, 6)]

    def test_complexity_ignores_case(self):
        analyzer = SequenceAnalyzer()
        assert analyzer.calculate_complexity("ACGTacgt") == analyzer.calculate_complexity("ACGTACGT")

    def test_dinucleotide_bias_ignores_case(self):
        analyzer = SequenceAnalyzer()
        assert analyzer.calculate_dinucleotide_bias("ACGTacgtAC") == pytest.approx(
            analyzer.calculate_dinucleotide_bias("ACGTACGTAC"))

    def test_masked_count_keeps_case(self):
        assert SequenceAnalyzer().count_masked_bases("ACGTacgt") == 4
