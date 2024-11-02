# tests/test_stats.py

import unittest
from unittest.mock import patch, mock_open
from io import StringIO
from baitUtils import stats
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import os

class TestStats(unittest.TestCase):
    def setUp(self):
        # Sample sequences for testing
        self.seq1 = SeqRecord(Seq("ATGCGC"), id="bait1")
        self.seq2 = SeqRecord(Seq("ATGCGCGC"), id="bait2")
        self.seq3 = SeqRecord(Seq("AAATTTCCCGGG"), id="bait3")
        self.sequences = [self.seq1, self.seq2, self.seq3]
    
    def test_calculate_gc_percentage(self):
        seq = "ATGCGC"
        expected_gc = 50.0
        gc = stats.calculate_gc_percentage(seq)
        self.assertEqual(gc, expected_gc)
    
    def test_calculate_tm(self):
        seq = "ATGCGC"
        # AT=2, GC=2: TM = 2*2 + 4*2 = 4 + 8 = 12
        expected_tm = 12
        tm = stats.calculate_tm(seq)
        self.assertEqual(tm, expected_tm)
    
    def test_calculate_masked_percentage(self):
        # Assuming masked as lowercase
        seq = "ATGcGc"
        # Number of lowercase letters: 2
        expected_masked = (2 / 6) * 100  # 33.333...
        masked = stats.calculate_masked_percentage(seq)
        self.assertAlmostEqual(masked, 33.3333, places=4)
    
    def test_max_homopolymer(self):
        seq = "AAACCCGGGTTT"
        expected_max = 3
        max_run = stats.max_homopolymer(seq)
        self.assertEqual(max_run, expected_max)
    
    def test_calculate_sequence_complexity(self):
        # Shannon entropy for "AABB"
        seq = "AABB"
        # P(A)=0.5, P(B)=0.5
        expected_entropy = - (0.5 * (-1) + 0.5 * (-1))  # 1.0
        entropy = stats.calculate_sequence_complexity(seq)
        self.assertEqual(entropy, 1.0)
    
        # Check with higher entropy
        seq = "ABCD"
        # P(A)=0.25, P(B)=0.25, P(C)=0.25, P(D)=0.25
        expected_entropy = 2.0  # -4*(0.25*log2(0.25)) = 2
        entropy = stats.calculate_sequence_complexity(seq)
        self.assertEqual(entropy, 2.0)
    
    def test_calculate_mfe(self):
        # Since it's a placeholder, we expect the dummy value
        seq = "ATGCGC"
        expected_mfe = -10.0
        mfe = stats.calculate_mfe(seq)
        self.assertEqual(mfe, expected_mfe)
    
    def test_calculate_self_score(self):
        # Placeholder
        seq = "ATGCGC"
        expected_score = 5.0
        score = stats.calculate_self_score(seq)
        self.assertEqual(score, expected_score)
    
    def test_calculate_entropy(self):
        # Placeholder
        seq = "AABBCC"
        expected_entropy = 1.5
        entropy = stats.calculate_entropy(seq)
        self.assertEqual(entropy, expected_entropy)
    
    def test_count_n(self):
        seq = "AANNCCTT"
        expected_n = 2
        count = stats.count_n(seq)
        self.assertEqual(count, expected_n)
    
    def test_count_gaps(self):
        seq = "AATG--CC"
        expected_gaps = 2
        count = stats.count_gaps(seq)
        self.assertEqual(count, expected_gaps)
    
    def test_should_keep(self):
        # Should keep: Length >= 100 and 40 <= GC% <= 60
        seq = "ATGCGC" * 20  # Length=120, GC% = 50%
        self.assertTrue(stats.should_keep(seq))
    
        # Should not keep: Length < 100
        seq_short = "ATGC" * 10  # Length=40, GC% = 50%
        self.assertFalse(stats.should_keep(seq_short))
    
        # Should not keep: GC% < 40
        seq_low_gc = "ATATAT" * 20  # Length=120, GC% = 0%
        self.assertFalse(stats.should_keep(seq_low_gc))
    
        # Should not keep: GC% > 60
        seq_high_gc = "GCGCGC" * 20  # Length=120, GC% = 100%
        self.assertFalse(stats.should_keep(seq_high_gc))
    
    def test_read_fasta(self):
        # Mock reading an uncompressed FASTA file
        with patch('builtins.open', mock_open(read_data=">bait1\nATGCGC\n>bait2\nATGCGCGC\n")):
            sequences = stats.read_fasta("dummy.fasta")
            self.assertEqual(len(sequences), 2)
            self.assertEqual(sequences[0].id, "bait1")
            self.assertEqual(str(sequences[0].seq), "ATGCGC")
            self.assertEqual(sequences[1].id, "bait2")
            self.assertEqual(str(sequences[1].seq), "ATGCGCGC")
    
        # Mock reading a gzipped FASTA file
        with patch('gzip.open', mock_open(read_data=">bait3\nAAATTTCCCGGG\n>bait4\nGCGCGCGC\n")):
            sequences = stats.read_fasta("dummy.fasta.gz")
            self.assertEqual(len(sequences), 2)
            self.assertEqual(sequences[0].id, "bait3")
            self.assertEqual(str(sequences[0].seq), "AAATTTCCCGGG")
            self.assertEqual(sequences[1].id, "bait4")
            self.assertEqual(str(sequences[1].seq), "GCGCGCGC")
    
    def test_calculate_parameters(self):
        # Using the sample sequences
        params_df = stats.calculate_parameters(self.sequences)
        self.assertEqual(len(params_df), 3)
        self.assertListEqual(list(params_df['Bait']), ['bait1', 'bait2', 'bait3'])
        self.assertListEqual(list(params_df['Length']), [6, 8, 12])
        self.assertAlmostEqual(params_df.loc[0, 'GC%'], 50.0)
        self.assertAlmostEqual(params_df.loc[1, 'GC%'], 75.0)
        self.assertAlmostEqual(params_df.loc[2, 'GC%'], 50.0)
        self.assertEqual(params_df.loc[0, 'Kept'], 'Yes')
        self.assertEqual(params_df.loc[1, 'Kept'], 'No')
        self.assertEqual(params_df.loc[2, 'Kept'], 'Yes')
    
    def test_filter_baits(self):
        # Create a DataFrame with varying parameters
        data = {
            'Bait': ['bait1', 'bait2', 'bait3', 'bait4'],
            'Length': [120, 80, 150, 90],
            'GC%': [50.0, 70.0, 45.0, 55.0],
            'Kept': ['Yes', 'No', 'Yes', 'No']
        }
        df = pd.DataFrame(data)
        filtered_df = stats.filter_baits(df)
        # Should keep only bait1 and bait3
        self.assertEqual(len(filtered_df), 2)
        self.assertListEqual(list(filtered_df['Bait']), ['bait1', 'bait3'])
    
    def test_main_without_filter(self):
        # Mock arguments for the main function without filtering
        args = argparse.Namespace(
            input='dummy.fasta',
            output='params.txt',
            filter=False,
            log=False
        )
        # Mock reading the FASTA file
        with patch('builtins.open', mock_open(read_data=">bait1\nATGCGC\n>bait2\nATGCGCGC\n")):
            with patch('pandas.DataFrame.to_csv') as mock_to_csv:
                stats.main(args)
                mock_to_csv.assert_called_once_with('params.txt', sep='\t', index=False)
    
    def test_main_with_filter(self):
        # Mock arguments for the main function with filtering
        args = argparse.Namespace(
            input='dummy.fasta',
            output='params_filtered.txt',
            filter=True,
            log=False
        )
        # Mock reading the FASTA file
        with patch('builtins.open', mock_open(read_data=">bait1\nATGCGC" * 20 + "\n>bait2\nATGCGCGC\n>bait3\nAAATTTCCCGGG\n")):
            with patch('pandas.DataFrame.to_csv') as mock_to_csv:
                stats.main(args)
                mock_to_csv.assert_called_once_with('params_filtered.txt', sep='\t', index=False)
    
    def test_main_with_logging(self):
        # Mock arguments with logging enabled
        args = argparse.Namespace(
            input='dummy.fasta',
            output='params.txt',
            filter=False,
            log=True
        )
        # Mock reading the FASTA file
        with patch('builtins.open', mock_open(read_data=">bait1\nATGCGC\n>bait2\nATGCGCGC\n")):
            with patch('pandas.DataFrame.to_csv') as mock_to_csv:
                with patch('logging.basicConfig') as mock_logging:
                    stats.main(args)
                    mock_logging.assert_called_once_with(
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S'
                    )
                mock_to_csv.assert_called_once_with('params.txt', sep='\t', index=False)
    
    def test_calculate_parameters_with_empty_sequences(self):
        # Test calculate_parameters with empty sequence list
        sequences = []
        params_df = stats.calculate_parameters(sequences)
        self.assertTrue(params_df.empty)
    
    def test_filter_baits_with_all_baits_kept(self):
        # All baits meet the criteria
        data = {
            'Bait': ['bait1', 'bait2'],
            'Length': [120, 150],
            'GC%': [50.0, 45.0],
            'Kept': ['Yes', 'Yes']
        }
        df = pd.DataFrame(data)
        filtered_df = stats.filter_baits(df)
        self.assertEqual(len(filtered_df), 2)
    
    def test_filter_baits_with_no_baits_kept(self):
        # No baits meet the criteria
        data = {
            'Bait': ['bait1', 'bait2'],
            'Length': [80, 90],
            'GC%': [30.0, 70.0],
            'Kept': ['No', 'No']
        }
        df = pd.DataFrame(data)
        filtered_df = stats.filter_baits(df)
        self.assertEqual(len(filtered_df), 0)
    
    def test_read_fasta_invalid_format(self):
        # Mock reading an invalid FASTA file
        with patch('builtins.open', mock_open(read_data="Invalid FASTA content")):
            with self.assertRaises(SystemExit) as cm:
                stats.read_fasta("invalid.fasta")
            self.assertEqual(cm.exception.code, 1)
    
    def test_should_keep_edge_cases(self):
        # Exactly at threshold
        seq = "ATGCGC" * 17  # Length=102, GC% = 50.0
        self.assertTrue(stats.should_keep(seq))
    
        # Just below length threshold
        seq = "ATGCGC" * 16 + "A"  # Length=97, GC% = 50.0
        self.assertFalse(stats.should_keep(seq))
    
        # Just below GC% threshold
        seq = "ATATAT" * 17  # Length=102, GC% = 0.0
        self.assertFalse(stats.should_keep(seq))
    
        # Just above GC% threshold
        seq = "GCGCGC" * 17  # Length=102, GC% = 100.0
        self.assertFalse(stats.should_keep(seq))
    
    # Additional tests can be added for more comprehensive coverage
    
if __name__ == '__main__':
    unittest.main()