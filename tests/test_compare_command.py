#!/usr/bin/env python3

"""
test_compare_command.py

Unit tests for the compare command (main Phase 3 entry point).
"""

import unittest
import tempfile
import shutil
import sys
import os
from pathlib import Path
from unittest.mock import patch, MagicMock, call
import argparse

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from baitUtils.compare import add_arguments, main, validate_inputs, parse_oligo_sets


class TestCompareCommand(unittest.TestCase):
    """Test compare command functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.reference_file = self.test_dir / "reference.fasta"
        self.output_dir = self.test_dir / "output"
        
        # Create mock reference file
        with open(self.reference_file, 'w') as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCG" * 100 + "\n")
        
        # Create mock oligo files
        self.oligo_files = {}
        for i, name in enumerate(['set1', 'set2', 'set3']):
            oligo_file = self.test_dir / f"{name}.fasta"
            with open(oligo_file, 'w') as f:
                f.write(f">{name}_oligo_1\n")
                f.write("ATCGATCGATCG" * 10 + "\n")
                f.write(f">{name}_oligo_2\n") 
                f.write("GCTAGCTAGCTA" * 10 + "\n")
            self.oligo_files[name] = str(oligo_file)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_add_arguments(self):
        """Test argument parser setup."""
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        
        # Test required arguments
        self.assertTrue(any(action.dest == 'reference' for action in parser._actions))
        self.assertTrue(any(action.dest == 'output' for action in parser._actions))
        self.assertTrue(any(action.dest == 'sets' for action in parser._actions))
        
        # Test optional arguments
        self.assertTrue(any(action.dest == 'min_identity' for action in parser._actions))
        self.assertTrue(any(action.dest == 'enable_statistical_analysis' for action in parser._actions))
        self.assertTrue(any(action.dest == 'plot_format' for action in parser._actions))
    
    def test_parse_oligo_sets(self):
        """Test oligo set specification parsing."""
        # Valid specifications
        valid_specs = [
            f"Set1:{self.oligo_files['set1']}",
            f"Set2:{self.oligo_files['set2']}",
            f"Set3:{self.oligo_files['set3']}"
        ]
        
        result = parse_oligo_sets(valid_specs)
        
        self.assertEqual(len(result), 3)
        self.assertIn('Set1', result)
        self.assertIn('Set2', result)
        self.assertIn('Set3', result)
        self.assertEqual(result['Set1'], self.oligo_files['set1'])
    
    def test_parse_oligo_sets_duplicate_names(self):
        """Test handling of duplicate oligo set names."""
        duplicate_specs = [
            f"Set1:{self.oligo_files['set1']}",
            f"Set1:{self.oligo_files['set2']}"  # Duplicate name
        ]
        
        with self.assertRaises(SystemExit):
            parse_oligo_sets(duplicate_specs)
    
    def test_validate_inputs_valid(self):
        """Test input validation with valid inputs."""
        args = MagicMock()
        args.reference = str(self.reference_file)
        args.sets = [
            f"Set1:{self.oligo_files['set1']}",
            f"Set2:{self.oligo_files['set2']}"
        ]
        args.min_identity = 90.0
        args.min_coverage = 1.0
        args.target_coverage = 10.0
        args.significance_level = 0.05
        
        # Should not raise exception
        with patch('subprocess.run'):
            validate_inputs(args)
    
    def test_validate_inputs_missing_reference(self):
        """Test validation with missing reference file."""
        args = MagicMock()
        args.reference = str(self.test_dir / "nonexistent.fasta")
        args.sets = [f"Set1:{self.oligo_files['set1']}"]
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_validate_inputs_insufficient_sets(self):
        """Test validation with insufficient oligo sets."""
        args = MagicMock()
        args.reference = str(self.reference_file)
        args.sets = [f"Set1:{self.oligo_files['set1']}"]  # Only one set
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_validate_inputs_invalid_set_format(self):
        """Test validation with invalid oligo set format."""
        args = MagicMock()
        args.reference = str(self.reference_file)
        args.sets = [
            f"Set1:{self.oligo_files['set1']}",
            "InvalidFormat"  # Missing colon
        ]
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_validate_inputs_missing_oligo_file(self):
        """Test validation with missing oligo file."""
        args = MagicMock()
        args.reference = str(self.reference_file)
        args.sets = [
            f"Set1:{self.oligo_files['set1']}",
            "Set2:nonexistent.fasta"
        ]
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    
    def test_validate_inputs_invalid_parameters(self):
        """Test validation with invalid parameter values."""
        args = MagicMock()
        args.reference = str(self.reference_file)
        args.sets = [
            f"Set1:{self.oligo_files['set1']}",
            f"Set2:{self.oligo_files['set2']}"
        ]
        
        # Test invalid identity
        args.min_identity = 150.0  # Invalid
        args.min_coverage = 1.0
        args.target_coverage = 10.0
        args.significance_level = 0.05
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
        
        # Test invalid coverage
        args.min_identity = 90.0
        args.min_coverage = -1.0  # Invalid
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
        
        # Test invalid significance level
        args.min_coverage = 1.0
        args.significance_level = 1.5  # Invalid
        
        with self.assertRaises(SystemExit):
            validate_inputs(args)
    

class TestCompareCommandArguments(unittest.TestCase):
    """Test compare command argument parsing."""
    
    def test_argument_defaults(self):
        """Test default argument values."""
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        
        # Parse minimal arguments
        args = parser.parse_args([
            '--reference', 'ref.fasta',
            '--output', 'output/',
            '--sets', 'Set1:set1.fasta', 'Set2:set2.fasta'
        ])
        
        # Check defaults
        self.assertEqual(args.min_identity, 90.0)
        self.assertEqual(args.min_length, 100)
        self.assertEqual(args.min_coverage, 1.0)
        self.assertEqual(args.target_coverage, 10.0)
        self.assertEqual(args.threads, 1)
        self.assertTrue(args.enable_statistical_analysis)
        self.assertEqual(args.significance_level, 0.05)
        self.assertEqual(args.multiple_comparison_correction, 'fdr')
        self.assertEqual(args.plot_format, 'png')
        self.assertEqual(args.plot_dpi, 300)
        self.assertFalse(args.keep_intermediates)
        self.assertFalse(args.quiet)
        self.assertFalse(args.log)
    
    def test_argument_overrides(self):
        """Test argument value overrides."""
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        
        # Parse with custom values
        args = parser.parse_args([
            '--reference', 'ref.fasta',
            '--output', 'output/',
            '--sets', 'Set1:set1.fasta', 'Set2:set2.fasta',
            '--min-identity', '95.0',
            '--threads', '4',
            '--significance-level', '0.01',
            '--multiple-comparison-correction', 'bonferroni',
            '--plot-format', 'pdf',
            '--keep-intermediates',
            '--quiet'
        ])
        
        # Check overridden values
        self.assertEqual(args.min_identity, 95.0)
        self.assertEqual(args.threads, 4)
        self.assertEqual(args.significance_level, 0.01)
        self.assertEqual(args.multiple_comparison_correction, 'bonferroni')
        self.assertEqual(args.plot_format, 'pdf')
        self.assertTrue(args.keep_intermediates)
        self.assertTrue(args.quiet)
    
    def test_multiple_sets_parsing(self):
        """Test parsing multiple oligo sets."""
        parser = argparse.ArgumentParser()
        add_arguments(parser)
        
        args = parser.parse_args([
            '--reference', 'ref.fasta',
            '--output', 'output/',
            '--sets', 'Set1:set1.fasta', 'Set2:set2.fasta', 'Set3:set3.fasta', 'Set4:set4.fasta'
        ])
        
        self.assertEqual(len(args.sets), 4)
        self.assertEqual(args.sets[0], 'Set1:set1.fasta')
        self.assertEqual(args.sets[1], 'Set2:set2.fasta')
        self.assertEqual(args.sets[2], 'Set3:set3.fasta')
        self.assertEqual(args.sets[3], 'Set4:set4.fasta')


if __name__ == '__main__':
    unittest.main()