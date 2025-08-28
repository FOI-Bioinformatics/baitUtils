#!/usr/bin/env python3

"""
statistical_plots.py

Main command for generating statistical plots from sequence statistics files.
Refactored from plot.py for better organization and maintainability.
"""

import argparse
import logging
import os
import sys
import gzip
from pathlib import Path
from typing import List, Optional

import pandas as pd

from baitUtils._version import __version__
from baitUtils.plotting_utils import StatisticalPlotter, validate_columns


class PlotFileHandler:
    """Handles file operations for plotting command."""
    
    @staticmethod
    def open_params_file(filename: str) -> pd.DataFrame:
        """
        Open and read the parameters file into a pandas DataFrame.
        Supports both uncompressed and gzipped files.
        
        Args:
            filename: Path to the parameters file
            
        Returns:
            DataFrame containing the parameters
        """
        logging.info(f"Reading parameters file: {filename}")
        try:
            if filename.endswith('.gz'):
                df = pd.read_csv(filename, sep='\t', compression='gzip')
            else:
                df = pd.read_csv(filename, sep='\t')
            logging.info(f"Parameters file contains {len(df)} records.")
            return df
        except Exception as e:
            logging.error(f"Failed to read parameters file: {e}")
            raise
    
    @staticmethod
    def create_output_directory(directory: str) -> None:
        """
        Create the output directory if it doesn't exist.
        
        Args:
            directory: Path to the output directory
        """
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
                logging.info(f"Created output directory: {directory}")
            except Exception as e:
                logging.error(f"Failed to create output directory '{directory}': {e}")
                raise
        else:
            logging.info(f"Output directory already exists: {directory}")


class PlotCommandProcessor:
    """Processes the plot command with all its options."""
    
    def __init__(self):
        """Initialize the plot command processor."""
        self.plotter = StatisticalPlotter()
        self.file_handler = PlotFileHandler()
    
    def process_command(self, args) -> None:
        """
        Process the plot command with given arguments.
        
        Args:
            args: Parsed command-line arguments
        """
        # Set up logging
        self._setup_logging(args.log)
        
        # Create output directory
        self.file_handler.create_output_directory(args.outdir)
        
        # Load data
        df = self.file_handler.open_params_file(args.input)
        
        # Validate color column
        if args.color:
            self._validate_color_column(df, args.color)
        
        # Validate and select columns to plot
        selected_columns = validate_columns(df, args.columns, args.plot_type)
        
        # Generate plots
        self.plotter.plot_statistics(
            df=df,
            columns=selected_columns,
            plot_types=args.plot_type,
            color_col=args.color,
            outdir=args.outdir,
            file_format=args.format
        )
        
        logging.info("All plots have been generated successfully.")
    
    def _setup_logging(self, enable_debug: bool) -> None:
        """Configure logging settings."""
        log_level = logging.DEBUG if enable_debug else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    
    def _validate_color_column(self, df: pd.DataFrame, color_col: str) -> None:
        """Validate and prepare color column for plotting."""
        if color_col not in df.columns:
            logging.error(f"Color column '{color_col}' not found in parameters file.")
            raise ValueError(f"Color column '{color_col}' not found")
        
        if not pd.api.types.is_categorical_dtype(df[color_col]) and \
           not pd.api.types.is_object_dtype(df[color_col]):
            logging.warning(f"Color column '{color_col}' is not categorical. "
                           f"Converting to category.")
            df[color_col] = df[color_col].astype('category')


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add command-line arguments for statistical plotting."""
    parser.add_argument('--version', action='version', 
                       version=f"baitUtils statistical-plots {__version__}")
    
    # Input/Output
    parser.add_argument('-i', '--input', required=True,
                       help='Input parameters file (can be .txt or .txt.gz)')
    parser.add_argument('-o', '--outdir', default='plots',
                       help='Output directory for plots (default: plots)')
    
    # Plot configuration
    parser.add_argument('--columns', nargs='+',
                       help='List of columns to plot. If not specified, '
                            'all numeric columns will be plotted.')
    parser.add_argument('--plot_type', nargs='+',
                       choices=['histogram', 'boxplot', 'scatterplot', 
                               'violinplot', 'densityplot', 'heatmap', 
                               'pairplot', 'swarmplot', 'jointplot', 
                               'boxenplot', 'pca'],
                       default=['histogram'],
                       help='Type of plots to generate (default: histogram)')
    parser.add_argument('--color', choices=['Kept'],
                       help='Column to use for coloring the plots. '
                            'Currently supports "Kept".')
    parser.add_argument('--format', choices=['png', 'pdf', 'svg'], 
                       default='png',
                       help='Output file format for the plots (default: png)')
    
    # Logging
    parser.add_argument('-l', '--log', action='store_true',
                       help='Enable detailed logging')


def main(args) -> None:
    """Main function for statistical plotting command."""
    processor = PlotCommandProcessor()
    try:
        processor.process_command(args)
    except Exception as e:
        logging.error(f"Plot generation failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate statistical plots from sequence statistics file.'
    )
    add_arguments(parser)
    args = parser.parse_args()
    main(args)