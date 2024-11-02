#!/usr/bin/env python3

"""
plot.py

A script to generate statistical plots from the parameters file produced by stats.py.

Features:
- Generate plots for selected columns or all available columns.
- Color plots based on the 'Kept' column to distinguish between kept and filtered baits.
- Support for various plot types including histogram, boxplot, scatterplot, violinplot, densityplot, heatmap, pairplot, swarmplot, jointplot, boxenplot, and PCA.
- For scatterplots and jointplots, generate all pairwise combinations of selected columns if two or more are provided.
- Warn if the number of scatterplot or jointplot combinations exceeds 100 to prevent excessive plot generation.
- Generate PCA plots to visualize data in reduced dimensions, including both scores and loadings plots.
- Automatic handling of compressed (gzip) and uncompressed parameters files.

Usage:
    baitUtils plot -i params.txt -o plots/ --columns GC% Tm MFE --plot_type histogram boxplot scatterplot --color Kept
"""

import argparse
from baitUtils._version import __version__  # Import version from _version.py
import logging
import os
import sys
import gzip
from typing import List, Optional, Union, Dict, Any, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

def add_arguments(parser):
    """
    Add command-line arguments specific to the 'plot' subcommand.

    Args:
        parser (argparse.ArgumentParser): The argument parser to which arguments will be added.
    """
    parser.add_argument('--version', action='version', version=f"baitUtils plot {__version__}")
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input parameters file (can be .txt or .txt.gz)'
    )
    parser.add_argument(
        '-o', '--outdir',
        default='plots',
        help='Output directory for plots (Default: plots)'
    )
    parser.add_argument(
        '--columns',
        nargs='+',
        help='List of columns to plot. If not specified, all numeric columns will be plotted.'
    )
    parser.add_argument(
        '--plot_type',
        nargs='+',
        choices=['histogram', 'boxplot', 'scatterplot', 'violinplot', 'densityplot', 'heatmap', 'pairplot', 'swarmplot', 'jointplot', 'boxenplot', 'pca'],
        default=['histogram'],
        help='Type of plots to generate. Choose from histogram, boxplot, scatterplot, violinplot, densityplot, heatmap, pairplot, swarmplot, jointplot, boxenplot, pca. (Default: histogram)'
    )
    parser.add_argument(
        '--color',
        choices=['Kept'],
        help='Column to use for coloring the plots. Currently supports "Kept".'
    )
    parser.add_argument(
        '--format',
        choices=['png', 'pdf', 'svg'],
        default='png',
        help='Output file format for the plots (Default: png)'
    )
    parser.add_argument(
        '-l', '--log',
        action='store_true',
        help='Enable detailed logging'
    )

def main(args):
    """
    Main function for the 'plot' subcommand.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    setup_logging(args.log)
    create_output_directory(args.outdir)
    df = open_params_file(args.input)
    # If coloring is based on 'Kept', ensure the column exists
    if args.color:
        if args.color not in df.columns:
            logging.error(f"Color column '{args.color}' not found in parameters file.")
            sys.exit(1)
        if not pd.api.types.is_categorical_dtype(df[args.color]) and not pd.api.types.is_object_dtype(df[args.color]):
            logging.warning(f"Color column '{args.color}' is not categorical. Converting to category.")
            df[args.color] = df[args.color].astype('category')
    # Validate and select columns to plot
    selected_columns = validate_columns(df, args.columns, args.plot_type)
    # Generate plots
    plot_statistics(
        df=df,
        columns=selected_columns,
        plot_types=args.plot_type,
        color_col=args.color,
        outdir=args.outdir,
        file_format=args.format
    )
    logging.info("All plots have been generated successfully.")

def setup_logging(enable_debug: bool) -> None:
    """
    Configure logging settings.

    Args:
        enable_debug (bool): If True, set logging level to DEBUG. Otherwise, INFO.
    """
    log_level = logging.DEBUG if enable_debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def open_params_file(filename: str) -> pd.DataFrame:
    """
    Open and read the parameters file into a pandas DataFrame.
    Supports both uncompressed and gzipped files.

    Args:
        filename (str): Path to the parameters file.

    Returns:
        pd.DataFrame: DataFrame containing the parameters.
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
        sys.exit(1)

def validate_columns(df: pd.DataFrame, columns: Optional[List[str]], plot_types: List[str]) -> List[str]:
    """
    Validate and return the list of columns to plot.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        columns (Optional[List[str]]): User-specified columns to plot.
        plot_types (List[str]): List of plot types selected.

    Returns:
        List[str]: Validated list of columns to plot.
    """
    if columns:
        missing_cols = [col for col in columns if col not in df.columns]
        if missing_cols:
            logging.error(f"Specified columns not found in parameters file: {missing_cols}")
            sys.exit(1)
        # Select only numeric columns
        numeric_cols = df[columns].select_dtypes(include='number').columns.tolist()
        if len(numeric_cols) < len(columns):
            non_numeric = set(columns) - set(numeric_cols)
            logging.warning(f"Non-numeric columns will be skipped for plotting: {non_numeric}")
        validated_cols = numeric_cols
    else:
        # If no columns specified, select all numeric columns excluding 'Kept' if used for coloring
        validated_cols = df.select_dtypes(include='number').columns.tolist()
        if 'Kept' in validated_cols:
            validated_cols.remove('Kept')
    # If scatterplot, jointplot, PCA, or heatmap is selected, ensure at least two columns are present
    required_plots = {'scatterplot', 'jointplot', 'pca', 'heatmap', 'pairplot'}
    if any(plot in plot_types for plot in required_plots):
        if len(validated_cols) < 2:
            logging.error("At least two numeric columns are required for scatterplots, jointplots, PCA, heatmaps, and pairplots.")
            sys.exit(1)
    if not validated_cols:
        logging.error("No valid numeric columns available for plotting.")
        sys.exit(1)
    logging.info(f"Columns to plot: {validated_cols}")
    return validated_cols

def create_output_directory(directory: str) -> None:
    """
    Create the output directory if it doesn't exist.

    Args:
        directory (str): Path to the output directory.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info(f"Created output directory: {directory}")
        except Exception as e:
            logging.error(f"Failed to create output directory '{directory}': {e}")
            sys.exit(1)
    else:
        logging.info(f"Output directory already exists: {directory}")

def generate_histogram(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a histogram for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.histplot(data=df, x=column, hue=color_col, multiple="stack", palette="Set1", kde=True)
    else:
        sns.histplot(data=df, x=column, kde=True, color='blue')
    plt.title(f'Histogram of {column}')
    plt.xlabel(column)
    plt.ylabel('Frequency')
    plt.tight_layout()
    filename = os.path.join(outdir, f'histogram_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved histogram: {filename}")

def generate_boxplot(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a boxplot for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.boxplot(data=df, y=column, x=color_col, palette="Set1")
        plt.xlabel(color_col)
    else:
        sns.boxplot(data=df, y=column, color='green')
        plt.xlabel('')
    plt.title(f'Boxplot of {column}')
    plt.ylabel(column)
    plt.tight_layout()
    filename = os.path.join(outdir, f'boxplot_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved boxplot: {filename}")

def generate_violinplot(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a violin plot for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.violinplot(data=df, y=column, x=color_col, palette="Set1", inner="quartile")
        plt.xlabel(color_col)
    else:
        sns.violinplot(data=df, y=column, color='orange', inner="quartile")
        plt.xlabel('')
    plt.title(f'Violin Plot of {column}')
    plt.ylabel(column)
    plt.tight_layout()
    filename = os.path.join(outdir, f'violinplot_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved violin plot: {filename}")

def generate_densityplot(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a density plot for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.kdeplot(data=df, x=column, hue=color_col, palette="Set1", fill=True)
    else:
        sns.kdeplot(data=df, x=column, fill=True, color='red')
    plt.title(f'Density Plot of {column}')
    plt.xlabel(column)
    plt.ylabel('Density')
    plt.tight_layout()
    filename = os.path.join(outdir, f'densityplot_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved density plot: {filename}")

def generate_boxenplot(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a boxen plot for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.boxenplot(data=df, y=column, x=color_col, palette="Set1")
        plt.xlabel(color_col)
    else:
        sns.boxenplot(data=df, y=column, color='brown')
        plt.xlabel('')
    plt.title(f'Boxen Plot of {column}')
    plt.ylabel(column)
    plt.tight_layout()
    filename = os.path.join(outdir, f'boxenplot_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved boxen plot: {filename}")

def generate_swarmplot(df: pd.DataFrame, column: str, color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a swarm plot for the specified column.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        column (str): Column to plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.swarmplot(data=df, y=column, x=color_col, palette="Set1")
        plt.xlabel(color_col)
    else:
        sns.swarmplot(data=df, y=column, color='cyan')
        plt.xlabel('')
    plt.title(f'Swarm Plot of {column}')
    plt.ylabel(column)
    plt.tight_layout()
    filename = os.path.join(outdir, f'swarmplot_{column}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved swarm plot: {filename}")

def generate_scatterplot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: Optional[str],
    outdir: str,
    file_format: str
) -> None:
    """
    Generate and save a scatterplot for the specified x and y columns.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        x_col (str): Column to plot on the x-axis.
        y_col (str): Column to plot on the y-axis.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(10, 6))
    if color_col:
        sns.scatterplot(data=df, x=x_col, y=y_col, hue=color_col, palette="Set1")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
    else:
        sns.scatterplot(data=df, x=x_col, y=y_col, color='purple')
        plt.xlabel(x_col)
        plt.ylabel(y_col)
    plt.title(f'Scatterplot of {x_col} vs {y_col}')
    plt.tight_layout()
    filename = os.path.join(outdir, f'scatterplot_{x_col}_vs_{y_col}.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved scatterplot: {filename}")

def generate_heatmap(
    df: pd.DataFrame,
    columns: List[str],
    color_col: Optional[str],
    outdir: str,
    file_format: str
) -> None:
    """
    Generate and save a correlation heatmap for the specified columns.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        columns (List[str]): List of columns to include in the heatmap.
        color_col (Optional[str]): Not used in heatmap but kept for consistency.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    plt.figure(figsize=(12, 10))
    corr_matrix = df[columns].corr()
    sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', square=True)
    plt.title('Correlation Heatmap')
    plt.tight_layout()
    filename = os.path.join(outdir, f'heatmap_correlation.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved correlation heatmap: {filename}")

def generate_pairplot(df: pd.DataFrame, columns: List[str], color_col: Optional[str], outdir: str, file_format: str) -> None:
    """
    Generate and save a pair plot for the specified columns.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        columns (List[str]): List of columns to include in the pair plot.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    sns.pairplot(df, vars=columns, hue=color_col, palette="Set1" if color_col else None)
    plt.suptitle('Pair Plot', y=1.02)
    plt.tight_layout()
    filename = os.path.join(outdir, f'pairplot.{file_format}')
    plt.savefig(filename)
    plt.close()
    logging.info(f"Saved pair plot: {filename}")

def generate_jointplot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: Optional[str],
    outdir: str,
    file_format: str
) -> None:
    """
    Generate and save a joint plot for the specified x and y columns.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        x_col (str): Column to plot on the x-axis.
        y_col (str): Column to plot on the y-axis.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plot.
    """
    if color_col:
        g = sns.jointplot(data=df, x=x_col, y=y_col, hue=color_col, palette="Set1", kind='scatter')
    else:
        g = sns.jointplot(data=df, x=x_col, y=y_col, kind='scatter', color='magenta')
    g.fig.suptitle(f'Joint Plot of {x_col} vs {y_col}', y=1.02)
    plt.tight_layout()
    filename = os.path.join(outdir, f'jointplot_{x_col}_vs_{y_col}.{file_format}')
    g.savefig(filename)
    plt.close()
    logging.info(f"Saved joint plot: {filename}")

def generate_pca(
    df: pd.DataFrame,
    columns: List[str],
    color_col: Optional[str],
    outdir: str,
    file_format: str
) -> None:
    """
    Generate and save PCA plots (scores and loadings) for the specified columns.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        columns (List[str]): List of columns to include in PCA.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plots.
    """
    if len(columns) < 2:
        logging.error("At least two columns are required for PCA.")
        sys.exit(1)
    
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[columns])
    
    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    
    # Create a DataFrame for loadings
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loadings_df = pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=columns)
    
    # Concatenate with the original DataFrame for coloring
    if color_col:
        pca_df = pd.concat([pca_df, df[color_col].reset_index(drop=True)], axis=1)
    
    # Plot PCA Scores
    plt.figure(figsize=(10, 8))
    if color_col:
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=color_col, palette="Set1")
    else:
        sns.scatterplot(data=pca_df, x='PC1', y='PC2', color='darkorange')
    plt.title('PCA Scores Plot')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')
    plt.legend(title=color_col if color_col else None)
    plt.tight_layout()
    filename_scores = os.path.join(outdir, f'pca_scores_plot.{file_format}')
    plt.savefig(filename_scores)
    plt.close()
    logging.info(f"Saved PCA scores plot: {filename_scores}")
    
    # Plot PCA Loadings
    plt.figure(figsize=(10, 8))
    for i, (col, pc1, pc2) in enumerate(zip(loadings_df.index, loadings_df['PC1'], loadings_df['PC2'])):
        plt.arrow(0, 0, pc1, pc2, color='r', alpha=0.5)
        plt.text(pc1*1.15, pc2*1.15, col, color='g', ha='center', va='center')
    plt.xlabel('PC1 Loadings')
    plt.ylabel('PC2 Loadings')
    plt.title('PCA Loadings Plot')
    plt.grid()
    plt.tight_layout()
    filename_loadings = os.path.join(outdir, f'pca_loadings_plot.{file_format}')
    plt.savefig(filename_loadings)
    plt.close()
    logging.info(f"Saved PCA loadings plot: {filename_loadings}")

def plot_statistics(
    df: pd.DataFrame,
    columns: List[str],
    plot_types: List[str],
    color_col: Optional[str],
    outdir: str,
    file_format: str
) -> None:
    """
    Generate and save plots for the specified columns and plot types.

    Args:
        df (pd.DataFrame): DataFrame containing the parameters.
        columns (List[str]): List of columns to plot.
        plot_types (List[str]): List of plot types to generate.
        color_col (Optional[str]): Column to use for coloring.
        outdir (str): Output directory for plots.
        file_format (str): File format for the plots.
    """
    # Generate histograms, boxplots, violinplots, densityplots, swarmplots, boxenplots
    for column in columns:
        for plot_type in plot_types:
            if plot_type == 'histogram':
                generate_histogram(df, column, color_col, outdir, file_format)
            elif plot_type == 'boxplot':
                generate_boxplot(df, column, color_col, outdir, file_format)
            elif plot_type == 'violinplot':
                generate_violinplot(df, column, color_col, outdir, file_format)
            elif plot_type == 'densityplot':
                generate_densityplot(df, column, color_col, outdir, file_format)
            elif plot_type == 'swarmplot':
                generate_swarmplot(df, column, color_col, outdir, file_format)
            elif plot_type == 'boxenplot':
                generate_boxenplot(df, column, color_col, outdir, file_format)
    
    # Handle scatterplots separately to accommodate pairwise combinations
    if 'scatterplot' in plot_types:
        # Determine which columns to use for scatterplots
        scatter_cols = columns if columns else df.select_dtypes(include='number').columns.tolist()
        if len(scatter_cols) < 2:
            logging.error("At least two numeric columns are required for scatterplots.")
            sys.exit(1)
        # Generate all pairwise combinations
        scatter_pairs = list(combinations(scatter_cols, 2))
        num_pairs = len(scatter_pairs)
        if num_pairs > 100:
            logging.warning(f"Number of scatterplot combinations ({num_pairs}) exceeds 100. Only the first 100 will be plotted.")
            scatter_pairs = scatter_pairs[:100]
        for x_col, y_col in scatter_pairs:
            generate_scatterplot(df, x_col, y_col, color_col, outdir, file_format)
    
    # Handle jointplots separately to accommodate pairwise combinations
    if 'jointplot' in plot_types:
        # Determine which columns to use for jointplots
        joint_cols = columns if columns else df.select_dtypes(include='number').columns.tolist()
        if len(joint_cols) < 2:
            logging.error("At least two numeric columns are required for jointplots.")
            sys.exit(1)
        # Generate all pairwise combinations
        joint_pairs = list(combinations(joint_cols, 2))
        num_joint_pairs = len(joint_pairs)
        if num_joint_pairs > 100:
            logging.warning(f"Number of jointplot combinations ({num_joint_pairs}) exceeds 100. Only the first 100 will be plotted.")
            joint_pairs = joint_pairs[:100]
        for x_col, y_col in joint_pairs:
            generate_jointplot(df, x_col, y_col, color_col, outdir, file_format)
    
    # Handle heatmaps if selected
    if 'heatmap' in plot_types:
        # Heatmaps are generally useful for correlation matrices
        # We'll assume all numeric columns are relevant
        heatmap_cols = columns if columns else df.select_dtypes(include='number').columns.tolist()
        if len(heatmap_cols) < 2:
            logging.error("At least two numeric columns are required for heatmaps.")
            sys.exit(1)
        generate_heatmap(df, heatmap_cols, color_col, outdir, file_format)
    
    # Handle pairplots if selected
    if 'pairplot' in plot_types:
        generate_pairplot(df, columns, color_col, outdir, file_format)
    
    # Handle PCA if selected
    if 'pca' in plot_types:
        generate_pca(df, columns, color_col, outdir, file_format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate statistical plots from bait parameters file.')
    add_arguments(parser)
    args = parser.parse_args()
    main(args)