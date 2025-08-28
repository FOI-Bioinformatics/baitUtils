#!/usr/bin/env python3

"""
plotting_utils.py

Core plotting utilities for generating statistical visualizations.
Provides various plot types including histograms, boxplots, scatterplots,
correlation heatmaps, PCA plots, and more.
"""

import logging
import os
from typing import List, Optional
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class PlotGenerator:
    """Generates various types of statistical plots for sequence data."""
    
    def __init__(self, figure_size: tuple = (10, 6)):
        """
        Initialize plot generator with default figure size.
        
        Args:
            figure_size: Default figure size for plots
        """
        self.figure_size = figure_size
        
    def generate_histogram(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                          outdir: str, file_format: str) -> None:
        """Generate and save a histogram for the specified column."""
        plt.figure(figsize=self.figure_size)
        if color_col:
            sns.histplot(data=df, x=column, hue=color_col, multiple="stack", 
                        palette="Set1", kde=True)
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

    def generate_boxplot(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                        outdir: str, file_format: str) -> None:
        """Generate and save a boxplot for the specified column."""
        plt.figure(figsize=self.figure_size)
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

    def generate_violinplot(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                           outdir: str, file_format: str) -> None:
        """Generate and save a violin plot for the specified column."""
        plt.figure(figsize=self.figure_size)
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

    def generate_densityplot(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                            outdir: str, file_format: str) -> None:
        """Generate and save a density plot for the specified column."""
        plt.figure(figsize=self.figure_size)
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

    def generate_boxenplot(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                          outdir: str, file_format: str) -> None:
        """Generate and save a boxen plot for the specified column."""
        plt.figure(figsize=self.figure_size)
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

    def generate_swarmplot(self, df: pd.DataFrame, column: str, color_col: Optional[str], 
                          outdir: str, file_format: str) -> None:
        """Generate and save a swarm plot for the specified column."""
        plt.figure(figsize=self.figure_size)
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

    def generate_scatterplot(self, df: pd.DataFrame, x_col: str, y_col: str, 
                           color_col: Optional[str], outdir: str, file_format: str) -> None:
        """Generate and save a scatterplot for the specified x and y columns."""
        plt.figure(figsize=self.figure_size)
        if color_col:
            sns.scatterplot(data=df, x=x_col, y=y_col, hue=color_col, palette="Set1")
        else:
            sns.scatterplot(data=df, x=x_col, y=y_col, color='purple')
        plt.title(f'Scatterplot of {x_col} vs {y_col}')
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.tight_layout()
        filename = os.path.join(outdir, f'scatterplot_{x_col}_vs_{y_col}.{file_format}')
        plt.savefig(filename)
        plt.close()
        logging.info(f"Saved scatterplot: {filename}")

    def generate_jointplot(self, df: pd.DataFrame, x_col: str, y_col: str, 
                          color_col: Optional[str], outdir: str, file_format: str) -> None:
        """Generate and save a joint plot for the specified x and y columns."""
        if color_col:
            g = sns.jointplot(data=df, x=x_col, y=y_col, hue=color_col, 
                             palette="Set1", kind='scatter')
        else:
            g = sns.jointplot(data=df, x=x_col, y=y_col, kind='scatter', color='magenta')
        g.fig.suptitle(f'Joint Plot of {x_col} vs {y_col}', y=1.02)
        plt.tight_layout()
        filename = os.path.join(outdir, f'jointplot_{x_col}_vs_{y_col}.{file_format}')
        g.savefig(filename)
        plt.close()
        logging.info(f"Saved joint plot: {filename}")

    def generate_heatmap(self, df: pd.DataFrame, columns: List[str], 
                        color_col: Optional[str], outdir: str, file_format: str) -> None:
        """Generate and save a correlation heatmap for the specified columns."""
        plt.figure(figsize=(12, 10))
        corr_matrix = df[columns].corr()
        sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap='coolwarm', square=True)
        plt.title('Correlation Heatmap')
        plt.tight_layout()
        filename = os.path.join(outdir, f'heatmap_correlation.{file_format}')
        plt.savefig(filename)
        plt.close()
        logging.info(f"Saved correlation heatmap: {filename}")

    def generate_pairplot(self, df: pd.DataFrame, columns: List[str], 
                         color_col: Optional[str], outdir: str, file_format: str) -> None:
        """Generate and save a pair plot for the specified columns."""
        sns.pairplot(df, vars=columns, hue=color_col, 
                    palette="Set1" if color_col else None)
        plt.suptitle('Pair Plot', y=1.02)
        plt.tight_layout()
        filename = os.path.join(outdir, f'pairplot.{file_format}')
        plt.savefig(filename)
        plt.close()
        logging.info(f"Saved pair plot: {filename}")

    def generate_pca(self, df: pd.DataFrame, columns: List[str], 
                    color_col: Optional[str], outdir: str, file_format: str) -> None:
        """Generate and save PCA plots (scores and loadings) for the specified columns."""
        if len(columns) < 2:
            logging.error("At least two columns are required for PCA.")
            return
        
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
        for i, (col, pc1, pc2) in enumerate(zip(loadings_df.index, 
                                               loadings_df['PC1'], loadings_df['PC2'])):
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


class StatisticalPlotter:
    """High-level interface for generating statistical plots."""
    
    def __init__(self):
        """Initialize statistical plotter."""
        self.plot_generator = PlotGenerator()
    
    def plot_statistics(self, df: pd.DataFrame, columns: List[str], plot_types: List[str],
                       color_col: Optional[str], outdir: str, file_format: str) -> None:
        """
        Generate and save plots for the specified columns and plot types.
        
        Args:
            df: DataFrame containing the data
            columns: List of columns to plot
            plot_types: List of plot types to generate
            color_col: Column to use for coloring
            outdir: Output directory for plots
            file_format: File format for the plots
        """
        # Generate single-column plots
        single_column_plots = ['histogram', 'boxplot', 'violinplot', 'densityplot', 
                              'swarmplot', 'boxenplot']
        
        for column in columns:
            for plot_type in plot_types:
                if plot_type in single_column_plots:
                    self._generate_single_column_plot(
                        df, column, plot_type, color_col, outdir, file_format
                    )
        
        # Handle multi-column plots
        self._generate_multi_column_plots(df, columns, plot_types, color_col, 
                                        outdir, file_format)
    
    def _generate_single_column_plot(self, df: pd.DataFrame, column: str, plot_type: str,
                                   color_col: Optional[str], outdir: str, 
                                   file_format: str) -> None:
        """Generate plots that require only one column."""
        plot_methods = {
            'histogram': self.plot_generator.generate_histogram,
            'boxplot': self.plot_generator.generate_boxplot,
            'violinplot': self.plot_generator.generate_violinplot,
            'densityplot': self.plot_generator.generate_densityplot,
            'swarmplot': self.plot_generator.generate_swarmplot,
            'boxenplot': self.plot_generator.generate_boxenplot
        }
        
        if plot_type in plot_methods:
            plot_methods[plot_type](df, column, color_col, outdir, file_format)
    
    def _generate_multi_column_plots(self, df: pd.DataFrame, columns: List[str], 
                                   plot_types: List[str], color_col: Optional[str],
                                   outdir: str, file_format: str) -> None:
        """Generate plots that require multiple columns."""
        # Handle scatterplots
        if 'scatterplot' in plot_types:
            self._generate_pairwise_plots(df, columns, 'scatter', color_col, 
                                        outdir, file_format)
        
        # Handle jointplots
        if 'jointplot' in plot_types:
            self._generate_pairwise_plots(df, columns, 'joint', color_col, 
                                        outdir, file_format)
        
        # Handle heatmaps
        if 'heatmap' in plot_types:
            if len(columns) >= 2:
                self.plot_generator.generate_heatmap(df, columns, color_col, 
                                                   outdir, file_format)
        
        # Handle pairplots
        if 'pairplot' in plot_types:
            if len(columns) >= 2:
                self.plot_generator.generate_pairplot(df, columns, color_col, 
                                                    outdir, file_format)
        
        # Handle PCA
        if 'pca' in plot_types:
            if len(columns) >= 2:
                self.plot_generator.generate_pca(df, columns, color_col, 
                                               outdir, file_format)
    
    def _generate_pairwise_plots(self, df: pd.DataFrame, columns: List[str], 
                               plot_type: str, color_col: Optional[str],
                               outdir: str, file_format: str) -> None:
        """Generate pairwise combination plots (scatter or joint plots)."""
        if len(columns) < 2:
            logging.error(f"At least two numeric columns are required for {plot_type}plots.")
            return
        
        # Generate all pairwise combinations
        pairs = list(combinations(columns, 2))
        num_pairs = len(pairs)
        
        if num_pairs > 100:
            logging.warning(f"Number of {plot_type}plot combinations ({num_pairs}) "
                          f"exceeds 100. Only the first 100 will be plotted.")
            pairs = pairs[:100]
        
        for x_col, y_col in pairs:
            if plot_type == 'scatter':
                self.plot_generator.generate_scatterplot(df, x_col, y_col, color_col, 
                                                       outdir, file_format)
            elif plot_type == 'joint':
                self.plot_generator.generate_jointplot(df, x_col, y_col, color_col, 
                                                     outdir, file_format)


def validate_columns(df: pd.DataFrame, columns: Optional[List[str]], 
                    plot_types: List[str]) -> List[str]:
    """
    Validate and return the list of columns to plot.
    
    Args:
        df: DataFrame containing the data
        columns: User-specified columns to plot
        plot_types: List of plot types selected
        
    Returns:
        Validated list of columns to plot
    """
    if columns:
        missing_cols = [col for col in columns if col not in df.columns]
        if missing_cols:
            logging.error(f"Specified columns not found in parameters file: {missing_cols}")
            raise ValueError(f"Missing columns: {missing_cols}")
        
        # Select only numeric columns
        numeric_cols = df[columns].select_dtypes(include='number').columns.tolist()
        if len(numeric_cols) < len(columns):
            non_numeric = set(columns) - set(numeric_cols)
            logging.warning(f"Non-numeric columns will be skipped: {non_numeric}")
        validated_cols = numeric_cols
    else:
        # Select all numeric columns excluding 'Kept' if used for coloring
        validated_cols = df.select_dtypes(include='number').columns.tolist()
        if 'Kept' in validated_cols:
            validated_cols.remove('Kept')
    
    # Check if we have enough columns for multi-column plots
    required_plots = {'scatterplot', 'jointplot', 'pca', 'heatmap', 'pairplot'}
    if any(plot in plot_types for plot in required_plots):
        if len(validated_cols) < 2:
            logging.error("At least two numeric columns are required for "
                         "scatterplots, jointplots, PCA, heatmaps, and pairplots.")
            raise ValueError("Insufficient columns for multi-column plots")
    
    if not validated_cols:
        logging.error("No valid numeric columns available for plotting.")
        raise ValueError("No valid numeric columns")
    
    logging.info(f"Columns to plot: {validated_cols}")
    return validated_cols