import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def process_gene_counts(df):
    
    """
    
    """
    
    gene_grouped = df.groupby('gene')
    gene_summary = gene_grouped.describe()
    
    return gene_summary

def process_cell_counts(df):
    
    """
    
    """
    
    cell_grouped = df.groupby('cell')
    cell_summary = cell_grouped.describe()
    cell_summary['umi_counts_per_cell'] = cell_summary['count']['count'].values*cell_summary['count']['mean'].values
    cell_summary['log_umi'] = np.log10(cell_summary['umi_counts_per_cell'])
    cell_summary = cell_summary.sort_values('log_umi', ascending=False)
    
    return cell_summary

def drop_N_barcodes(df):

    """Identify cell barcodes containing an N and drop"""

    N_idx = df.loc[df.cell.str.contains("N")].index.astype(int)
    df = df.drop(index=N_idx)

    return df

def _summarize_experiment(counts_path):
    
    df = pd.read_csv(counts_path, sep="\t", header=0)
    
    cell_summary = process_cell_counts(df)
    gene_summary = process_gene_counts(df)
    
    return df, cell_summary, gene_summary

def _filter_counts(df, cell_summary, gene_summary, min_genes=200, min_cells=3, save_path=None):
    
    """
    
    """

    num_genes_unfiltered = df.gene.nunique()
    num_cells_unfiltered = df.cell.nunique()

    filtered_genes = gene_summary['count'].loc[gene_summary['count']['count'] > min_cells]
    filtered_cells = cell_summary['count'].loc[cell_summary['count']['count'] > min_genes ]
    
    df_filtered = df.loc[df.cell.isin(filtered_cells.index.values)].loc[df.gene.isin(filtered_genes.index.values)]
    
    num_genes = df_filtered.gene.nunique()
    num_cells = df_filtered.cell.nunique()
    
    gene_ratio = (num_genes/num_genes_unfiltered)*100
    cell_ratio = (num_cells/num_cells_unfiltered)*100
    
    print((num_cells_unfiltered - num_cells), "cells filtered. Remaining:", num_cells, "(",np.round(cell_ratio, 2), "%)")
    print((num_genes_unfiltered - num_genes), "genes filtered. Remaining:", num_genes, "(",np.round(gene_ratio, 2), "%)")
    
    df_filtered = drop_N_barcodes(df_filtered)
    
    if save_path != None:
        df_filtered.to_csv(save_path, index=False)
    
    return df, df_filtered