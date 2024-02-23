import scanpy as sc
import pandas as pd
import numpy as np
import os
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from collections import Counter
import glob
import alphashape
import geopandas as gpd
import seaborn as sns
from shapely.ops import transform
import imageio as io

def prepare_transcripts(input_file):
    transcripts_cellpose = pd.read_csv(os.path.join(input_file, 'transcripts_cellpose.csv'), index_col=0)
    transcripts = pd.read_csv(os.path.join(input_file, 'transcripts.csv'), index_col=0)
    transcripts = transcripts.dropna(subset=['cell'])
    transcripts_nucleus = transcripts[transcripts['overlaps_nucleus'] == 1]
    transcripts_cellpose.index = transcripts_cellpose.transcript_id.values
    transcripts = transcripts[np.isin(transcripts.index.values, transcripts_cellpose.index.values)]
    cell_values = transcripts.cell
    cell_values.fillna('Not_assigned-0', inplace=True)
    transcripts['cell_number'] = cell_values
    return transcripts, transcripts_cellpose
    
def assign_nuclei_to_cells(transcripts, transcripts_cellpose):
    overlap = transcripts[transcripts.overlaps_nucleus == 1]
    nuclei_associated = transcripts_cellpose.loc[overlap.index.values]
    overlap['associated_nucleus'] = nuclei_associated.cell_id.values 
    cell_numbers = overlap.cell_number.values
    associated_nuclei = overlap.associated_nucleus.values
    
    # Create a dictionary to store the most common value for each unique cell number
    most_common_values = {}
    
    # Iterate through the pairs of cell numbers and associated nuclei
    for cell_number, nucleus in zip(cell_numbers, associated_nuclei):
        if cell_number not in most_common_values:
            most_common_values[cell_number] = Counter()
    
        most_common_values[cell_number][nucleus] += 1
    
    # Calculate the most common nucleus for each unique cell number
    result = {cell_number: counter.most_common(1)[0][0] if counter.most_common(1)[0][1] > 1 else -1 for cell_number, counter in most_common_values.items()}
    return result
    
def find_main_nucleus(transcripts, transcripts_cellpose, result):
    keys = list(result.keys())
    values = list(result.values())
    index = [i for i in range(len(keys))]
    
    keydf = pd.DataFrame(zip(keys, values, index), columns=['cell_number', 'nucleus', 'inds'])
    keydf['nucleus'] = keydf['nucleus'].replace(-1, np.nan)
    
    transcripts_with_gt_nucleus = transcripts.merge(transcripts_cellpose['cell_id'], left_index=True, right_index=True, how='left')
    transcripts_with_gt_and_main_nucleus = transcripts_with_gt_nucleus.merge(keydf, left_on='cell', right_on='cell_number')
    transcripts_with_gt_and_main_nucleus_filtered = transcripts_with_gt_and_main_nucleus[transcripts_with_gt_and_main_nucleus.nucleus != -1]
    transcripts_with_gt_and_main_nucleus_filtered['indexer'] = [i for i in range(len(transcripts_with_gt_and_main_nucleus_filtered.index))]
    groupby_most_common_nucleus = transcripts_with_gt_and_main_nucleus_filtered.groupby('nucleus')
    return transcripts_with_gt_and_main_nucleus_filtered,  groupby_most_common_nucleus

def reassign_multiple_nuclei(transcripts_with_gt_and_main_nucleus_filtered, groupby_most_common_nucleus):
    reassignments = np.zeros(len(transcripts_with_gt_and_main_nucleus_filtered))
    for group_name, group_data in groupby_most_common_nucleus:
        unique = np.unique(group_data.cell_id.values)
        unique_true = unique[unique > 0]
        cluster_means = []
        for i in unique_true:
            cluster_mean = np.mean(group_data[group_data.cell_id == i][['x', 'y']].values, axis=0)
            cluster_means.append(cluster_mean)
        cluster_means = np.array(cluster_means)
        #Find which the index of the point in cluster means that is closest to each point in group_data[['x', 'y']].values
        # Calculate the Euclidean distance between each point in group_data and each cluster mean
        distances = np.linalg.norm(cluster_means[:, np.newaxis, :] - group_data[['x', 'y']].values, axis=2)
        
        # Find the index of the cluster mean with the minimum distance for each point
        closest_indices = unique_true[np.argmin(distances, axis=0)]
        reassignments[group_data.indexer.values] = closest_indices
        
    transcripts_with_gt_and_main_nucleus_filtered['split_cell'] = reassignments
    return transcripts_with_gt_and_main_nucleus_filtered

def make_adata(transcripts_with_gt_and_main_nucleus_filtered):
    cell_by_gene = transcripts_with_gt_and_main_nucleus_filtered.groupby(['split_cell', 'gene']).size().unstack(fill_value=0)
    transcripts_nucleus = transcripts_with_gt_and_main_nucleus_filtered[transcripts_with_gt_and_main_nucleus_filtered['overlaps_nucleus'] == 1]
    cell_by_gene_nucleus = transcripts_nucleus.groupby(['split_cell', 'gene']).size().unstack(fill_value=0)
    cell_by_gene = cell_by_gene.loc[cell_by_gene_nucleus.index.values]
    cell_by_gene_nucleus = cell_by_gene_nucleus.loc[cell_by_gene.index.values]
    
    new_cyto_nuc = pd.DataFrame(zip(np.sum(cell_by_gene, axis=1), np.sum(cell_by_gene_nucleus, axis=1)), index= cell_by_gene.index.astype(int), columns=['total_transcripts', 'nuclear_transcripts'])    
    anndata = sc.AnnData(cell_by_gene.values, var=pd.DataFrame(index=cell_by_gene.columns), obs=new_cyto_nuc)
    anndata.layers['raw'] = anndata.X
    anndata.obs['cytoplasmic_transcripts'] = anndata.obs['total_transcripts'] - anndata.obs['nuclear_transcripts'] 
    anndata.obs['nuclear_transcript_percentage'] = anndata.obs['nuclear_transcripts']/anndata.obs['total_transcripts']
    anndata.var['gene'] = anndata.var.index.values
    anndata.obs['cell'] = anndata.obs.index.values
    cell_spatial = transcripts_with_gt_and_main_nucleus_filtered.groupby('split_cell')[['x', 'y']].mean()
    cell_spatial =cell_spatial.set_index(pd.Index(cell_spatial.index.values.astype(int).astype(str)))
    anndata.uns['points'] = transcripts_with_gt_and_main_nucleus_filtered
    anndata.obs = anndata.obs.set_index(pd.Index(anndata.obs.index.values.astype(int).astype(str)))
    anndata.obs = anndata.obs.merge(cell_spatial, how='left', left_index=True, right_index=True)
    anndata.obsm['X_spatial'] = anndata.obs[['x', 'y']].values
    anndata = anndata[:, ~((anndata.var.index.str.contains('BLANK')) | (anndata.var.index.str.contains('NegControl')))]
    return anndata