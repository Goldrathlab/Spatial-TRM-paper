import tensorflow as tf
import scanpy as sc
import os
from scipy.spatial import KDTree
from tensorflow import keras
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from tqdm.notebook import tqdm
import glob


def create_grid_bins(spatial_points, n):
    """
    Create a grid of bins to assign spatial points to

    Parameters:
    spatial_points (np.array): An array of spatial points
    n (int): The number of bins to create

    Returns:
    grid_bins (np.array): An array of bins
    bin_centers (np.array): An array of bin centers
    """
    xmin, ymin = np.min(spatial_points, axis=0)
    xmax, ymax = np.max(spatial_points, axis=0)

    xbins = np.linspace(xmin, xmax, n + 1)
    ybins = np.linspace(ymin, ymax, n + 1)

    grid_bins = [[[] for _ in range(n)] for _ in range(n)]
    bin_centers = []

    for i in range(n):
        for j in range(n):
            bin_center_x = (xbins[i] + xbins[i + 1]) / 2
            bin_center_y = (ybins[j] + ybins[j + 1]) / 2
            bin_centers.append([bin_center_x, bin_center_y])

    for point in tqdm(range(len(spatial_points))):

        x, y = spatial_points[point]
        xi = np.searchsorted(xbins, x, side="right") - 1
        yi = np.searchsorted(ybins, y, side="right") - 1

        try:
            grid_bins[xi][yi].append(point)
        except:
            None

    return grid_bins, bin_centers


def create_binned_data(adata, bins, centers, unique_bins):
    """
    Create a binned AnnData object from the original AnnData object

    Parameters:
    adata (AnnData): An AnnData object containing the original data
    bins (np.array): An array of bins
    centers (np.array): An array of bin centers
    unique_bins (np.array): An array of unique bins

    Returns:
    adata_filtered (AnnData): An AnnData object containing the binned data
    """
    expression_matrix = []
    for b in tqdm(range(len(unique_bins))):
        where_bin = np.where(bins == b)[0]
        try:
            bin_expression = np.array(
                np.mean(adata.X[where_bin, :], axis=0).flatten()
            ).squeeze()
        except:
            bin_expression = np.array([float(0) for i in range(len(adata.var.index))])
        expression_matrix.append(bin_expression)
    expression_matrix = np.array(expression_matrix)
    ad = sc.AnnData(
        X=expression_matrix,
        obs=pd.DataFrame(index=unique_bins),
        var=pd.DataFrame(index=adata.var.index.tolist()),
    )
    ad.obsm["spatial"] = np.array(centers)
    nan_obs_indices = np.where(np.isnan(ad.X.sum(axis=1)))[0]

    # Filter out observations with NaN values
    adata_filtered = ad[~np.isin(ad.obs_names, ad.obs_names[nan_obs_indices])].copy()
    return adata_filtered


def plot_top_words(model, feature_names, n_top_words, title):
    """
    Plot the top words for each topic in the NMF model.

    Parameters:
    model (NMF): An NMF model
    feature_names (list): A list of feature names
    n_top_words (int): The number of top words to plot
    title (str): The title of the plot

    Returns:
    None
    """
    fig, axes = plt.subplots(
        3, int(len(model.components_) / 3) + 1, figsize=(15, 7), sharex=False
    )
    axes = axes.flatten()
    for topic_idx, topic in enumerate(model.components_):
        top_features_ind = topic.argsort()[: -n_top_words - 1 : -1]
        top_features = [feature_names[i] for i in top_features_ind]
        weights = topic[top_features_ind] / np.sum(topic[top_features_ind])

        ax = axes[topic_idx]
        ax.barh(top_features, weights, height=0.7, color="r")

        ax.set_xlim(0, np.max(weights))
        ax.set_title(f"Neighborhood {topic_idx +1}", fontdict={"fontsize": 15})
        ax.invert_yaxis()
        ax.tick_params(axis="both", which="major", labelsize=10)
        for i in "top right left".split():
            ax.spines[i].set_visible(False)
        fig.suptitle(title, fontsize=20)

    plt.subplots_adjust(top=0.90, bottom=0.05, wspace=0.90, hspace=0.3)
    plt.show()
    # fig.savefig(r'D:\Alex\MERSCOPE_reanalysis_output\24hr_kt56.pdf')
    plt.close()
