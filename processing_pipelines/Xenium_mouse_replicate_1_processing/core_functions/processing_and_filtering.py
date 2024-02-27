import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from anndata import AnnData
import os


def plot_qc_feature(cell_by_gene, cell_meta, path, control_probes=False):
    """
    Plot quality control features

    Parameters:
    cell_by_gene (pd.DataFrame): A DataFrame containing the cell by gene matrix
    cell_meta (pd.DataFrame): A DataFrame containing the cell metadata
    path (str): The path to save the plot
    control_probes (bool): Whether or not to plot control probes

    Returns:
    None
    """

    qc_info = cell_meta

    if control_probes:
        colors = ["#98FB98", "#E6E6FA", "#ADD8E6", "#FFB6C1", "#FFE4E1"]
        metrics = [
            "total_transcripts",
            "nuclear_transcripts",
            "cytoplasmic_transcripts",
            "nuclear_transcript_percentage",
            "control_probe_counts",
        ]
        ncols = 5
    else:
        colors = ["#98FB98", "#E6E6FA", "#ADD8E6", "#FFB6C1"]
        metrics = [
            "total_transcripts",
            "nuclear_transcripts",
            "cytoplasmic_transcripts",
            "nuclear_transcript_percentage",
        ]
        ncols = 4

    fig, axes = plt.subplots(
        figsize=(15, 3), dpi=200, ncols=ncols, constrained_layout=True
    )
    # Loop through each axis and create violin and strip plots
    for i, ax in enumerate(axes):
        sns.violinplot(y=qc_info[metrics[i]], color=colors[i], ax=ax)
        sns.stripplot(
            y=qc_info[metrics[i]], jitter=True, color="black", ax=ax, size=0.1
        )
        ax.set_title(metrics[i])  # Set title for each subplot
    try:
        os.mkdir(os.path.join(path, "figures", "quality"))
    except:
        print("quality directory already made")
    # Show the plots
    fig.tight_layout()
    fig.savefig(os.path.join(path, "figures", "quality", "qc.png"))
    plt.show()


def qc_before_clustering(
    adata,
    min_transcript_threshold=20,
    max_transcript_threshold=800,
    min_nuclear_transcripts=8,
    max_nuclear_transcripts=600,
    min_cyto_transcripts=-1,
    max_cyto_transcripts=700,
    min_nuc_pct=0,
    max_nuc_pct=1.01,
):
    """
    Perform quality control filtering

    Parameters:
    adata (AnnData): An AnnData object containing the original data

    Returns:
    adata (AnnData): An AnnData object containing the filtered data
    """
    print(f"{len(adata.obs.index)} cells before QC filtering")
    adata = adata[
        (adata.obs["total_transcripts"] > min_transcript_threshold)
        & (adata.obs["total_transcripts"] < max_transcript_threshold),
        :,
    ]
    adata = adata[
        (adata.obs["nuclear_transcripts"] > min_nuclear_transcripts)
        & (adata.obs["nuclear_transcripts"] < max_nuclear_transcripts),
        :,
    ]
    adata = adata[
        (adata.obs["cytoplasmic_transcripts"] > min_cyto_transcripts)
        & (adata.obs["cytoplasmic_transcripts"] < max_cyto_transcripts),
        :,
    ]
    adata = adata[
        (adata.obs["nuclear_transcript_percentage"] > min_nuc_pct)
        & (adata.obs["nuclear_transcript_percentage"] < max_nuc_pct),
        :,
    ]
    print(f"{len(adata.obs.index)} cells after QC filtering")

    return adata
