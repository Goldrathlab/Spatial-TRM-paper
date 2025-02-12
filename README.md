# Tissue-resident memory CD8 T cell diversity is spatiotemporally imprinted
![Project Banner](images/connections.png)

> [!NOTE]
> In this repository, we show how to reproduce the figures from our 2025 manuscript: **Tissue-resident memory CD8 T cell diversity is spatiotemporally imprinted**[^1][^2]. Additionally, we provide our data processing pipelines to create fully-processed Anndata objects containing all of the spatial data used to construct figures.

## Table of Contents

- [Abstract](#abstract)
- [Setup](#setup)
- [Download](#download)
- [Preprocessing](#preprocessing)
- [Figures](#figures)
- [Submitting changes](#submitting-changes)

## Abstract

Tissue-resident memory CD8 T (T<sub>RM</sub>) cells provide protection from infection at barrier sites. In the small intestine, T<sub>RM</sub> cells are found in at least two distinct subpopulations: one with higher expression of effector molecules and another with greater memory potential1. However, the origins of this diversity remain unknown. Here we proposed that distinct tissue niches drive the phenotypic heterogeneity of T<sub>RM</sub> cells. To test this, we leveraged spatial transcriptomics of human samples, a mouse model of acute systemic viral infection and a newly established strategy for pooled optically encoded gene perturbations to profile the locations, interactions and transcriptomes of pathogen-specific T<sub>RM</sub> cell differentiation at single-transcript resolution. We developed computational approaches to capture cellular locations along three anatomical axes of the small intestine and to visualize the spatiotemporal distribution of cell types and gene expression. Our study reveals that the regionalized signalling of the intestinal architecture supports two distinct T<sub>RM</sub> cell states: differentiated T<sub>RM</sub> cells and progenitor-like T<sub>RM</sub> cells, located in the upper villus and lower villus, respectively. This diversity is mediated by distinct ligand–receptor activities, cytokine gradients and specialized cellular contacts. Blocking TGFβ or CXCL9 and CXCL10 sensing by antigen-specific CD8 T cells revealed a model consistent with anatomically delineated, early fate specification. Ultimately, our framework for the study of tissue immune networks reveals that T cell location and functional state are fundamentally intertwined.

## Setup

This repository provides a `devcontainer` setup to ensure that scripts run consistently and reproducibly. For detailed instructions on using devcontainers, please refer to the [official documentation](https://code.visualstudio.com/docs/devcontainers/containers).

The `devcontainer` is built using a minimal Docker container and includes the [pixi](https://prefix.dev/) package manager, which automatically installs the dependencies specified in the `pixi.toml` file.

If you prefer, you can manually install the dependencies using `pixi install`.

Please note that the analysis was conducted on a Linux machine running Ubuntu 22.04 LTS. Many dependencies require a x64 architecture and may not be available for macOS with ARM64 architecture.


## Download

To reproduce the main figures, you will need additional files located in the `data` directory. A script to download the processed (and raw) data is provided [here](/data/download.ipynb).

List of required files:

```text
data
├── adata
│   ├── day8_r2_with_transcripts.h5ad
│   ├── human_09_r2_with_transcripts.h5ad
│   ├── human.h5ad
│   ├── perturb.h5ad
│   ├── tgfb.h5ad
│   ├── timecourse.h5ad
│   ├── uninfected.h5ad
│   └── visium_hd.h5ad
├── IF
│   └── timecourse
│       ├── day 005.txt
│       ├── day 006.txt
│       ├── day 007.txt
│       ├── day 060.txt
│       └── day 120.txt
├── images
│   ├── day8_r2_h_and_e_alignment_gan.npy
│   ├── day8_r2_IF_alignment.npy
│   ├── human_09_r2_h_and_e_alignment_gan.npy
│   └── human_09_r2_IF_alignment.npy
├── kegg_cytokines.csv
├── signatures
│   ├── blimp.txt
│   ├── Core Trm signature_Milner et al Nature 2017_vIL.txt
│   ├── id3.txt
│   ├── kurd.xlsx
│   └── TGFbeta.txt
├── transcripts
│   └── transcripts_figure_5c.csv
└── xenium_output
    ├── day8_r2
    │   ├── experiment.xenium
    │   └── morphology_mip.ome.tif
    └── human_09_r2
        ├── experiment.xenium
        └── morphology_mip.ome.tif
```

The complete data repository that accompanies the manuscript submission can be accessed [here](https://2024-spatial-trm.data.heeg.io/).


## Preprocessing

Data from 10x Xenium and Vizgen MERSCOPE were preprocessed using a custom segmentation and annotation pipeline. These pipelines can be found on seperate GitHub repositories

- [Xenium nextflow segmentation pipeline](https://github.com/maximilian-heeg/xenium-segmentation) v0.1.2.
  The following parameters were set:

  ```text
  tile.minimal_transcripts = 300000
  baysor.prior_segmentation_confidence = 0.95
  ```

- [Merscope nextflow segmentation pipeline](https://github.com/maximilian-heeg/vizgen-segmentation/) v0.1.0

  ```text
  tile.minimal_transcripts = 5000000
  baysor.prior_segmentation_confidence = 0.9
  ```

In this project, we had similar but separate processing pipelines for processing data from different settings.

1. Processing of Xenium mouse small intestine timecourse first replicates.
   [Xenium mouse rep 1 processing](/processing_pipelines/Xenium_mouse_replicate_1_processing)
2. Processing of Xenium mouse small intestine timecourse second replicates.
   [Xenium mouse rep 2 processing](/processing_pipelines/Xenium_mouse_replicate_2_processing)
3. Processing of MERSCOPE mouse small intestine WT vs TGFBR2 KO conditions.
   [MERSCOPE mouse processing](/processing_pipelines/MERSCOPE_mouse_processing)
4. Processing of Xenium human terminal ileum both replicates.
   [Xenium human processing](/processing_pipelines/Xenium_human_processing)
5. Processing of VisiumHD mouse small intestine.
   [VisiumHD processing](/processing_pipelines/visiumHD)
6. Processing of Uninfected mouse small intestine.
   [Uninfected Processing](/processing_pipelines/Uninfected_Processing)
7. Processing of pooled CRISPR KO spatial.
   [Spatial Perturb Processing](/processing_pipelines/perturb_spatial)

> [!NOTE]
> Part of our workflow included cell type annotation. This made use of manual exploration and is not well reflected in the code. We have provided excel sheets and csvs used to assign cell type annotations to cell clusters in each pipeline folder. Additionally, we labeled images manually in several parts of our pipelines. We have provided these labels in json format within each pipeline folder. Please [contact](#contact) us if you need us to give you intermediary objects at any point in the processing pipelines.

We also performed Immunofluorescence staining and H&E after our Xenium runs.
We show pipelines for aligning IF and H&E images with our Xenium data in:

1. Mouse data [Mouse histology alignment](/processing_pipelines/alignment/mouse_histology)
2. Human data [Human histology alignment](/processing_pipelines/alignment/human_histology)

## Figures

This section contains the scripts to reproduce the figures in the paper.

### Figure 1

| Figure | Link                                                  |
| ------ | ----------------------------------------------------- |
| 1b     | [Immunofluorescence IMAP](/Figure_1/1b.ipynb)         |
| 1c     | [Xenium, H&E, IF and transcripts](/Figure_1/1c.ipynb) |
| 1d     | [Overview plots](/Figure_1/1d.ipynb)                  |

### Figure 2

| Figure                 | Link                                                |
| ---------------------- | --------------------------------------------------- |
| 2a <br /> 2b <br /> 2c | [Axes in the small intestine](/Figure_2/2abc.ipynb) |
| 2d                     | [IMAP of P14](/Figure_2/2d.ipynb)                   |
| 2e                     | [Correlation along the axes](/Figure_2/2e.ipynb)    |
| 2f <br /> 2i           | [Convolved expression heatmap](/Figure_2/2fi.ipynb)  |
| 2g <br /> 2h           | [Expression IMAP of P14](/Figure_2/2gh.ipynb)       |
| 2j                     | [Signature IMAP of P14](/Figure_2/2j.ipynb)         |

### Figure 3

| Figure       | Link                                                                                                                                                                                  |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 3a           | [Spatial graph of the small intestine](/Figure_3/3a.ipynb)                                                                                                                            |
| 3b <br /> 3c | [Co-localization of cells and P14 subsets](/Figure_3/3bc.ipynb)                                                                                                                       |
| 3d           | [Convolved expression heatmap of cytokines](/Figure_3/3d.ipynb)                                                                                                                       |
| 3e           | [TGFb isoforms in the small intestine](/Figure_3/3e_5a.ipynb)    <br />   [Xenium picture](/Figure_3/3e_picture.ipynb)                                                                |
| 3f           | Cellchat <br /> [Part 1: Preprocessing python](/Figure_3/3g_part1.ipynb) <br /> [Part 2: Run Cellchat](/Figure_3/3g_part2.ipynb) <br /> [Part 3: Make plot](/Figure_3/3g_part3.ipynb) |

### Figure 4

| Figure                            | Link                                                           |
| --------------------------------- | -------------------------------------------------------------- |
| Preprocessing for Figure 4 object | [Figure 4 Preprocessing](/Figure_4/figure_object_creation)     |
| 4a                                | [Cxcl9 and Cxcl10 expression](/Figure_4/4a.ipynb)              |
| 4d                                | [Gene depletion](/Figure_4/4d.ipynb)                           |
| 4e <br /> 4f                      | [IMAPs and gates](/Figure_4/4ef_Extended_8e.ipynb)             |

### Figure 5

| Figure       | Link                                                                                                                                                                                                                                                                                                                                                 |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 5a <br /> 5b | [MDE plot and cell frequencies](/Figure_5/5ab.ipynb)                                                                                                                                                                                                                                                                                                 |
| 5c           | [Xenium, H&E, IF and transcripts](/Figure_5/5c.ipynb)                                                                                                                                                                                                                                                                                                |
| 5d           | [Mouse signature on human IMAP](/Figure_5/5d.ipynb)                                                                                                                                                                                                                                                                                                  |
| 5e <br /> 5f | [Convolved expression heatmap](/Figure_5/5ef.ipynb)                                                                                                                                                                                                                                                                                                  |
| 5g           | [Correlation between gene expression and cell types](/Figure_5/5g.ipynb)                                                                                                                                                                                                                                                                             |
| 5h           | Cellchat <br /> [Part 1: Dataset export](/Figure_5/5h_part1_dataset_export_signature.ipynb) <br /> [Part 2: Cellchat preparation](/Figure_5/5h_part2_cellchat_preparation%20signature.ipynb) <br /> [Part 3: Run Cellchat](/Figure_5/5h_part3_spatial_cellchat_signature.ipynb) <br /> [Part 4: Make plot](/Figure_5/5h_part4_output_cellchat.ipynb) |

### Extended Data Figure 7

| Figure       | Link                                                          |
| ------------ | ------------------------------------------------------------- |
| ED 7a           | [Overview plots](/Extended_Data_Figure_7/ED7a.ipynb)                          |
| ED 7b <br /> ED 7e | [IMAP and expression IMAP](/Extended_Data_Figure_7/ED7be.ipynb)               |
| ED 7d           | [Differenitally expressed genes](/Extended_Data_Figure_7/ED7d.ipynb)          |
| ED 7g <br /> ED 7h | [Correlation of expression to celltypes](Extended_Data_Figure_7/ED7gh.ipynb) |

### Other extended data figures

| Figure             | Link                                       |
| ------------------ | ------------------------------------------ |
| ED 6a              | [Notebook](/Figure_3/3g_part3.ipynb)       |
| ED 6i              | [Notebook](/Extended_Data_Figure_7/ED7be.ipynb)            |
| ED 9b <br /> ED 9d | [Notebook](/Figure_5/5ab.ipynb)            |
| ED 8e              | [Notebook](/Figure_4/4ef_Extended_8e.ipynb)|
| ED 8b <br /> ED 8d | [Notebook](/Figure_4/Extended_4bd.ipynb)   |
| ED 8c              | [Notebook](/Figure_4/Extended_4c.ipynb)    |

## Submitting changes

> [!IMPORTANT]
> Make sure to run `pixi run pre-commit` before commiting

## Contact

- Miguel Reina-Campos: :envelope: miguel@lji.org
- Alexander Monell: :envelope: amonell@ucsd.edu
- Maximilian Heeg: :envelope: maximilian.heeg@alleninstitute.org
- Ananda Goldrath: :envelope: agoldrath@alleninstitute.org

[^1]: [Preprint at bioRxiv](https://www.biorxiv.org/content/10.1101/2024.03.20.585130v1)
[^2]: [Publication at _Nature_](https://www.nature.com/articles/s41586-024-08466-x)
