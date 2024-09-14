# Functional Diversity of Memory CD8 T Cells is Spatiotemporally Imprinted
![Project Banner](images/connections.png)

> [!NOTE]
> In this repository, we show how to reproduce the figures from our 2024 manuscript: **Functional Diversity of Memory CD8 T Cells is Spatiotemporally Imprinted**[^1]. Additionally, we provide our data processing pipelines to create fully-processed Anndata objects containing all of the spatial data used to construct figures.

## Table of Contents

- [Abstract](#abstract)
- [Setup](#setup)
- [Download](#download)
- [Preprocessing](#preprocessing)
- [Figures](#figures)
- [Submitting changes](#submittingchanges)

## Abstract

Tissue-resident memory CD8 T cells (T<sub>RM</sub>) kill infected cells and recruit additional immune cells to limit pathogen invasion at barrier sites. Small intestinal (SI) T<sub>RM</sub> cells consist of distinct subpopulations with higher expression of effector molecules or greater memory potential. We hypothesized that occupancy of diverse anatomical niches imprints these distinct T<sub>RM</sub> transcriptional programs. We leveraged human samples and a murine model of acute systemic viral infection to profile the location and transcriptome of pathogen-specific T<sub>RM</sub> cell differentiation at single-transcript resolution. We developed computational approaches to capture cellular locations along three anatomical axes of the small intestine and to visualize the spatiotemporal distribution of cell types and gene expression. T<sub>RM</sub> populations were spatially segregated: with more effector- and memory-like T<sub>RM</sub> preferentially localized at the villus tip or crypt, respectively. Modeling ligand-receptor activity revealed patterns of key cellular interactions and cytokine signaling pathways that initiate and maintain T<sub>RM</sub> differentiation and functional diversity, including different TGFβ sources. Alterations in the cellular networks induced by loss of TGFβRII expression revealed a model consistent with TGFβ promoting progressive T<sub>RM</sub> maturation towards the villus tip. Ultimately, we have developed a framework for the study of immune cell interactions with the spectrum of tissue cell types, revealing that T cell location and functional state are fundamentally intertwined.

## Setup

This repository contains a `devcontainer` to allow to run the scripts in a reproducible manner. Please see the [documentation](https://code.visualstudio.com/docs/devcontainers/containers) for further information on how to use devcontainers.

## Download

The preprocessed files will be provided on Zenodo upon publication.
To reproduce the main figures, the following files in the `data` directory are needed.

> [!IMPORTANT]
> A script to download the processed data is included [here](/data/download.ipynb). In order to run this script, a base url is needed. This url will be provided with the submission of the manuscript.

```text
data
├── adata
│   ├── day8_r2_with_transcripts.h5ad
│   ├── human_09_r2_with_transcripts.h5ad
│   ├── human.h5ad
│   ├── perturb.h5ad
│   ├── tgfb.h5ad
│   ├── timecourse.h5ad
│   └── uninfected.h5ad
├── download.ipynb
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

In this project, we had four similar but separate processing pipelines for processing data from different settings.

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

| Figure       | Link                                                          |
| ------------ | ------------------------------------------------------------- |
| 4a           | [Overview plots](/Figure_4/4a.ipynb)                          |
| 4b <br /> 4e | [IMAP and expression IMAP](/Figure_4/4be.ipynb)               |
| 4d           | [Differenitally expressed genes](/Figure_4/4d.ipynb)          |
| 4g <br /> 4h | [Correlation of expression to celltypes](/Figure_4/4gh.ipynb) |

### Figure 5

| Figure                            | Link                                                           |
| --------------------------------- | -------------------------------------------------------------- |
| Preprocessing for Figure 5 object | [Figure 5 Preprocessing](/Figure_5/figure_object_creation)     |
| 5a                                | [Cxcl9 and Cxcl10 expression](/Figure_5/5a.ipynb)              |
| 5d                                | [Gene depletion](/Figure_5/5d.ipynb)                           |
| 5e <br /> 5f                      | [IMAPs and gates](/Figure_5/5ef_Extended_7e.ipynb)             |

### Figure 6

| Figure       | Link                                                                                                                                                                                                                                                                                                                                                 |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 6a <br /> 6b | [MDE plot and cell frequencies](/Figure_6/6ab.ipynb)                                                                                                                                                                                                                                                                                                 |
| 6c           | [Xenium, H&E, IF and transcripts](/Figure_6/6c.ipynb)                                                                                                                                                                                                                                                                                                |
| 6d           | [Mouse signature on human IMAP](/Figure_6/6d.ipynb)                                                                                                                                                                                                                                                                                                  |
| 6e <br /> 6f | [Convolved expression heatmap](/Figure_6/6ef.ipynb)                                                                                                                                                                                                                                                                                                  |
| 6g           | [Correlation between gene expression and cell types](/Figure_6/6g.ipynb)                                                                                                                                                                                                                                                                             |
| 6h           | Cellchat <br /> [Part 1: Dataset export](/Figure_6/6h_part1_dataset_export_signature.ipynb) <br /> [Part 2: Cellchat preparation](/Figure_6/6h_part2_cellchat_preparation%20signature.ipynb) <br /> [Part 3: Run Cellchat](/Figure_6/6h_part3_spatial_cellchat_signature.ipynb) <br /> [Part 4: Make plot](/Figure_6/6h_part4_output_cellchat.ipynb) |

### Extended data figures

| Figure             | Link                                       |
| ------------------ | ------------------------------------------ |
| ED 6a              | [Notebook](/Figure_3/3g_part3.ipynb)       |
| ED 6f              | [Notebook](/Figure_4/4be.ipynb)            |
| ED 8b <br /> ED 8d | [Notebook](/Figure_6/6ab.ipynb)            |
| ED 7e              | [Notebook](/Figure_5/5ef_Extended_7e.ipynb)|
| ED 7b <br /> ED 7d | [Notebook](/Figure_5/Extended_7bd.ipynb)   |
| ED 7c              | [Notebook](/Figure_5/Extended_7c.ipynb)    |

## Submitting changes

> [!IMPORTANT]
> Make sure to run `pixi run pre-commit` before commiting

## Contact

- Miguel Reina-Campos: :envelope: miguel@lji.org
- Alexander Monell: :envelope: amonell@ucsd.edu
- Maximilian Heeg: :envelope: maximilian.heeg@alleninstitute.org
- Ananda Goldrath: :envelope: agoldrath@alleninstitute.org

[^1]: [Preprint at bioRxiv](https://www.biorxiv.org/content/10.1101/2024.03.20.585130v1)
