# Paper title here
![Project Banner](path/to/banner_image.png)

> :information_source:  In this repository, we show how to reproduce the figures from our 2024 manuscript (link to paper). Additionally, we provide our data processing pipelines to create fully-processed Anndata objects containing all of the spatial data used to construct figures.[^1].

## Table of Contents

- [Abstract](#abstract)
- [Setup](#setup)
- [Download](#download)
- [Preprocessing](#preprocessing)
- [Figures](#figures)
- [Submitting changes](#submittingchanges)

## Abstract

## Setup

This repository contains a `devcontainer` to allow to run the scripts in a reproducible manner. Please see the [documentation](https://code.visualstudio.com/docs/devcontainers/containers) for further information on how to use devcontainers.

## Download

TODO: Add Script to download the data

## Preprocessing

Data from 10x Xenium and Vizgen MERSCOPE were preprocessed using a custom segmentation and annotation pipeline. These pipelines can be found on seperate GitHub repositories

- [Xenium nextflow segmentation pipeline](https://github.com/maximilian-heeg/xenium-segmentation) v0.1.2
  The following parameters were set:

  ```text
  tile.minimal_transcripts = 300000
  baysor.prior_segmentation_confidence = 0.95
  ```

- [Merscope nextflow segmentation pipeline](https://github.com/maximilian-heeg/vizgen-segmentation/) v0.1.0

  ```text
  tile.minimal_transcripts = 500000
  baysor.prior_segmentation_confidence = 0.9
  ```

In this project, we had four separate downstream processing pipelines for processing data from different settings.
1. Processing of Xenium mouse small intestine timecourse first replicates.
   [Pipeline](/processing_pipelines/Xenium_mouse_replicate_1_processing)
3. Processing of Xenium mouse small intestine timecourse second replicates.
4. Processing of MERSCOPE mouse small intestine WT vs TGFBR2 KO conditions.
5. Processing of Xenium human terminal ileum both replicates.
## Figures

This section contains the scripts to reproduce the figures in the paper.

### Figure 1

| Figure | Link                           |
|--------|--------------------------------|
| 1b     | TODO                           |
| 1c     | TODO                           |
| 1d     | [Notebook](/Figure_1/1d.ipynb) |

### Figure 2

| Figure                | Link                             |
|-----------------------|----------------------------------|
| 2a <br /> 2b<br /> 2c | [Notebook](/Figure_2/2abc.ipynb) |
| 2d                    | [Notebook](/Figure_2/2d.ipynb)   |
| 2e <br /> 2f          | [Notebook](/Figure_2/2ef.ipynb)  |
| 2g                    | [Notebook](/Figure_2/2g.ipynb)   |
| 2h <br /> 2i          | [Notebook](/Figure_2/2hi.ipynb)  |
| 2j                    | [Notebook](/Figure_2/2j.ipynb)   |

### Extended data figures

| Figure | Link                           |
|--------|--------------------------------|
| ED 2b  | [Notebook](/Figure_2/2g.ipynb) |
| ED 2d  | [Notebook](/Figure_2/2j.ipynb) |

## Submitting changes

Make sure to run `pre-commit run --all-files` before commiting

[^1]: TODO add link to paper
