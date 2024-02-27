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
import json


def make_alphashape(points: pd.DataFrame, alpha: float):
    """
    Create an alpha shape from a set of points

    Parameters:
    points (pd.DataFrame): A dataframe of points
    alpha (float): The alpha parameter to use for the alpha shape

    Returns:
    shape (shapely.geometry.Polygon): The alpha shape
    """
    points = np.array(points)
    shape = alphashape.alphashape(points, alpha=alpha)
    return shape


def get_pixel_size(path: str) -> float:
    """
    Get the pixel size from the experiment file

    Parameters:
    path (str): The path to the experiment file

    Returns:
    pixel_size (float): The pixel size
    """
    file = open(os.path.join(path, "experiment.xenium"))
    experiment = json.load(file)
    pixel_size = experiment["pixel_size"]
    return pixel_size


def subset_transcripts_file(transcripts, pixel_size, minx, maxx, miny, maxy):
    """
    Subset the transcripts file to only include those within the field of view

    Parameters:
    transcripts (pd.DataFrame): A DataFrame containing the transcripts
    pixel_size (float): The pixel size
    minx (float): The minimum x coordinate of the field of view
    maxx (float): The maximum x coordinate of the field of view
    miny (float): The minimum y coordinate of the field of view
    maxy (float): The maximum y coordinate of the field of view

    Returns:
    transcript_subset_fov (pd.DataFrame): The subset of the transcripts within the field of view
    """
    transcript_subset_fov = transcripts[
        (minx < transcripts.y * (1 / pixel_size))
        & (transcripts.y * (1 / pixel_size) < maxx)
        & (miny < transcripts.x * (1 / pixel_size))
        & (transcripts.x * (1 / pixel_size) < maxy)
    ]
    return transcript_subset_fov


def import_image(path: str):
    """
    Import the max-projected DAPI stain from the Xenium

    Parameters:
    path (str): The path to the experiment directory

    Returns:
    img (np.ndarray): The max-projected DAPI stain
    """
    file = os.path.join(path, "morphology_mip.ome.tif")
    img = io.imread(file)
    return img
