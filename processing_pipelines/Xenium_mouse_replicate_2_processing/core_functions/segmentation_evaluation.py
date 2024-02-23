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
from fast_alphashape import alphashape as alphafast
import imageio as io
import json

def make_alphashape(points: pd.DataFrame, alpha: float):
    points = np.array(points)
    shape = alphashape.alphashape(points, alpha=alpha)
    return shape

def make_alphashape_fast(points: pd.DataFrame, alpha: float):
    points = np.array(points)
    shape = alphafast(points, alpha=alpha)
    return shape

def get_pixel_size(path: str) -> float:
    file = open(os.path.join(path, "experiment.xenium"))
    experiment = json.load(file)
    pixel_size = experiment['pixel_size']
    return pixel_size

def subset_transcripts_file(transcripts, pixel_size, minx, maxx, miny, maxy):
    transcript_subset_fov = transcripts[(minx < transcripts.y*(1/pixel_size)) & (transcripts.y*(1/pixel_size) < maxx)& (miny < transcripts.x*(1/pixel_size)) & (transcripts.x*(1/pixel_size) < maxy)]
    return transcript_subset_fov

def import_image(path: str):
    file = os.path.join(path, "morphology_mip.ome.tif")
    img = io.imread(file)
    return img