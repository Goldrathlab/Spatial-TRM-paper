import scanpy as sc
import numpy as np
from tqdm.notebook import tqdm
import scipy.stats as stats
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import shapely
import glob
from sklearn.neighbors import NearestNeighbors
from PIL import Image, ImageDraw
import json
from scipy.spatial import cKDTree


def remove_outliers(spatial_points, thresh=99):
    """
    Remove outliers from the set of spatial points in the spiral.

    Parameters:
    spatial_points (np.ndarray): An array of spatial points
    thresh (int): The percentile threshold to use for removing outliers

    Returns:
    np.ndarray: The spatial points with outliers removed
    """
    ### Removing outliers
    # Step 1: Compute distances between each point
    nn = 100
    nbrs = NearestNeighbors(n_neighbors=nn, algorithm="kd_tree").fit(spatial_points)
    distances, _ = nbrs.kneighbors(spatial_points)

    # Step 2: For each point, remove the distance to itself (which will be 0)
    distances = distances[:, 1:]

    # Step 3: Calculate the average distance to the 5 nearest neighbors for each point
    avg_distances = np.mean(distances, axis=1)

    cutoff = np.percentile(avg_distances, 99)
    spatial_points = spatial_points[avg_distances < cutoff]
    return spatial_points


def create_base_image(spatial_points, other_spatial, downsize=20):
    """
    Draw and save image for labeling

    Parameters:
    spatial_points (np.ndarray): The array of spatial points from the experiment to plot
    other_spatial (np.ndarray): An array of spatial points not belonging to neighborhoods in the basal membrane.
    downsize (int): The amount to downsize the image

    Returns:
    Image: The base image
    """
    # Draw and save image

    points = spatial_points / downsize
    other_points = other_spatial / downsize

    # Define the size of the image (adjust as needed)
    image_width = 2000
    image_height = 2000

    # Create a white canvas as the base image
    base_image = Image.new("RGB", (image_width, image_height), (255, 255, 255))

    # Draw the points on the image
    draw = ImageDraw.Draw(base_image)
    point_size = 1  # Size of the points
    ct = 0
    for point in points:
        draw.ellipse(
            (
                point[0] - point_size,
                point[1] - point_size,
                point[0] + point_size,
                point[1] + point_size,
            ),
            fill="blue",
        )
        ct += 1

    ct = 0
    for point in other_points:
        draw.ellipse(
            (
                point[0] - point_size,
                point[1] - point_size,
                point[0] + point_size,
                point[1] + point_size,
            ),
            fill="red",
        )
        ct += 1

    return base_image


def extract_json_info(json_file_path):
    """
    Extract relevant information from the JSON file of labels.

    Parameters:
    json_file_path (str): The path to the JSON file with the image labels

    Returns:
    removals: A list of the points to remove
    points: A label going along the base of the basal membrane.
    top_points: A label going along the tips of villi.
    mid_points: A label going through the middle of the villi.
    """
    # Load the JSON data from the file
    with open(json_file_path, "r") as json_file:
        data = json.load(json_file)

    # Extract relevant information from the JSON data
    image_height = data["imageHeight"]
    image_width = data["imageWidth"]
    image_path = data["imagePath"]
    shapes = data["shapes"]

    # Process the shapes (annotations)
    removals = []
    points = []
    top_points = []
    mid_points = []
    for shape in shapes:
        label = shape["label"]
        if (label == "roll") or (label == "point"):
            points.append(shape["points"])
        elif label == "top_point":
            top_points.append(shape["points"])
        elif label == "mid":
            mid_points.append(shape["points"])
        else:
            removals.append(shape["points"])
    return removals, points, top_points, mid_points


def get_removal_indices(adata, removals, all_spatial):
    """
    Get the indices of the points to remove from the spatial data.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the spatial data.
    removals (list): A list of the points to remove.
    all_spatial (np.ndarray): The total list of spatial points.

    Returns:
    dont_remove: The indices of the points to keep
    set: A set of the indices of the points to keep
    """
    total_indices = []
    for ir in removals:
        ir_ = np.array(ir) * adata.uns["unrolling_downsize"]
        poly = shapely.Polygon(ir_)
        indices = []
        for i in tqdm(range(len(all_spatial))):
            pt = shapely.Point(all_spatial[i])
            if pt.within(poly):
                indices.append(i)
        total_indices.append(indices)

    total_indices = list(
        set([element for sublist in total_indices for element in sublist])
    )
    index_set = set(total_indices)
    dont_remove = [i for i in tqdm(range(len(all_spatial))) if i not in index_set]
    return dont_remove, index_set


def identify_spiral(adata, points, top_points, mid_points, num_points=100000):
    """
    Get the x and y coordinates evenly spaced along the entire length of the spiral.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the spatial data.
    points (list): A list of the points along the base of the basal membrane.
    top_points (list): A list of the points along the tips of villi.
    mid_points (list): A list of the points through the middle of the villi.
    num_points (int): The number of points to evenly space along the spiral.

    Returns:
    x_points_bottom: The x coordinates of the points along the base of the basal membrane.
    y_points_bottom: The y coordinates of the points along the base of the basal membrane.
    x_points_top: The x coordinates of the points along the tips of villi.
    y_points_top: The y coordinates of the points along the tips of villi.
    x_points_mid: The x coordinates of the points through the middle of the villi.
    y_points_mid: The y coordinates of the points through the middle of the villi.
    """
    spiral_main_bottom = [np.array(i) * adata.uns["unrolling_downsize"] for i in points]
    spiral_main_top = [
        np.array(i) * adata.uns["unrolling_downsize"] for i in top_points
    ]
    spiral_main_mid = [
        np.array(i) * adata.uns["unrolling_downsize"] for i in mid_points
    ]

    x_points_bottom = []
    y_points_bottom = []
    x_points_top = []
    y_points_top = []
    x_points_mid = []
    y_points_mid = []
    base_num_points = num_points
    top_num_points = base_num_points
    cumu_dist = 0
    for k in spiral_main_bottom:
        x = k[:, 0]
        y = k[:, 1]
        # Number of points you want to evenly space
        num_points = int(
            base_num_points
            * (len(k) / np.sum([np.shape(l)[0] for l in spiral_main_bottom]))
        )

        # Calculate the distances between consecutive points on the line
        distances = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)

        # Calculate the cumulative sum of distances and normalize it to [0, 1]
        cumulative_distances = np.cumsum(distances)
        cumu_dist += cumulative_distances[-1]
        normalized_distances = cumulative_distances / cumulative_distances[-1]

        # Create evenly spaced values from 0 to 1
        evenly_spaced_values = np.linspace(0, 1, num_points)

        # Use linear interpolation to calculate the x and y coordinates of the points
        x_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, x[:-1])
        )
        y_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, y[:-1])
        )
        x_points_bottom += x_points_sub
        y_points_bottom += y_points_sub

    for k in spiral_main_mid:
        x = k[:, 0]
        y = k[:, 1]

        num_points = int(
            base_num_points
            * (len(k) / np.sum([np.shape(l)[0] for l in spiral_main_mid]))
        )

        # Calculate the distances between consecutive points on the line
        distances = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)

        # Calculate the cumulative sum of distances and normalize it to [0, 1]
        cumulative_distances = np.cumsum(distances)
        normalized_distances = cumulative_distances / cumulative_distances[-1]

        # Create evenly spaced values from 0 to 1
        evenly_spaced_values = np.linspace(0, 1, num_points)

        # Use linear interpolation to calculate the x and y coordinates of the points
        x_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, x[:-1])
        )
        y_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, y[:-1])
        )
        x_points_mid += x_points_sub
        y_points_mid += y_points_sub

    for k in spiral_main_top:
        x = k[:, 0]
        y = k[:, 1]
        # Number of points you want to evenly space
        num_points = int(
            base_num_points
            * (len(k) / np.sum([np.shape(l)[0] for l in spiral_main_top]))
        )

        # Calculate the distances between consecutive points on the line
        distances = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)

        # Calculate the cumulative sum of distances and normalize it to [0, 1]
        cumulative_distances = np.cumsum(distances)
        normalized_distances = cumulative_distances / cumulative_distances[-1]

        # Create evenly spaced values from 0 to 1
        evenly_spaced_values = np.linspace(0, 1, num_points)

        # Use linear interpolation to calculate the x and y coordinates of the points
        x_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, x[:-1])
        )
        y_points_sub = list(
            np.interp(evenly_spaced_values, normalized_distances, y[:-1])
        )
        x_points_top += x_points_sub
        y_points_top += y_points_sub
    return (
        x_points_bottom,
        y_points_bottom,
        x_points_mid,
        y_points_mid,
        x_points_top,
        y_points_top,
    )


def get_distances_and_indices(
    adata,
    dont_remove,
    x_points_bottom,
    y_points_bottom,
    x_points_mid,
    y_points_mid,
    x_points_top,
    y_points_top,
    base_num_points=100000,
):
    """
    Get the distances and indices of the closest points in the top and bottom sets.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the spatial data.
    dont_remove (list): The indices of the points to keep.
    x_points_bottom (list): The x coordinates of the points along the base of the basal membrane.
    y_points_bottom (list): The y coordinates of the points along the base of the basal membrane.
    x_points_mid (list): The x coordinates of the points through the middle of the villi.
    y_points_mid (list): The y coordinates of the points through the middle of the villi.
    x_points_top (list): The x coordinates of the points along the tips of villi.
    y_points_top (list): The y coordinates of the points along the tips of villi.
    base_num_points (int): The number of points to evenly space along the spiral.

    Returns:
    all_points: The spatial points to keep.
    distances_top: The distances of the closest points in the top set.
    indices_top: The indices of the closest points in the top set.
    distances_bottom: The distances of the closest points in the bottom set.
    indices_bottom: The indices of the closest points in the bottom set.
    distances_mid: The distances of the closest points in the mid set.
    indices_mid: The indices of the closest points in the mid set.
    """
    x_points_top = np.array(x_points_top)
    y_points_top = np.array(y_points_top)
    x_points_bottom = np.array(x_points_bottom)
    y_points_bottom = np.array(y_points_bottom)
    x_points_mid = np.array(x_points_mid)
    y_points_mid = np.array(y_points_mid)

    all_points = adata.obsm["X_spatial"][dont_remove]

    # Create KD-trees for both sets of points
    tree_top = cKDTree(np.column_stack((x_points_top, y_points_top)))
    tree_mid = cKDTree(np.column_stack((x_points_mid, y_points_mid)))
    tree_bottom = cKDTree(np.column_stack((x_points_bottom, y_points_bottom)))

    # Initialize arrays to store distances and indices
    distances_top, indices_top = [], []
    distances_bottom, indices_bottom = [], []
    distances_mid, indices_mid = [], []

    # Iterate through each point in all_points
    for point in tqdm(all_points):
        # Query the top KD-tree for the closest point
        distance_top, index_top = tree_top.query(point)
        distances_top.append(distance_top)
        indices_top.append(index_top)

        # Query the bottom KD-tree for the closest point
        distance_bottom, index_bottom = tree_bottom.query(point)
        distances_bottom.append(distance_bottom)
        indices_bottom.append(index_bottom)

        # Query the bottom KD-tree for the closest point
        distance_mid, index_mid = tree_mid.query(point)
        distances_mid.append(distance_mid)
        indices_mid.append(index_mid)

    # Convert lists to numpy arrays for convenience
    distances_top = np.array(distances_top)
    indices_top = np.array(indices_top)
    distances_bottom = np.array(distances_bottom)
    indices_bottom = np.array(indices_bottom)
    distances_mid = np.array(distances_mid)
    indices_mid = np.array(indices_mid)

    # Now, distances_top, indices_top contain distances and indices of closest points in the top set,
    # and distances_bottom, indices_bottom contain distances and indices of closest points in the bottom set.
    indices_top = indices_top / base_num_points
    indices_bottom = indices_bottom / base_num_points
    indices_mid = indices_mid / base_num_points
    return (
        all_points,
        distances_top,
        indices_top,
        distances_bottom,
        indices_bottom,
        distances_mid,
        indices_mid,
    )
