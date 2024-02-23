# import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import os
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from collections import Counter
import glob
from fast_alphashape import alphashape
import geopandas as gpd
import seaborn as sns
from shapely.ops import transform
from shapely.affinity import scale
import imageio as io
import tifffile as tiff
import imagecodecs
import shapely.affinity as sa
import cv2
import json
import plotly.graph_objs as go
import matplotlib.patches as patches
from core_functions.segmentation_evaluation import *

#### Plotting

def run_iExplorer(image_plot, if_channel, image_cmap, plot_transcripts, transcript_colors, gene_subset, pt_size, plot_segmentation, segmentation_face_color, segmentation_boundary_colors, inside_alpha, outside_alpha, continuous_vals, recalculate_shapes, finalized_adata, transcripts_df, pixel_size, max_x, min_x, max_y, min_y, xenium_dapi=None, h_an_e = None, IF_image = None):
    if plot_segmentation == 'mask':
        if continuous_vals == True:
            celltypes = []
            ids = np.array([i.split('_')[-1] for i in finalized_adata.obs.index.values]).astype(int)
            id_df = pd.DataFrame(zip(ids, finalized_adata.obs[segmentation_face_color].values), columns=['id', segmentation_face_color])
            transcripts_with_obs = transcripts_df.merge(id_df, left_on='split_cell', right_on='id', how='left')
            transcripts_with_obs = transcripts_with_obs.dropna(axis=0)
            

            print('Making Shapes')
            gby = transcripts_with_obs[(transcripts_with_obs.split_cell != 0) & (transcripts_with_obs.split_cell != -1)].groupby("split_cell")
            if type(recalculate_shapes) == bool:
                shapes = []
                pts = []
                for group in tqdm(gby):
                    shapes.append(make_alphashape(group[1][['x', 'y']].values, alpha=0.05))
                    pts.append(np.mean(group[1][['x', 'y']].values, axis=0))
                    celltypes.append(group[1][segmentation_face_color].values[0])
                shapes = gpd.GeoSeries(shapes)
                pts = np.array(pts)
            else:
                for group in tqdm(gby):
                    celltypes.append(group[1][segmentation_face_color].values[0])
                shapes = recalculate_shapes
                print('If you get error try setting recalculate_shapes to True')

            from matplotlib.colors import Normalize
            from matplotlib.cm import coolwarm

            # Generate an example array of numbers (replace this with your own data)
            data = np.array(celltypes)

            # Define the colormap and normalization
            cmap = coolwarm
            norm = Normalize(vmin=data.min(), vmax=data.max())

            # Create a colormap object
            mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

            # Map the data to colors
            colors = mappable.to_rgba(data)
        else:
            celltypes = []
            ids = np.array([i.split('_')[-1] for i in finalized_adata.obs.index.values]).astype(int)
            id_df = pd.DataFrame(zip(ids, finalized_adata.obs[segmentation_face_color].values), columns=['id', segmentation_face_color])
            transcripts_with_obs = transcripts_df.merge(id_df, left_on='split_cell', right_on='id', how='left')
            transcripts_with_obs = transcripts_with_obs.dropna(axis=0)
            

            print('Making Shapes')
            gby = transcripts_with_obs[(transcripts_with_obs.split_cell != 0) & (transcripts_with_obs.split_cell != -1)].groupby("split_cell")
            

            if type(recalculate_shapes) == bool:
                shapes = []
                for group in tqdm(gby):
                    shapes.append(make_alphashape(group[1][['x', 'y']].values, alpha=0.05))
                    ctype = group[1][segmentation_face_color].values[0]
                    cell_location = np.where(finalized_adata.obs[segmentation_face_color].cat.categories == ctype)[0] 
                    try:
                        celltypes.append(finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location][0])
                    except:
                        celltypes.append(finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location[0]])
                shapes = gpd.GeoSeries(shapes)
            else:
                shapes = recalculate_shapes
                for group in tqdm(gby):
                    ctype = group[1][segmentation_face_color].values[0]
                    cell_location = np.where(finalized_adata.obs[segmentation_face_color].cat.categories == ctype)[0] 
                    try:
                        celltypes.append(finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location][0])
                    except:
                        celltypes.append(finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location[0]])
                print('If you get error try setting recalculate_shapes to True')
            colors = celltypes

        def scale_to_image(x, y):
            return(x/pixel_size, y/pixel_size)
            
        
        print('Plotting')
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        if image_plot == 'xenium_dapi':
            img_cropped = xenium_dapi[min_x:max_x, min_y:max_y]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
        elif image_plot == 'h_and_e':
            img_cropped = h_an_e[min_x:max_x, min_y:max_y]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
        elif image_plot == 'IF':
            img_cropped = IF_image[min_x:max_x, min_y:max_y, if_channel]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
        # Create an empty GeoDataFrame to store adjusted polygons
        adjusted_shapes = []

        # Iterate through the shapes DataFrame and adjust each polygon
        for original_polygon in shapes:
            scaled_polygon = sa.translate(original_polygon, -min_y, -min_x)
            adjusted_shapes.append(scaled_polygon)

        adjusted_shapes = gpd.GeoSeries(adjusted_shapes)

        for geometry, color in zip(adjusted_shapes, colors):
            if geometry.geom_type == 'Polygon':
                patch = plt.Polygon(list(zip(*geometry.exterior.xy)), facecolor=color, edgecolor='none', alpha=inside_alpha, zorder=1)
                ax.add_patch(patch)
            elif geometry.geom_type == 'MultiPolygon':
                for poly in geometry:
                    patch = plt.Polygon(list(zip(*poly.exterior.xy)), facecolor=color, edgecolor='none', alpha=inside_alpha, zorder=1)
                    ax.add_patch(patch)

        # Plot polygon edges with edgecolor based on data values
        for geometry, color in zip(adjusted_shapes, colors):
            if geometry.geom_type == 'Polygon':
                ax.plot(*geometry.exterior.xy, color=color, alpha=outside_alpha)
            elif geometry.geom_type == 'MultiPolygon':
                for poly in geometry:
                    ax.plot(*poly.exterior.xy, color=color, alpha=outside_alpha)


        if plot_transcripts == True:
            transcripts_genes_only = transcripts_df[transcripts_df['gene'].isin(gene_subset)]
            col_ct = 0
            for i in gene_subset:
                transcripts_genes_only_current = transcripts_genes_only[transcripts_genes_only['gene'] == i]
                for x, y in zip(transcripts_genes_only_current.x.values, transcripts_genes_only_current.y.values):
                    circle = patches.Circle((x - min_y, y - min_x), radius=pt_size, edgecolor='black', linewidth=1, facecolor=transcript_colors[col_ct], alpha=1,zorder=2)
                    ax.add_patch(circle)
                col_ct += 1

            col_ct = 0
            for i in gene_subset:
                plt.scatter([], [], c=transcript_colors[col_ct], label=i)
                col_ct += 1
        plt.legend()

        plt.xlim(0, max_y - min_y)
        plt.ylim(0, max_x - min_x)
        plt.axis('equal')
        plt.show()
        return shapes
    elif plot_segmentation == 'transcripts':
        if continuous_vals == True:
            celltypes = []
            ids = np.array([i.split('_')[-1] for i in finalized_adata.obs.index.values]).astype(int)
            id_df = pd.DataFrame(zip(ids, finalized_adata.obs[segmentation_face_color].values), columns=['id', segmentation_face_color])
            transcripts_with_obs = transcripts_df.merge(id_df, left_on='split_cell', right_on='id', how='left')
            transcripts_with_obs = transcripts_with_obs.dropna(axis=0)
            

            print('Making Shapes')
            gby = transcripts_with_obs[(transcripts_with_obs.split_cell != 0) & (transcripts_with_obs.split_cell != -1)].groupby("split_cell")
            

            shapes = []
            xpts = []
            ypts = []
            for group in tqdm(gby):
                xpts.extend(group[1][['x', 'y']].values[:, 0])
                ypts.extend(group[1][['x', 'y']].values[:, 1])
                celltypes.extend(group[1][segmentation_face_color].values)
            xpts = np.array(xpts)
            ypts = np.array(ypts)

            from matplotlib.colors import Normalize
            from matplotlib.cm import coolwarm

            # Generate an example array of numbers (replace this with your own data)
            data = np.array(celltypes)

            # Define the colormap and normalization
            cmap = coolwarm
            norm = Normalize(vmin=data.min(), vmax=data.max())

            # Create a colormap object
            mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

            # Map the data to colors
            colors = mappable.to_rgba(data)
        else:
            if continuous_vals:
                celltypes = []
                ids = np.array([i.split('_')[-1] for i in finalized_adata.obs.index.values]).astype(int)
                id_df = pd.DataFrame(zip(ids, finalized_adata.obs[segmentation_face_color].values), columns=['id', segmentation_face_color])
                transcripts_with_obs = transcripts_df.merge(id_df, left_on='split_cell', right_on='id', how='left')
                transcripts_with_obs = transcripts_with_obs.dropna(axis=0)
                

                print('Making Shapes')
                gby = transcripts_with_obs[(transcripts_with_obs.split_cell != 0) & (transcripts_with_obs.split_cell != -1)].groupby("split_cell")
                

                xpts = []
                ypts = []
                for group in tqdm(gby):
                    xpts.extend(group[1][['x', 'y']].values[:, 0])
                    ypts.extend(group[1][['x', 'y']].values[:, 1])
                    celltypes.extend(group[1][segmentation_face_color].values)
                xpts = np.array(xpts)
                ypts = np.array(ypts)
                colors = celltypes
            else:
                celltypes = []
                ids = np.array([i.split('_')[-1] for i in finalized_adata.obs.index.values]).astype(int)
                id_df = pd.DataFrame(zip(ids, finalized_adata.obs[segmentation_face_color].values), columns=['id', segmentation_face_color])
                transcripts_with_obs = transcripts_df.merge(id_df, left_on='split_cell', right_on='id', how='left')
                transcripts_with_obs = transcripts_with_obs.dropna(axis=0)
                

                print('Making Shapes')
                gby = transcripts_with_obs[(transcripts_with_obs.split_cell != 0) & (transcripts_with_obs.split_cell != -1)].groupby("split_cell")
                
            
                xpts = []
                ypts = []
                for group in tqdm(gby):
                    xpts.extend(group[1][['x', 'y']].values[:, 0])
                    ypts.extend(group[1][['x', 'y']].values[:, 1])
                    ctype = group[1][segmentation_face_color].values[0]
                    cell_location = np.where(finalized_adata.obs[segmentation_face_color].cat.categories == ctype)[0] 
                    try:
                        to_append = finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location][0]
                    except:
                        to_append = finalized_adata.uns[f'{segmentation_face_color}_colors'][cell_location[0]]
                    celltypes.extend([to_append for j in range(len(group[1][segmentation_face_color].values))])
                xpts = np.array(xpts)
                ypts = np.array(ypts)
                colors = celltypes

        def scale_to_image(x, y):
            return(x/pixel_size, y/pixel_size)
            
        
        print('Plotting')
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        if image_plot == 'dapi':
            img_cropped = xenium_dapi[min_x:max_x, min_y:max_y]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
        elif image_plot == 'h_and_e':
            img_cropped = h_an_e[min_x:max_x, min_y:max_y]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
        elif image_plot == 'IF':
            img_cropped = IF_image[min_x:max_x, min_y:max_y, if_channel]
            ax.imshow(
                img_cropped,
                vmax=np.percentile(img_cropped, 99.9),
                cmap=image_cmap
            )
            
        plt.scatter(xpts - min_y, ypts - min_x, c=colors, s=1, linewidths=0.1)


        if plot_transcripts == True:
            transcripts_genes_only = transcripts_df[transcripts_df['gene'].isin(gene_subset)]
            col_ct = 0
            for i in gene_subset:
                transcripts_genes_only_current = transcripts_genes_only[transcripts_genes_only['gene'] == i]
                for x, y in zip(transcripts_genes_only_current.x.values, transcripts_genes_only_current.y.values):
                    circle = patches.Circle((x - min_y, y - min_x), radius=pt_size, edgecolor='black', linewidth=1, facecolor=transcript_colors[col_ct], alpha=1,zorder=2)
                    ax.add_patch(circle)
                col_ct += 1

            col_ct = 0
            for i in gene_subset:
                plt.scatter([], [], c=transcript_colors[col_ct], label=i)
                col_ct += 1
        plt.legend()

        plt.xlim(0, max_y - min_y)
        plt.ylim(0, max_x - min_x)
        plt.axis('equal')
        plt.show()
    