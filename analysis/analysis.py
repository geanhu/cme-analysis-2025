import pandas as pd
import numpy as np
import os
from scipy import ndimage
from skimage.measure import regionprops_table
from skimage.segmentation import find_boundaries

def cells(
        cell_mask: np.ndarray,
        image: np.ndarray,
        name: str,
        result_dir: str,
        background_median = True
):
    #size
    ids, counts = np.unique(cell_mask, return_counts=True)
    cells_df = pd.DataFrame({
        'id': ids,
        'area (px)': counts
    })
    cells_df = cells_df.set_index('id')

    #intensity
    integrated_intensities = ndimage.sum_labels(image, labels=cell_mask, index=ids)
    mean_intensities = ndimage.mean(image, labels=cell_mask, index=ids)

    #background (median)
    background_mask = (cell_mask == 0)
    background = image[background_mask]
    #take median, because some bright spots are mistakenly segemented out -> inflates mean
    background_intensity = np.median(background)

    #save to df
    cells_df['integrated_intensity'] = integrated_intensities
    if background_median:
        cells_df['mean_intensity'] = [background_intensity] + list(mean_intensities[1:])
    else:
        cells_df['mean_intensity'] = mean_intensities

    #save dataframe
    cells_df.to_csv(os.path.join(result_dir, f'{name}-cells.csv'))

    return cells_df

def puncta(
        puncta_mask: np.ndarray,
        image: np.ndarray,
        name: str,
        result_dir: str,
        background_median = True,
        puncta_threshold = 2
):
    '''
    Size
    '''
    ids, counts = np.unique(puncta_mask, return_counts=True)
    puncta_df = pd.DataFrame({
        'id': ids,
        'area': counts
    })
    puncta_df = puncta_df.set_index('id')

    '''
    Intensity
    '''
    integrated_intensities = ndimage.sum_labels(image, labels=puncta_mask, index=ids)
    mean_intensities = ndimage.mean(image, labels=puncta_mask, index=ids)

    #background (median)
    background_mask = (puncta_mask == 0)
    background = image[background_mask]
    #take median, because some bright spots are mistakenly segemented out -> inflates mean
    background_intensity = np.median(background)

    #save to df
    puncta_df['integrated_intensity'] = integrated_intensities
    if background_median:
        puncta_df['mean_intensity'] = [background_intensity] + list(mean_intensities[1:])
    else:
        puncta_df['mean_intensity'] = mean_intensities

    '''
    Properties
    '''
    properties = (
        'label', 
        'perimeter', 
        'major_axis_length', 
        'minor_axis_length',
        'eccentricity'
    )
    props_dict = regionprops_table(
        puncta_mask,
        properties = properties,
    )
    props_df = pd.DataFrame(props_dict).set_index('label')
    puncta_df = pd.concat((puncta_df, props_df), axis=1)
    # Area vs perfect circle by perimeter (punishes fuzzy/spiky)
    puncta_df['circularity'] = 4 * np.pi * puncta_df['area'] / (puncta_df['perimeter']**2)
    # Area vs perfect circle by major axis (punishes ellipsity/elongation)
    puncta_df['roundness'] = 4 * puncta_df['area'] / (np.pi * puncta_df['major_axis_length']**2)

    #save
    puncta_df.to_csv(os.path.join(result_dir, f'{name}-puncta.csv'))

    return puncta_df

def puncta_cells(
    puncta_mask: np.ndarray,
    cell_mask: np.ndarray,
    puncta_df: pd.DataFrame,
    cells_df: pd.DataFrame,
    image: np.ndarray,
    name: str,
    result_dir: str
):
    #match puncta to cell
    cell_ids = np.unique(cell_mask)
    puncta_ids = np.unique(puncta_mask)
    #creates list for each puncta, where index in list = matching cell, value = overlapping px
    cell_puncta_hists = ndimage.histogram(
        cell_mask,
        min=0,
        max=cell_ids.max(),
        bins=cell_ids.max() + 1,
        labels=puncta_mask,
        index=puncta_ids
    )
    cell_puncta_hists = np.array([list(hist) for hist in cell_puncta_hists]) #list of lists -> 2D array for indexing
    cell_puncta = np.argmax(cell_puncta_hists[:, :], axis=1) #WARNING: can partition into background
    puncta_df['cell_id'] = cell_puncta

    #C_D
    cells_df['c_d'] = cells_df.apply(
        lambda row: find_c_d(row, puncta_df),
        axis = 1
    )
    
    #C_L
    cells_df['c_l'] = cells_df.apply(
        lambda row: find_c_l(row, puncta_mask, cell_mask, puncta_df, image),
        axis = 1
    )

    '''
    Distance to edge
    '''
    #find cell edges
    cell_boundaries = find_boundaries(cell_mask, mode='inner')

    #find distance of all pixels -> nearest boundary
    distance_map = ndimage.distance_transform_edt(~cell_boundaries) #flip so boundaries are 0

    #find cell centroids
    centroids = regionprops_table(puncta_mask, properties=('label', 'centroid'))
    centroids_y = centroids['centroid-0'].astype(int)
    centroids_x = centroids['centroid-1'].astype(int)

    #get distances of centroids to boundary
    distances = distance_map[centroids_y, centroids_x] # type: ignore
    distances = [np.nan] + list(distances) #undefined for background
    puncta_df['distance_to_edge'] = distances

    puncta_df.to_csv(os.path.join(result_dir, f'{name}-puncta.csv'))
    cells_df.to_csv(os.path.join(result_dir, f'{name}-cells.csv'))
    return 

def find_c_d(
    row: pd.Series,
    puncta_df: pd.DataFrame,
):
    #find puncta in cell
    cell_id = row.name
    if cell_id == 0: #skip background
        return np.nan
    puncta_in_cell = puncta_df[puncta_df['cell_id'] == cell_id]
    if puncta_in_cell.empty: #if no puncta -> exit
        return np.nan
    
    #mean
    total_puncta_int = puncta_in_cell['integrated_intensity'].sum()
    total_puncta_area = puncta_in_cell['area'].sum()
    if total_puncta_area == 0:
        return np.nan #should never do this, but prevent zero division just in case

    return total_puncta_int / total_puncta_area

def find_c_l(
    row: pd.Series,
    puncta_mask: np.ndarray,
    cell_mask: np.ndarray,
    puncta_df: pd.DataFrame,
    image: np.ndarray,
):
    #find all puncta in cell
    cell_id = row.name
    if cell_id == 0: #skip background
        return np.nan
    puncta_in_cell = puncta_df[puncta_df['cell_id'] == cell_id]
    if puncta_in_cell.empty:
        return np.nan 
    puncta_ids = list(puncta_in_cell.index)

    #mask
    cell_mask = (cell_mask == cell_id)
    puncta_mask = np.isin(puncta_mask, puncta_ids)
    cytoplasm_mask = (cell_mask & ~(puncta_mask))

    #find value
    cytoplasm = image[cytoplasm_mask]
    if cytoplasm.size == 0:
        return np.nan #prevent zero division
    cytoplasm_int = float(np.mean(cytoplasm))
    
    return cytoplasm_int