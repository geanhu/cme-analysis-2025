import argparse
from pathlib import Path
import roifile
from matplotlib import path
import zipfile

import imagej
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scyjava as sj

def exp_model(t, a, k, c):
    return a * np.exp(-k * t) + c

def roi_to_mask(shape, roi_data):
    roi = roifile.roiread(roi_data)
    vertices = np.vstack((roi.rois[0].x, roi.rois[0].y)).T #type:ignore
    y, x = np.mgrid[:shape[0], :shape[1]]
    points = np.vstack((x.flatten(), y.flatten())).T
    mask_path = path.Path(vertices)
    mask = mask_path.contains_points(points).reshape(shape) #type:ignore
    return mask

def zip_to_mask(shape, path):
    mask_list = []
    with zipfile.ZipFile(path, 'r') as z:
        for filename in z.namelist():
            with z.open(filename) as f:
                roi_data = f.read()
                mask_list.append(roi_to_mask(shape, roi_data))
    return mask_list

def frap(
    files: list,
    ij,
    bleach_frame: int,
    frame_time: float
):
    importer_options = sj.jimport('loci.plugins.in.ImporterOptions')
    options = importer_options()
    options.setSplitChannels(True)
    bf = sj.jimport('loci.plugins.BF')

    #to store parameters
    df = pd.DataFrame()
    mobile_fraction = []
    t_half = []

    '''
    Get ROIs & image data
    '''
    for file in files:
        print(f'Analyzing file {file} ...')
        #Open
        options.setId(file)
        channels = bf.openImagePlus(options)
        image = None
        for channel in channels:
            if "C=1" in channel.getTitle():
                image = channel
        image_data = ij.py.from_java(ij.convert().convert(image, ij.py.get_array_type(image)))
        
        #Get ROIs
        name = Path(file).name
        shape = image_data.shape
        bleach_rois = zip_to_mask(shape, f'{name}-bleach.zip')
        cell_rois = zip_to_mask(shape, f'{name}-cells.zip')
        background_roi = zip_to_mask(shape, f'{name}-background.zip')[0]

        '''
        Analyze & fit model
        '''
        roi_count = 0
        for bleach_mask, cell_mask in zip(bleach_rois, cell_rois):
            image_data = np.transpose(image_data, (2, 0, 1)) #(t, y, x) for easier iterating
            num_frames = image_data.shape[0]

            #Average ROIs
            cell_intensity = np.array([np.mean(frame[cell_mask]) for frame in image_data])
            bleach_intensity = np.array([np.mean(frame[bleach_mask]) for frame in image_data])
            background_intensity = np.array([np.mean(frame[background_roi])] for frame in image_data)
            cell_intensity = cell_intensity - background_intensity
            bleach_intensity = bleach_intensity - background_intensity

            #Calculate normalization constants
            bleach_pre = np.mean(bleach_intensity[:bleach_frame])
            cell_pre = np.mean(cell_intensity[:bleach_frame])

            #Normalize
            bleach_intensity = bleach_intensity / bleach_pre
            cell_intensity = cell_intensity / cell_pre
            fluor_norm = bleach_intensity / cell_intensity

            #Curve fit
            time_series = np.arange(num_frames) * frame_time
            time_aligned = time_series - (bleach_frame * frame_time)
            init = [fluor_norm[bleach_frame], 0.1, fluor_norm[-1]]
            params, _ = curve_fit(
                exp_model,
                time_aligned[bleach_frame:],
                fluor_norm[bleach_frame:],
                p0 = init,
                maxfev = 5000
            )
            _, k, c = params
            mobile_fraction.append(c)
            t_half.append(np.log(2) / k)

            '''Store data'''
            name = f'{Path(file).name}_{roi_count}'
            #Curve
            df[name] = pd.Series(fluor_norm, index=time_aligned)

            roi_count += 1
    
    print(
        f'Mobile fraction: {np.mean(mobile_fraction)} +/- {np.std(mobile_fraction)}'
    )
    print(
        f't_half: {np.mean(t_half)} +/- {np.std(t_half)}'
    )
    return df

def main():
    '''
    Parse Arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input',
        help = 'Input directory',
        type = str
    )
    parser.add_argument(
        'output',
        help = 'Output directory',
        type = str
    )
    args = parser.parse_args()

    '''
    Initialize ImageJ
    '''
    print('Initializing ImageJ ...')
    ij = imagej.init(
        ij_dir_or_version_or_endpoint='/Applications/Fiji.app',
        mode = 'headless',
        add_legacy = True,
    )

    '''
    Parse input files
    '''
    frap_files = []
    for file in Path(args.input).iterdir():
        if file.suffix == '.czi':
            frap_files.append(str(file))

    '''
    Params
    '''
    bleach_frame = 11
    frame_time = 1
    df = frap(frap_files, ij, bleach_frame, frame_time)
    df.to_csv('{args.output}/frap.csv')
    print('Saved analysis results')

    #kill
    sj.jimport('java.lang.System').exit(0)

    return 0

if __name__ == "__main__":
    main()