#internal
from utils import preprocess
from utils import io
from analysis import analysis

#external
import os
import skimage.io
from scipy.ndimage import label
import warnings
import numpy as np
from cellpose import models
#---for stand-alone behavior
import argparse
from pathlib import Path
import imagej
import scyjava as sj
import pandas as pd

def cellpose(
        input: str,
        result_dir: str,
        ij,
        name: str,
        save_intermediate: bool,
        exclude_slices: int = 0,
):

    '''
    Open image and read channel colors to find channel to segment
    '''
    #Find channel to segment
    channel_order = io.parse_channels(input)
    #3. fall back to first channel
    channel_to_segment_idx = 0 

    subtract_background = True
    photobleaching_correction = False #TODO: fix
    median = True
    erosion = False #TODO: fix

    #1. Use transilluminated channel
    if 'T' in channel_order: 
        channel_to_segment_idx = channel_order.index('T')
        subtract_background = False
        photobleaching_correction = False
        median = False
        erosion = False
    #2. Use GFP channel
    elif 'G' in channel_order: #TODO: temp
        channel_to_segment_idx = channel_order.index('G')
    
    #Open image and extract target channel
    channel_to_segment = io.open_bioformats(input)[channel_to_segment_idx]

    intermediate_dir = None
    if save_intermediate:
        intermediate_dir = os.path.join(result_dir, 'cellpose-preprocessing')
        os.makedirs(intermediate_dir, exist_ok=True)

    '''
    Preprocess channel image for better Cellpose segmentation
    '''
    if intermediate_dir:
        channel_to_segment = preprocess.cellpose(
            channel = channel_to_segment,
            ij = ij,
            output_dir = intermediate_dir,
            save_intermediate = True,
            run_subtract_background = subtract_background,
            run_photobleaching_correction = photobleaching_correction,
            exclude_slices=exclude_slices,
            run_median = median,
            run_erosion = erosion
        )
    else:
        channel_to_segment = preprocess.cellpose(
            channel = channel_to_segment,
            ij = ij,
            save_intermediate = False,
            run_subtract_background = subtract_background,
            run_photobleaching_correction = photobleaching_correction,
            exclude_slices=exclude_slices,
            run_median = median,
            run_erosion = erosion
        )

    '''
    Segment channel
    '''
    #convert back to python formats
    channel_to_segment = ij.py.from_java(channel_to_segment).values
    assert len(channel_to_segment.shape) == 2, channel_to_segment.shape

    #setup Cellpose
    model = models.Cellpose(gpu=True, model_type='cyto')

    #segment
    diameter = 50
    masks, flows, styles, diams = model.eval([channel_to_segment], diameter=diameter, channels=[[0, 0]])
    mask = preprocess.fill_holes(masks[0], ij)
    num_cells = np.unique(mask).shape[0]
    
    #save mask
    #Suppress low contrast warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    skimage.io.imsave(os.path.join(result_dir, f'{name}-cells-mask.tif'), mask)

    '''
    Return puncta image
    '''
    puncta_channel_idx = 0
    if 'G' in channel_order:
        puncta_channel_idx = channel_order.index('G')
    puncta_channel = io.open_bioformats(input)[puncta_channel_idx]
    puncta_channel = preprocess.z_proj(
        puncta_channel,
        exclude_slices = 3
    )
    puncta_channel = ij.py.from_java(puncta_channel).values
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mip.tif'), puncta_channel)
    
    return mask, puncta_channel, num_cells

def puncta(
        input: str,
        result_dir: str,
        name: str,
        ij,
        threshold: float = 0.0
):
    #Suppress low contrast warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    '''
    Open fluorescence channel containing puncta
    '''
    #find puncta channel
    channel_order = io.parse_channels(input)
    #2. fall back to first channel
    channel_to_segment_idx = 0
    #1. Use GFP
    if 'G' in channel_order:
        channel_to_segment_idx = channel_order.index('G')
    #Open image and extract target channel
    channel_to_segment = io.open_bioformats(input)[channel_to_segment_idx]

    '''
    Preprocess to enhance puncta contrast
    '''
    channel_to_segment = preprocess.z_proj(
        channel_to_segment,
        exclude_slices = 3
    )
    
    '''
    TODO: this is not fair comparison

    Problem:
    when finding threshold
    - median -> Renyi -- method too stringent
    - median -> Otsu -- creates weird shapes because not stringent enough
    when applying threshold
    - subtraction -> Z-score threshold -- method not stringent enough, admits noise
    '''
    if threshold == 0:
        channel_enhanced = preprocess.subtract_background(
            channel_to_segment
        )
        channel_enhanced = preprocess.gaussian_subtraction(
            channel_enhanced,
            ij
        )
    else:
        channel_enhanced = preprocess.median(
            channel_to_segment,
            ij,
            stack = False
        )
    
    '''
    Threshold puncta in each cell
    '''
    #Threshold -> creates 8-bit 0/255 mask
    channel_enhanced.show()
    if threshold == 0.0:
        print('Finding threshold ...')
        macro_str='''
        setAutoThreshold("RenyiEntropy dark no-reset");
        setOption("BlackBackground", true);
        ''' 
        ij.py.run_macro(macro_str)
        threshold = channel_enhanced.getProcessor().getMinThreshold()
    else:
        upper = channel_enhanced.getProcessor().getMaxThreshold()
        #calculate threshold from z-score
        stats = channel_enhanced.getStatistics()
        threshold = int(threshold * stats.stdDev + stats.mean)
        print(threshold) #sanity check setting threshold
        ij.IJ.setThreshold(channel_enhanced, threshold, upper)
    
    #find threshold z
    stats = channel_enhanced.getStatistics()
    if stats.stdDev > 0:
        threshold_z = (threshold - stats.mean) / stats.stdDev
    else:
        threshold_z = threshold

    #save label mask
    macro_str = 'run("Convert to Mask");'
    ij.py.run_macro(macro_str)
    channel_enhanced = ij.IJ.getImage()
    threshold_mask = ij.py.from_java(channel_enhanced).values
    channel_enhanced.close()
    labeled_array, num_objects = label(threshold_mask) # type: ignore
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mask.tif'), labeled_array)

    return labeled_array, num_objects, threshold_z

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input',
        help = "path to input file",
        type = str
    )
    parser.add_argument(
        'output_dir',
        help = "path to output directory",
        type = str
    )
    parser.add_argument(
        '--threshold',
        type = int,
        default = 0
    )
    args = parser.parse_args()

    #init ij
    print("Initializing ImageJ ...")
    ij = imagej.init(
        ij_dir_or_version_or_endpoint='/Applications/Fiji.app',
        mode = "headless",
        add_legacy = True,
    )

    #segment puncta
    print("Segmenting puncta ...")
    name = Path(args.input).stem.replace('.', '')
    mask, num_puncta, threshold = puncta(
        args.input,
        args.output_dir,
        name,
        ij,
        args.threshold
    )
    print(f'{num_puncta} puncta segmented using threshold {threshold}')

    #analyze
    print("Analyzing puncta ...")
    puncta_image = skimage.io.imread(f'{args.output_dir}/{name}-puncta-mip.tif')
    cells_df = pd.read_csv(f'{args.output_dir}/{name}-cells.csv')
    puncta_df = analysis.puncta(
        puncta_mask = mask, # type: ignore
        image = puncta_image,
        name = name,
        result_dir = args.output_dir
    )
    analysis.puncta_cells(
        puncta_mask = mask, # type: ignore
        cell_mask = skimage.io.imread(f'{args.output_dir}/{name}-cells-mask.tif'),
        image = puncta_image,
        puncta_df = puncta_df,
        cells_df = cells_df,
        name = name,
        result_dir = args.output_dir
    )

    #force shutdown to kill waiting workers
    #TODO: there is probably a cleaner solution to this
    sj.jimport('java.lang.System').exit(0)
    
    return

if __name__ == "__main__":
    main()