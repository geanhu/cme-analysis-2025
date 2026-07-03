#internal
from utils import preprocess
from utils import io
from analysis import analysis

#external
import os
import skimage.io
from scipy.ndimage import label
import scipy
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.morphology import remove_small_objects
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
        projection = 'sum',
        crop = 100
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
    erosion = False

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

    #crop
    channel_to_segment = preprocess.crop(channel_to_segment, amount = crop)

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
    puncta_channel_original = io.open_bioformats(input)[puncta_channel_idx]

    # Max project
    puncta_channel = preprocess.z_proj(
        puncta_channel_original,
        exclude_slices = 3
    )
    puncta_channel_max = ij.py.from_java(puncta_channel).values
    puncta_channel_max = preprocess.crop(puncta_channel_max, amount = crop)
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mip.tif'), puncta_channel_max)

    # Sum project
    puncta_channel = preprocess.z_proj(
        puncta_channel_original,
        method = 'sum',
        exclude_slices = 3
    )
    puncta_channel_sum = ij.py.from_java(puncta_channel).values
    puncta_channel_sum = preprocess.crop(puncta_channel_sum, amount = crop)
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-sum.tif'), puncta_channel_sum)

    if projection == 'sum':
        puncta_channel = puncta_channel_sum
    else:
        puncta_channel = puncta_channel_max
    
    return mask, puncta_channel, num_cells

'''
!! Not used
'''
def puncta(
        input: str,
        result_dir: str,
        name: str,
        ij,
        threshold: float = 0.0,
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
    channel_enhanced = channel_to_segment
    '''
    channel_enhanced = preprocess.subtract_background(
            channel_to_segment
        )
    channel_enhanced = preprocess.gaussian_subtraction(
        channel_enhanced,
        ij
    )
    '''
    
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

        #find threshold z
        stats = channel_enhanced.getStatistics()
        if stats.stdDev > 0:
            threshold_z = (threshold - stats.mean) / stats.stdDev
        else:
            threshold_z = threshold
        
        #mask
        macro_str = 'run("Convert to Mask");'
        ij.py.run_macro(macro_str)

        #convert
        threshold_mask = ij.py.from_java(channel_enhanced).values
        channel_enhanced.close()
    else:
        threshold_z = threshold
        upper = channel_enhanced.getProcessor().getMaxThreshold()
        #calculate threshold from z-score
        stats = channel_enhanced.getStatistics()
        threshold = int(threshold * stats.stdDev + stats.mean)
        print(f'Using threshold : {threshold}') #sanity check setting threshold
        #set
        threshold_mask = ij.py.from_java(channel_enhanced)
        threshold_mask = ((threshold_mask >= threshold) * 255).astype(np.uint8)
    
    #save label mask
    labeled_array, num_objects = label(threshold_mask) # type: ignore
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mask.tif'), labeled_array)

    return labeled_array, num_objects, threshold_z

def puncta_local(
        result_dir: str,
        name: str,
        cell_mask,
        log_threshold : float = 3.5,
        watershed_threshold : float = 3.0,
        projection = 'sum'
):
    #open image
    puncta_image = skimage.io.imread(f'{result_dir}/{name}-puncta-mip.tif')
    puncta_image = puncta_image.astype(float) #convert to float to prevent underflow during background subtraction
    
    #subtract background (cell 0) average
    background_med = np.median(puncta_image[cell_mask == 0]) #median in case any puncta in background
    puncta_image_sb = puncta_image - background_med
    
    #segment puncta only using area of each cell
    puncta_mask = np.zeros_like(cell_mask) #initialize puncta mask as all background
    puncta_num = 1

    #pre-calculate LoG for entire image
    max_response = np.zeros_like(puncta_image_sb)
    for sigma in range(2, 6): #estimate puncta radius 2-5 px
            response = -1 * scipy.ndimage.gaussian_laplace(puncta_image_sb, sigma = sigma)
            response *= (sigma ** 2) #normalization to compare between sigmas
            max_response = np.maximum(max_response, response) #keep best response

    #segment each cell
    for num_cell in np.unique(cell_mask)[1:]: #skip 0, is background
        #1. Get cell area
        num_cell_mask = (cell_mask == num_cell) #boolean mask
        num_cell_puncta_image_sb = puncta_image_sb[num_cell_mask]
        num_cell_puncta_image = puncta_image[num_cell_mask]
        #
        num_cell_med = np.median(num_cell_puncta_image)
        num_cell_mad = scipy.stats.median_abs_deviation(num_cell_puncta_image, scale='normal')

        #2. Identify puncta centers by LoG
        num_cell_response = max_response[num_cell_mask]
        log_threshold_cell = np.median(num_cell_response) + (log_threshold * scipy.stats.median_abs_deviation(num_cell_response, scale='normal'))
        log_peaks_mask = (max_response * num_cell_mask) > log_threshold_cell
        seeds = peak_local_max(
            max_response * num_cell_mask,
            min_distance = 5, #peaks must be 5 px apart
            labels=log_peaks_mask #must be at least >2 std response
        )

        #3. Seeded watershed
        markers = np.zeros_like(puncta_image, dtype=int)
        for i, (row, col) in enumerate(seeds):
            markers[row, col] = puncta_num + i
        watershed_limit = (puncta_image > (num_cell_med + watershed_threshold * num_cell_mad)) & num_cell_mask
        segmented_puncta = watershed(-puncta_image, markers, mask=watershed_limit)

        #5. Store into puncta_mask
        puncta_mask[num_cell_mask] = segmented_puncta[num_cell_mask]
        puncta_num += len(seeds)

    #4. Throw away small areas (noise)
    puncta_mask = remove_small_objects(puncta_mask, min_size=5)
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mask.tif'), puncta_mask)

    return puncta_mask, np.unique(puncta_mask).shape[0] - 1

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