from cellpose import models
from utils import preprocess
from utils import io
import os
import skimage.io
from scipy.ndimage import label
import warnings

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
    
    return mask, puncta_channel

def puncta(
        input: str,
        result_dir: str,
        name: str,
        mask,
        ij,
):
    #Suppress low contrast warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    #can also provide path to mask -> will open as np array
    if type(mask) == str:
        mask = skimage.io.imread(mask)

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

    channel_enhanced = preprocess.subtract_background(
        channel_to_segment
    )
    channel_enhanced = preprocess.gaussian_subtraction(
        channel_enhanced,
        ij
    )
        
    '''
    Threshold puncta in each cell
    '''
    #Threshold -> creates 8-bit 0/255 mask
    channel_enhanced.show()
    macro_str='''
    setAutoThreshold("RenyiEntropy dark no-reset");
    setOption("BlackBackground", true);
    run("Convert to Mask");
    '''
    ij.py.run_macro(macro_str)

    #save label mask
    threshold_mask = ij.py.from_java(channel_enhanced).values
    channel_enhanced.close()
    labeled_array, num_objects = label(threshold_mask) # type: ignore
    skimage.io.imsave(os.path.join(result_dir, f'{name}-puncta-mask.tif'), labeled_array)

    return labeled_array, num_objects

def main():
    return

if __name__ == "main":
    main()