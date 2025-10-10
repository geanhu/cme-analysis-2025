import os
import scyjava as sj
import numpy as np

def cellpose(
        channel,
        ij,
        output_dir = None,
        save_intermediate: bool = False,
        run_subtract_background: bool = True,
        run_photobleaching_correction: bool = True,
        exclude_slices = 0,
        run_median: bool = True,
        run_erosion: bool = True
):
    step_counter = 1

    #1. Subtract background
    if run_subtract_background:
        channel = subtract_background(channel)
        if save_intermediate and output_dir:
            ij.IJ.saveAsTiff(channel, os.path.join(output_dir, f"{step_counter}-subtract-background.tif"))
        step_counter += 1

    #2. Correct photobleaching
    if run_photobleaching_correction:
        corrected = photobleaching_correction(channel, ij)
        if corrected:
            channel = corrected
            if save_intermediate and output_dir:
                ij.IJ.saveAsTiff(channel, os.path.join(output_dir, f"{step_counter}-photobleaching-correction.tif"))
            step_counter += 1

    #3. Z Proj, sum
    channel = z_proj(channel, method='sum', exclude_slices=exclude_slices)
    if save_intermediate and output_dir:
        ij.IJ.saveAsTiff(channel, os.path.join(output_dir, f"{step_counter}-sum-projection.tif"))
    step_counter += 1

    #4. Median = 5
    if run_median:
        channel = median(channel, ij, size=5)
        if save_intermediate and output_dir:
            ij.IJ.saveAsTiff(channel, os.path.join(output_dir, f"{step_counter}-median.tif"))
        step_counter += 1

    #5. Erosion disk rad = 4
    if run_erosion:
        channel = erosion(channel, ij, radius=4)
        if save_intermediate and output_dir:
            ij.IJ.saveAsTiff(channel, os.path.join(output_dir, f"{step_counter}-erosion.tif"))
        step_counter += 1

    return channel

def erosion(image, ij, radius=4):
    image.show()

    macro_str = f"""
    run("Morphological Filters", "operation=Erosion element=Disk radius={radius}")
    """
    ij.py.run_macro(macro_str)

    erosion_image = ij.IJ.getImage()
    image.close()

    return erosion_image

def median(image, ij, size=5, stack=True):
    #register image
    image.show()

    #stack
    if stack:
        stack_str = ' stack'
    else:
        stack_str = ''

    #run macro
    macro_str = f"""
    run("Median...", "radius={size}{stack_str}");
    """
    ij.py.run_macro(macro_str)

    #get new image
    median_image = ij.IJ.getImage()
    image.close()

    return median_image

def subtract_background(image, radius=70.0, paraboloid=False):
    '''
    Rolling ball background subtraction on active image

    Set radius = widest cell diameter in pixels / 2
    '''
    for slice in range(1, image.getStackSize() + 1):
        #Set slice, to avoid "junk" data being dumped on first slice
        image.setSlice(slice)

        #Import BackgroundSubtracter() Java object
        BackgroundSubtracter = sj.jimport('ij.plugin.filter.BackgroundSubtracter')

        #Create instance
        bs = BackgroundSubtracter()

        #Get processor from ImagePlus for compatibility with function input
        ip = image.getProcessor()

        #Subtract background
        bs.rollingBallBackground(
            ip, #processor
            radius, #rolling ball radius
            False, #create background instead
            False, #background is light and objects are dark
            paraboloid, #use paraboloid function
            False, #pre-smooth image with 3x3 kernel
            False #prevent removing of objects in corners
        )

    return image

def photobleaching_correction(image, ij, method='Exponential Fitting'):
    #TODO: I don't know if the test image has no decay or if this doesn't work
    #register image
    image.show()
    
    image.setSlice(image.getStackSize())
    before = image.getStatistics().mean

    #run macro
    macro_str = f"""
    'run("Bleach Correction", "correction=[{method}]");
    """
    ij.py.run_macro(macro_str)

    #get new image
    corrected_image = ij.IJ.getImage()
    
    corrected_image.setSlice(corrected_image.getStackSize())
    after = corrected_image.getStatistics().mean
    image.close()
    if after - before < 0.1:
        return None

    return corrected_image

def z_proj(image, method = 'max', exclude_slices = 0):
    #Import ZProjector object
    ZProjector = sj.jimport('ij.plugin.ZProjector')

    #Create instance
    zp = ZProjector()

    #Z Project
    image = zp.run(
        image,
        method,
        int(1 + exclude_slices),
        int(image.getStackSize() - exclude_slices)
    )
    #TODO: if T > 1, Add " all" to 'method' to project all hyperstack time points.
    
    return image

def fill_holes(
    mask: np.ndarray,
    ij
) -> np.ndarray:
    #send to ij
    mask_image = ij.IJ.createImage(
        'mask',
        '16-bit black',
        mask.shape[1],
        mask.shape[0],
        1
    )
    mask_image.getProcessor().setIntArray(mask.astype(np.int32))

    #fill holes
    mask_image.show()
    ij.py.run_macro('run("Fill Label Holes", "background=4 labeling=[16 bits]");')
    mask_filled = ij.IJ.getImage()
    mask_image.close()

    #send to python
    mask_filled = ij.py.from_java(mask_filled).values
    mask_filled = mask_filled.T #fix axis issues

    return mask_filled

def gaussian_subtraction(image, ij, sigma=10):
    #duplicate
    image.show()

    #blur
    blurred = image.duplicate()
    blurred.getProcessor().blurGaussian(sigma)

    #subtract
    ImageCalculator = sj.jimport('ij.plugin.ImageCalculator')
    ic = ImageCalculator()
    enhanced = ic.run(
        image,
        blurred,
        'subtract'
    )
    image.close()
    blurred.close()

    return enhanced