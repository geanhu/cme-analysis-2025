#outside imports
import argparse
import time
import scyjava as sj
import os
import imagej
import warnings
from pathlib import Path
import numpy as np
from tqdm import tqdm

#module imports
from analysis import segment
from analysis import analysis

def main():
    #Suppress low contrast warnings
    warnings.filterwarnings("ignore", category=UserWarning)

    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_dir',
        help = "Input dir",
        type = str
    )
    parser.add_argument(
        '-s', '--save_intermediate',
        help = 'Save results from all intermediate steps',
        type = bool,
        default = False
    )
    parser.add_argument(
        '-l', '--log_threshold',
        help = 'Threshold (modified z-score) for puncta center detection after LoG filter',
        type = float,
        default = 3.5
    )
    parser.add_argument(
        '-w', '--watershed_threshold',
        help = 'Threshold (modified z-score) for outer limit of watershed segmentation',
        type = float,
        default = 3.0
    )
    parser.add_argument(
        '-p', '--puncta_projection',
        help = 'Projection mode for fluorescence channel (sum or maximum)',
        type = str,
        default = 'sum'
    )
    parser.add_argument(
        '-c', '--crop',
        help = 'Number of pixels to crop from input image (one side, i.e. will crop 2*this number in each dimension)',
        type = int,
        default = 50
    )
    args = parser.parse_args()

    #initialize ij
    print('Initializing ImageJ gateway ...')
    print(f'Projection mode = {args.puncta_projection}')
    start_time = time.perf_counter()
    ij = imagej.init(
        ij_dir_or_version_or_endpoint='/Applications/Fiji.app',
        mode = "headless",
        add_legacy = True,
    )
    end_time = time.perf_counter()
    print(f'ImageJ initialized in {end_time - start_time:.3f} seconds')

    '''Find files to process'''
    file_list = []
    input_dir = Path(args.input_dir)
    for subdir in input_dir.iterdir(): #timepoints
        if subdir.is_dir(): 

            for file_path in subdir.iterdir():
                if file_path.is_file(): #each image inside timepoint folder

                    #only nd2 files and not hidden 
                    if file_path.suffix == '.nd2' and not str(file_path).startswith('.'):
                        file_list.append(file_path)
    
    '''
    #process
    last_timepoint = [file_path for file_path in file_list if ('180' in file_path.parent.name and 'after' not in file_path.name)]
    #run pipeline, save z-score of auto threshold
    threshold_z = []
    for file_path in last_timepoint:
        threshold_z.append(pipeline(str(file_path), str(file_path.parent), ij, args.save_intermediate, analyze_puncta = True))
        file_list.remove(file_path)
    threshold_z = float(np.mean(threshold_z)) #avg z-score cutoffs if multiple 180 timepoint pictures
    
    #run pipeline, use z-score to threshold
    print(f'Applying threshold z-score {threshold_z:.3f} to all images ...')
    for file_path in file_list:
        print(f'Timepoint: {file_path.parent.name}')
        pipeline(str(file_path), str(file_path.parent), ij, args.save_intermediate, analyze_puncta = True,
                 threshold = threshold_z)
    '''

    #Run pipeline on all files
    print(f'Analyzing {len(file_list)} images')
    print(f'LoG threshold {args.log_threshold}, Watershed threshold {args.watershed_threshold}')
    for file_path in tqdm(sorted(file_list, reverse=True)):
        print(f'Timepoint: {str(file_path.parent.name)}')
        pipeline(str(file_path),
                 str(file_path.parent),
                 ij,
                 args.save_intermediate,
                 analyze_puncta = True,
                 log_threshold = args.log_threshold,
                 watershed_threshold = args.watershed_threshold,
                 puncta_projection = args.puncta_projection,
                 crop = args.crop)
    
    #force shutdown to kill waiting workers which hang sometimes / why??
    #TODO: there is probably a cleaner solution to this
    sj.jimport('java.lang.System').exit(0)
    print(f'Processed {len(file_list) + 1} files')

def pipeline(
        input_file: str,
        output_dir: str,
        ij,
        save_intermediate: bool = False,
        analyze_puncta: bool = False,
        threshold: float = 0.0,
        log_threshold: float = 3.5,
        watershed_threshold: float = 3.0,
        puncta_projection = 'sum',
        crop = 100
):
    #Create folder to save results
    name = ''.join(os.path.basename(input_file).split('.')[:-1]) #remove ext
    result_dir = os.path.join(output_dir, f'{name}-processed')
    os.makedirs(result_dir, exist_ok = True)
    print(f'Processing {name}')

    #segment cells
    start_time = time.perf_counter()
    mask, puncta_image, num_cells = segment.cellpose(
        input = input_file,
        result_dir = result_dir,
        ij = ij,
        save_intermediate = save_intermediate,
        exclude_slices=3,
        name = name,
        projection = puncta_projection,
        crop = crop
    )
    end_time = time.perf_counter()
    print(f'Segmented {num_cells} cells in {end_time - start_time:.3f} seconds')

    threshold_z = threshold
    if analyze_puncta:
        #segment puncta
        start_time = time.perf_counter()
        '''
        Global method
        puncta_mask, puncta_num, threshold_z = segment.puncta(
            input = input_file,
            result_dir = result_dir,
            name = name,
            ij = ij,
            threshold = threshold
        ) # type: ignore
        '''
        
        #Local method
        puncta_mask, puncta_num = segment.puncta_local(
            result_dir = result_dir,
            name = name,
            cell_mask = mask,
            log_threshold = log_threshold,
            watershed_threshold = watershed_threshold,
            projection = puncta_projection
        )

        end_time = time.perf_counter()
        print(f'Segmented {puncta_num} puncta in {end_time - start_time:.3f} seconds')

    #analyze
    start_time = time.perf_counter()
    #cells
    cells_df = analysis.cells(
        cell_mask = mask,
        image = puncta_image,
        name = name,
        result_dir = result_dir,
    )
    if analyze_puncta:
        #puncta
        puncta_df = analysis.puncta(
            puncta_mask = puncta_mask, # type: ignore
            image = puncta_image,
            name = name,
            result_dir = result_dir
        )
        #correlate puncta to cells
        analysis.puncta_cells(
            puncta_mask = puncta_mask, # type: ignore
            cell_mask = mask,
            image = puncta_image,
            puncta_df = puncta_df,
            cells_df = cells_df,
            name = name,
            result_dir = result_dir
        )
    end_time = time.perf_counter()
    print(f'Analysis completed in {end_time - start_time:.3f} seconds')

    if threshold_z != 0.0:
        return threshold_z
    return 0.0

if __name__ == "__main__":
    main()