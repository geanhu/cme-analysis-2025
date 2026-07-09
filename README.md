# CME Analysis 2025

Code for analysis of confocal microscope images of CME events in yeast.

Rewritten from [this version](https://github.com/geanhu/cme-movie-analysis-2024), which inconveniently uses ImageJ to launch Groovy scripts and bash scripts to launch Python processes. Current version aims to use PyImageJ to enable better program portability and odd file I/O and global variable errors when using ImageJ to launch external processes.

## Installation
Create a conda environment with the correct dependencies using the `.yml` file [here](/install/environment.yml)
```
conda env create -f /install/environment.yml
```
Activate the conda environment before running any code
```
conda activate cme_analysis_2025
```
Please also install the necessary ImageJ plugins to your local installation of ImageJ as listed [here](install/fiji-packages.txt). Packages can be installed through the ImageJ GUI through `Help > Update > Manage` update sites.
