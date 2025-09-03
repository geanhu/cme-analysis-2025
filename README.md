# CME Analysis 2025

Code for analysis of confocal microscope images of CME events in yeast.

Rewritten from [this version](https://github.com/geanhu/cme-movie-analysis-2024), which inconveniently uses ImageJ to launch Groovy scripts and bash scripts to launch Python processes. Current version aims to use PyImageJ to enable better program portability and odd file I/O and global variable errors when using ImageJ to launch external processes.