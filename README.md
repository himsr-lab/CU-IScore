![header](https://user-images.githubusercontent.com/19319377/116955473-e20f9d00-ac4f-11eb-91fc-56399caedeb4.png)
# CU-IScore
## IHC scoring of multi-channel images
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4599591.svg)](https://doi.org/10.5281/zenodo.4599591)
### Scoring methods
The CU-IScore macro implements two scoring methods applicable for immunofluorescence or DAB/hematoxylin multi-channel images: In short, the novel method returns higher score values for immunofluorescence images with higher average pixel values, while the classic method returns higher score values for DAB/hematoxylin images with lower average pixel values. The classic method is based on the [IHC Profiler](https://doi.org/10.1371/journal.pone.0096801) design, but extends its applicability to multi-channel images with various image formats, bit depths, or pixel value ranges. In addition, CU-IScore creates detailed scoring reports for individual channels as well as a summary table for image batches.

![histo](https://user-images.githubusercontent.com/19319377/116958138-d4114a80-ac56-11eb-896b-89e4d8bb0a12.png)

**Figure 1: Histogram example.** Distribution of pixel counts as a function of pixel values from the grayscale image (top). The histogram is partitioned into four equally-sized intervals for the I-Score calculation. Pixels with values below a certain threshold (black triangle, bottom) are disregarded as background (white peak, left). Pixels within the user-specified value range contribute to the I-Score weighted by three pixel intensity categories (black distribution, center). In addition, pixels with values above a certain threshold (white triangle, bottom) can be disregarded as outliers (white line, right).

### Software documentation
The documentation of our macros is in the corresponding source code: You can view the source code on GitHub by following the links to the macros.

### Software requirements
CU-IScore requires a recent version of the [Fiji](https://fiji.sc/) image processing package:
* ImageJ2 app (>= 1.52a)
* Bio-Formats plugin (>= 6.4.0)

Any multi-channel image that can be imported with the Bio-Formats plugin can be processed by CU-IScore. However, you might have to adjust the [`suffixes`](https://github.com/christianrickert/CU-IScore/blob/44a05ef2cfef58cbf6988ee03a8dbb64a2206076/CU-IScore.ijm#L81) variable to select the file extensions for your specific instrument. In addition, if the metadata extraction and therefore the slice labeling fails, you will have to identify individual channels by slice number.

CU-IScore also requires a recent version of [CU-MacroLibrary](https://github.com/christianrickert/CU-MacroLibrary/) to be installed.

### Example files
The [example folder](https://github.com/christianrickert/CU-IScore/tree/main/example) contains a single [Vectra® Polaris™ image](https://github.com/christianrickert/CU-IScore/blob/main/example/Polaris%20Pt%2012%20Point%2013.tif?raw=true) (1176x1080 px).
Running CU-IScore with the default [`Variables`](https://github.com/christianrickert/CU-IScore/blob/35f04fcf80ba537980315bac7216f839f54cc220/CU-IScore.ijm#L77), should yield results identical to the data in the [results subfolder](https://github.com/christianrickert/CU-IScore/tree/main/example/Polaris%20Pt%2012%20Point%2013) as well as in the [summary table](https://github.com/christianrickert/CU-IScore/blob/main/example/CU-IScore.csv).

CU-IScore produces distinct result files for every multi-channel image in the batch:
* `*.csv` - [summary table](https://github.com/christianrickert/CU-IScore/blob/main/example/CU-IScore.csv)
* `*.txt` - [log file](https://github.com/christianrickert/CU-IScore/blob/main/example/Polaris%20Pt%2012%20Point%2013/CU-IScore.txt)

### Copyright notices
The [Vectra® Polaris™ image](https://www.akoyabio.com/phenoptics/mantra-vectra-instruments/vectra-polaris/) was kindly provided by [Benjamin Bitler](https://medschool.cuanschutz.edu/ob-gyn/divisions/division-of-reproductive-sciences/our-faculty/benjamin-bitler-phd) for demonstration purposes.
