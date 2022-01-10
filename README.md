# R-PaleoAnomalies
 
The set of functions gathered under the hood of ***R-PaleoAnomalies*** is meant to be used for analyzing paleoecological records, when the goal is peak detection to reconstruct the occurrence, the return intervals, and the magnitude of distinct events.

  
  
## Rationale for developing *R-PaleoAnomalies*
*R-PaleoAnomalies* builds on *CharAnalysis* (https://github.com/phiguera/CharAnalysis), a software for analyzing sediment-charcoal records written in and compiled with Matlab 7.0 by Phil Higuera (Higuera et al., 2009), with significant input by (amongst others) Patrick Bartlein (U of OR), Daniel Gavin (U of OR), Jennifer Marlon, and Ryan Kelly.

Two main reasons led to the development of *R-PaleoAnomalies*. Firstly, as R is an open source product, modifying the program to suit individual needs may be more straigthforward. Secondly, an integration and inter-operability with other existing R-packages may allow using peak-detection analysis in conjunction with other workflows and types of paleoecological records (see for instance [Cagliero et al., 2021](https://doi.org/10.1007/s00334-021-00862-x)).

***
  
## Usage
A typical workflow of the peak-detection analysis as implemented in *CharAnalysis* includes the following steps (Higuera et al., 2011):
* 1.) *resample* the record to equally spaced sampling intervals in time (years);
* 2.) *decompose* the resampled record into a long-term trend (background component) and peaks (peak component);
* 3.) *screen* the peak component to distinguish signal from noise using: 
  + 3.1.1) a unique *global* 2-component Gaussian mixture model, or
  + 3.1.2) *local* 2-component Gaussian mixture models,
  + 3.2) and eventually also a minimum-count test;
* 4.) *evaluate* the suitability of the record for peak-detection analysis using the signal-to-noise index (Kelly et al., 2011).

*R-PaleoAnomalies* performs steps 1.) and 2.) for several variables of one dataset type (e.g. different estimates of charcoal abundance). Steps 3.) and 4.) are performed for one user-selected variable.

To run your own data, make a new folder within an umbrella folder, and save it under a name, e.g., `Data-In`. Then place a file (e.g., `MyData.csv`) in that folder. The input data file will contain the sample depths, sample ages, sample volume, and the variable(s). The file should have the following formatting: It has headers and at least six fields. The first five columns will report the metadata for the samples, the subsequent columns contain the variable(s) to be analysed (e.g., the abundance of charcoal pieces).

CmTop | CmBot | AgeTop | AgeBot | Volume | variable1 | variable2 | ... | nth-variable
------|-------|--------|--------|--------|-----------|-----------|-----|-------------
 0.5  |  1    | -42    | -24    | 3      | 8         | 0.01      | ... |    ...      
 1    |  1.5  | -24    | 5      | 3      | 18        | 0.005     | ... |    ...      
 

The depths and ages should be arranged in ascending order. Sample ages should thus be reported as *calendar ages BP*.


Load the data into the R environment
> MyData <- read.csv("./Data-In/MyData.csv")

Until the packaging is not finished, download the entire *R-PaleoAnomalies* program as a .zip or tar.gz archive [here](https://github.com/wfinsinger/R-PaleoAnomalies/archive/refs/heads/main.zip). Alternativly, download individual files by visiting the GitHub pages.

The `check_pretreat()` function can be used verify the input data is formatted correctly. If the samples in the input file are not contiguous, the `check_pretreat()` function will add rows for the missing samples. Should the dataset contain samples that were deposited in a very short amount of time (e.g., slumps, tephras), for which AgeTop = AgeBot, these samples will be flagged and removed, and a new depth scale will be created to replace the original one.
> MyData <- check_pretreat(Mydata)

The functions can either be run individually and stepwise, or the wrapper function `peak_detection()` can be used to perform an analysis including steps from 1.) to 4.) in one go.

For instance, to analyse the MyData dataset for variable1:
> MyData_peaks <- peak_detection(series = MyData, proxy = "variable1")

***

## Example
*R-PaleoAnomalies* comes with an example dataset, called `co_char_data` ([Code Lake](https://figshare.com/articles/dataset/Higuera_et_al_2009_lake_sediment_pollen_and_charcoal_data/984310/4), Higuera et al., 2009). It is placed in the folder 'Data-In'.

The dataset can be loaded into the R environment
> load("./Data-In/co_char_data.rda")

and analysed, for instance, with the following settings (leaving other arguments with their default values) 
> co_thresh_loc <- peak_detection(series = co_char_data, proxy = "char",
                                first = -51, last = 7500, yrInterp = 15,
                                detr_type = "mov.median", sens = F)

With these settings, the results obtained using *R-PaleoAnomalies* strikingly resemble those obtained with *CharAnalysis*.

![Code Lake: peak-detection outputs](/README_Figures/01_Code_Lake_peak_detection.jpg "Code Lake: peak-detection outputs")

![Code Lake: reconstructed fire-return intervals (FRI)](/README_Figures/02_Code_Lake_FRIs.jpg "Code Lake: reconstructed fire-return intervals (FRI)")


***

## Acknowledgements & Credits
The development of this set of functions would not have been possible without the Matlab-coded templates that were written and made open source by Philip Higuera (https://github.com/phiguera/CharAnalysis). I'm thankful to Dan Gavin for suggesting tweaks to accomodate for non-standard data-input formats. His suggestions led to the `check_pretreatment()` function. I'm extremely thankful to Petr Kunes and several CharAnalysis users who ignited the conversation over the past few years about getting that program into R.

If you use *R-PaleoAnomalies* in your publications, please cite *https://github.com/wfinsinger/R-PaleoAnomalies* and any non-default settings applied. The packaging of these functions is in progress.

***

## Issues & Contributions
If you are having problems running *R-PaleoAnomalies* or if you have any suggestions, please use the "Issues" tab.
Contributions to this work are more than welcome. You can just fork, make changes,and then file a pull request. Alternatively, you can get in touch with me to discuss how your improvements may fit with ongoing development of add-ons. Thanks!


## References
> Cagliero E, Morresi D, Paradis L, Curović M, Spalević V, Marchi N, Meloni F, Bentaleb I, Motta R, Garbarino M, Lingua E, Finsinger W (2022) Legacies of past human activities on one of the largest old-growth forests in south-east European mountains. *Vegetation History and Archaeobotany*, online available [link](https://doi.org/10.1007/s00334-021-00862-x).


> Higuera PE, LB Brubaker, PM Anderson, FS Hu, TA Brown (2009) Vegetation mediated the impacts of postglacial climatic change on fire regimes in the south-central Brooks Range, *Alaska Ecological Monographs* 79: 201-219 [link](https://doi.org/10.1890/07-2019.1)

> Higuera PE, Gavin DG, Bartlein PJ, Hallett DJ (2010) Peak detection in sediment–charcoal records: impacts of alternative data analysis methods on fire-history interpretations. *International Journal of Wildland Fire* 19: 996. [link](http://dx.doi.org/10.1071/WF09134)

> Kelly RF, PE Higuera, CM Barrett, FS Hu (2011) A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records *Quaternary Research* 75: 11-17 [link](http://dx.doi.org/10.1016/j.yqres.2010.07.011)