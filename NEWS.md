# tapas 0.1.5

- check_pretreat() function: added some lines to flag for the presence of 
  overlapping sample depths. Modified some other warnings as well.
  
- added the get_overlap_depths() function.

# tapas 0.1.4

- added a function to load results of charcoal analyses made with ImageJ


# tapas 0.1.3

- modified the pretreatment() function to calculate sediment-accumulation
  rates in a different way: now the code uses explicit sample thicknesses and
  sample-deposition times rather than diff(depth) and diff(age). Many thanks
  to Pat Bartlein for flagging this issue (26/02/2024).

- Added the pretreatment() function from the paleofire v1.2.4 R package
  (30/10/2023).

- Corrected the calculation of peak magnitude (06/06/2023).

- Added sample data sets to run examples for arco() function.

- Added the arco() function to screen charcoal-area records.

- Added two sample data sets (rand_peaks) and (red_noise) for users to
  play around.

- Added the cpts_ar() function that determines zone boundaries for
  single proxy accumulation-rate records based on a change-point
  analysis, and checks for the influence of sediment-accumulation rates
  on the detected change points.

# tapas 0.1.2

- Added the plot_raw() function to produce a polygon plot and optionally
  overlay that with a bar plot of the raw input data.

- Added two functions (tapas2mgcv and mgcv2tapas) that allow determine
  the ‘background trend’ with a GAM-fitted model. This allows to by-pass
  the SeriesDetrend() function.

# tapas 0.1.1

- Slightly modified the local_thresh() function. Now it also
  accommodates discontinuous and/or non-binned data sets. Though, this
  is currently at a work-in-progress stage, as the SeriesDetrend()
  function still requires continuous and binned temporal series (i.e. an
  output from the pretreatment_data() function).

- Got rid of trailing horizontal white spaces.

- Fixed issue in pretreatment_data() arising when only one variable is
  included.

- Fixed issue with the SNI as it did not call the correct data.

- Added a new function export_tapas() that gathers in one data.frame the
  output of the peak_detection(), local_thresh(), or global_thresh()
  functions. Many thanks to Kristen Beck for suggesting this.

- The selection of values to be evaluated with the local GMMs is now
  simply based on the “smoothing.yr” parameter. Previously, the
  selection was based on their index, which was more tedious.

- Fixed issue with “Plot.Anomalies()”, as y-axis labels were not printed
  correctly.

# tapas 0.1.0

- Initial release of the package
