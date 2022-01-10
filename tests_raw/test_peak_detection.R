
## Load functions ####
source("./R/check_pretreatment.R")
source("./R/pretreatment_full.r")
source("./R/pretreatment_data.r")
source("./R/SeriesDetrend.r")
source("./R/Global_Thresh.r")
source("./R/Local_Thresh.r")
source("./R/SNI.r")
source("./R/PlotAnomalies.r")
source("./R/Plot_ReturnIntervals.r")
source("./R/peak_detection.R")


## Load Code Lake (Higuera et al. 2009): charcoal data ####
load("./Data-In/co_char_data.rda")


## Use Code Lake charcoal data with R-PaleoAnomalies ####

## With a 'global' GMM-inferred model
co_thresh_glob <- peak_detection(series = co_char_data, proxy = "char",
                                 thresh_type = "global", out_dir = "Figures")

Plot.Anomalies(co_thresh_glob, plot.neg = F)

Plot_ReturnIntervals(co_thresh_glob)

## With a 'local' GMM-inferred models
co_thresh_loc <- peak_detection(series = co_char_data, proxy = "char", out = "accI",
                                first = -51, last = 7500, yrInterp = 15,
                                detr_type = "mov.median", smoothing_yr = 500,
                                thresh_type = "local", thresh_value = 0.95,
                                min_CountP = 0.95, MinCountP_window = 150,
                                series_name = "Code Lake", sens = F,
                                out_dir = "Figures")

co_thresh_loc <- peak_detection(series = co_char_data, proxy = "char",
                                first = -51, last = 7500, yrInterp = 15,
                                detr_type = "mov.median", sens = F, out_dir = "Figures")
