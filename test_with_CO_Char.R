
##########################################################################################
# R code to detect positive and negative anomalies in palaeoecological time series.
# Walter Finsinger, June 2020 - October 2021

# The functions are loaded at the beginning. Explanations concerning the
#   user-determined parameters are given in the functions files.

# The code does the following:
# 
# 1) the time series are virtually resampled to equal sampling intervals.
# To do this, I used the pretreatment() function from the paleofire() package
# (aka pretreatment as from CharAnalysis, Higuera et al., 2009), which I slightly
# modified in order to get all output data needed later.
# 
# 2) the resampled time series are detrended using a smoothing function, 
# and Figures are produced to visually evaluate the smoothed record and the detrended data.
# 
# 3) the code fits a 2-component Gaussian Mixture Model (GMM) to the
# detrended data. One of the components is centred around 0, and we set
# the nth quantile of the gaussian distribution of that component to determine
# the threshold separating the noise from the signal. Figures are produced to visually
# evaluate the GMMs.
# 
# 4) the code finally allows plotting the resampled data, the smoothed data,
# the thresholds, as well as the detected anomalies (values beyond the noise threshold).
##########################################################################################


rm(list = ls())


### 1) Read r functions ########################################################
source("./R/pretreatment_full.r")
source("./R/pretreatment_data.r")
source("./R/SeriesDetrend.r")
source("./R/Global_Thresh.r")
source("./R/Local_Thresh.r")

source("./R/SNI.r")
source("./R/PlotAnomalies.r")
source("./R/Plot_ReturnIntervals.r")



### Load data ############################################################################

### Charcoal data ###
co <- read.csv("./Data-In/CO_charData.csv", header = T)
co <- dplyr::rename(co, "AgeTop" = "ageTop..yr.BP.")
co <- dplyr::rename(co, "AgeBot" = "AgeBot..yr.BP.")
co <- dplyr::rename(co, "char" = "CharCount.....")
co$dummy_chararea <- co$char * 0.2421

## Resample data #########################################################################
min.ages <- min(co$AgeTop)
max.ages <- max(co$AgeTop)
yr.interp <- 15
co_i <- pretreatment_data(series = co, out = "accI", series.name = "co",
                             first = -50, last = 7450, yrInterp = yr.interp)
plot(co_i$raw$series$age, co_i$raw$series$char, type = "l")
lines(co_i$int$series.int$age, co_i$int$series.int$char, col = "red")


## Detrend the data ######################################################################
co_detr <- SeriesDetrend(series = co_i, smoothing.yr = 500,
                            detr.type = "rob.loess", out.dir = "Figures")
co_detr <- SeriesDetrend(series = co_i, smoothing.yr = 500,
                          detr.type = "mov.mode", out.dir = "Figures")
#rm(char_thresh_gl, co, co_i)

save.image(file = "./Data-Out/co_detr.rdata")



## Fit gaussian mixture models to determine the noise-signal threshold ###################

## With an overall (global threshold) & with minimum-count test
char_thresh_gl <- global_thresh(series = co_detr, proxy = "charAR",
                                thresh.value = 0.95, smoothing.yr = 500,
                                keep_consecutive = F, minCountP = 0.05,
                                out.dir = "Figures")

char_thresh_gl <- global_thresh(series = co_detr, proxy = "charAR",
                              thresh.value = 0.95, smoothing.yr = 500,
                              keep_consecutive = F, minCountP = 0.05,
                              out.dir = "Figures")

char_thresh_gl <- global_thresh(series = co_detr, proxy = "dummy_charareaAR",
                              thresh.value = 0.95, smoothing.yr = 500,
                              keep_consecutive = F, minCountP = NULL,
                              out.dir = "Figures")

## With overall threshold & w/out minimum-count test
char_thresh_gl <- global_thresh(series = co_detr, proxy = "charAR",
                               thresh.value = 0.95, smoothing.yr = 500,
                               keep_consecutive = F, minCountP = 0.05,
                               out.dir = "Figures")

char_thresh_gl <- global_thresh(series = co_detr, proxy = "charAR",
                              thresh.value = 0.95, smoothing.yr = 500,
                              keep_consecutive = F, minCountP = 0.05,
                              out.dir = "Figures")

## With locally-fitted thresholds
char_thresh_loc <- local_thresh(series = co_detr, proxy = "charAR", thresh.yr = 500,
                                thresh.value = 0.95, smoothing.yr = 500,
                                keep_consecutive = F, minCountP = 0.05,
                                out.dir = "Figures")



## Make final plots ######################################################################
layout(1)
par(mfrow = c(2,1), mar = c(0,5,2,0.5), oma = c(4,1,0,0))
Plot.Anomalies(series = char_thresh_gl, plot.crosses = T,
               plot.x = F, plot.neg = F)
Plot_ReturnIntervals(series = char_thresh_gl, plot.x = T, plot.neg = F)

Plot.Anomalies(series = char_thresh_loc, plot.crosses = T,
               plot.x = F, plot.neg = F)
Plot_ReturnIntervals(series = char_thresh_loc, plot.x = T, plot.neg = F)


