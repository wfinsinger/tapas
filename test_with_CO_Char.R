
##########################################################################################
# R code to detect positive and negative anomalies in palaeoecological time series.
# Walter Finsinger, June 2020 - )ctober 2021

# The functions are loaded at the beginning. Explanations concerning the
#   user-determined parameters are given in the functions files.

# The code does the following:
# 
# 1) the time series are virtually resampled to equal sampling intervals.
# To do this, I used the pretreatment() function from the paleofire() package
# (aka pretreatment as from CharAnalysis, Higuera et al., 200x), which I slightly
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
source("./R/SeriesThreshold_v5.r")
source("./R/SNI.r")
source("./R/PlotAnomalies.r")
source("./R/Plot_ReturnIntervals.r")



### Load data ##########################################################

### Charcoal data ###
co <- read.csv("./Data-In/CO_charData.csv", header = T)
co <- dplyr::rename(co, "AgeTop" = "ageTop..yr.BP.")
co <- dplyr::rename(co, "AgeBot" = "AgeBot..yr.BP.")
co <- dplyr::rename(co, "char" = "CharCount.....")

## Pretreatment
min.ages <- min(co$AgeTop)
max.ages <- max(co$AgeTop)
yr.interp <- 15
co.i1 <- pretreatment_data(series = co, out = "accI", series.name = "co",
                             first = -51, last = 7500, yrInterp = yr.interp)
plot(co.i1$raw$series$age, co.i1$raw$series$char, type = "l")
lines(co.i1$int$series.int$age, co.i1$int$series.int$char, col = "red")

## Detrending
co.detr1 <- SeriesDetrend(series = co.i1, smoothing.yr = 500,
                            detr.type = "rob.loess", out.dir = "Figures")
co.detr1 <- SeriesDetrend(series = co.i1, smoothing.yr = 500,
                          detr.type = "mov.mode", out.dir = "Figures")
#rm(char.thresh.gl1, co, co.i1)

## Fit gaussian mixture models to determine the noise-signal threshold
char.thresh.gl1 <- SeriesThresh(series = co.detr1, proxy = "charAR", thresh.yr = 500,
                                thresh.value = 0.99, smoothing.yr = 500, span.sm = 0.2,
                                gm.local = F, keep_consecutive = F,
                                out.dir = "Figures")


char.thresh.loc <- SeriesThresh(series = co.detr1, proxy = "charAR", thresh.yr = 500,
                                thresh.value = 0.99, smoothing.yr = 500,
                                gm.local = T, keep_consecutive = F,
                                out.dir = "Figures")



## Make final plots
layout(1)
par(mfrow = c(2,1), mar = c(1,5,2,0.5))
Plot.Anomalies(series = char.thresh.gl1, proxy = "charAR", plot.crosses = T,
               plot.x = F, plot.neg = F)
Plot_ReturnIntervals(series = char.thresh.gl1, proxy = "charAR", plot.neg = F)

Plot.Anomalies(series = char.thresh.loc, proxy = "charAR", plot.crosses = T,
               plot.x = T, plot.neg = F)
Plot_ReturnIntervals(series = char.thresh.loc, proxy = "charAR", plot.neg = F)


