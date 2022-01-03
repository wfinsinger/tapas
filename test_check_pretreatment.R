

## Load data ####
co <- read.csv("./Data-In/CO_short.csv")
co$chararea <- co$count * 0.01542 # just add a dummy variable to see if everything works

## Load source files ####
source("./R/check_pretreatment.r")
source("./R/pretreatment_full.r")


### Test with the paleofire::pretreatment() function
pippo <- paleofire::pretreatment(params = co[ ,1:5], serie = co[ ,6])
plot.CHAR(pippo)

### Test with the pretreatment_full() function ####
pippo <- pretreatment_full(params = co[ ,1:5], serie = co[ ,6])
plot.CHAR(pippo)



### Test when gaps are introduced in data set ####
co_gaps <- co[-c(2:3, 12), ]

# The pretreatment() function fails because the sample depths and the sample ages are not continuous
pippo2 <- paleofire::pretreatment(params = co_gaps[ ,1:5], serie = co_gaps[ ,6])
plot.CHAR(pippo2)

# Add missing samples, and use pretreatment 
co_gaps <- check_pretreat(series = co_gaps)
pippo2 <- paleofire::pretreatment(params = co_gaps[ ,1:5], serie = co_gaps[ ,6], Int = T)
plot.CHAR(pippo2)



### Test with a 'slump' sample (AgeTop[i] == AgeBot[i])
co_t <- read.csv("./Data-In/CO_short_tephra.csv")
co_t <- check_pretreat(series = co_t)
pippo3 <- paleofire::pretreatment(params = co_t[ ,1:5], serie = co_t[ ,6])
plot.CHAR(pippo3)
plot.CHAR(pippo)

