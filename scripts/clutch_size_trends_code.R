library(stats)
library(MASS)
library(nlme)
library(lmtest)
library(mgcv)
library(nlme)
library(ncdf4)
library(chron)
library(lattice)
library(nlstools)
library(MuMIn)
library(tidyverse)
library(corrplot)
library(voxel)
#?nlstools

########## Import data 
getwd()

dat <- read.csv("data/CB_FEMALE_CLUTCH_SIZE_POP.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)
names(dat)

plot(dat$TOTAL_FEMALE_MATURE_POP~dat$SURVEY_YEAR,type="l")
lines(dat$TOTAL_FEMALE_BARREN_POP~dat$SURVEY_YEAR)

pct_barren<-(dat$TOTAL_FEMALE_BARREN_POP/dat$TOTAL_FEMALE_MATURE_POP)*100

pct_barren
mean(pct_barren)

mean_dat<-as.data.frame(cbind(dat$SURVEY_YEAR,pct_barren))
colnames(mean_dat)<-c("Year","pct_barren")
plot(mean_dat$pct_barren~mean_dat$Year,type="l")
