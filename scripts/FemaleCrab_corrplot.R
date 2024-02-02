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

#======================== Import data ================================================================
getwd()

dat <- read.csv("data/Female_Tanner_Crab_Series_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)
cor(dat, use="complete.obs")
cor2 <- cor(dat, use="complete.obs")
dev.new()
par(mfrow=c(1,1))
corrplot(cor2)
plot(dat$EBS_SC3~dat$EBS_SC2,main = "Same yr SC2 estimates")
#======================================================================================================

dat <- read.csv("data/Female_Tanner_Crab_Series_for_analysis_laggedSC2s.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)
cor(dat, use="complete.obs")
cor2 <- cor(dat, use="complete.obs")
dev.new()
corrplot(cor2)
plot(dat$EBS_SC3~dat$EBS_SC2, main = "y-1 SC2 estimates")
#======================================================================================================

