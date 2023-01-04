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
#?nlstools

########## Import data 
getwd()

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)


########## Run models
mod1 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3)+ s(PDO_RA3, k=4)+ s(AO_RA3, k=4),
            data = dat)
summary(mod1)
plot(mod1, resid=T, pch=19, rug=F, se=F)

#Now run with era interaction
mod2 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3)+ s(PDO_RA3, k=4, by =Era_AICc)+ s(AO_RA3, k=4, by = Era_AICc),
            data = dat)
summary(mod2)
plot(mod1, resid=T, pch=19, rug=F, se=F)

MuMIn::AICc(mod1, mod2) 
