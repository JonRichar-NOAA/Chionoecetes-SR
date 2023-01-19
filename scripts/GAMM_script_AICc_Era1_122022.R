############ This script runs GAM models for time period of 1978 to 2005 ###########

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

dat0 <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed
names(dat0)

dat0 %>%
  subset(Era_AICc==1)->dat1

dat1
dat1 %>%
  select(releaseyear:NE.wind) ->dat

head(dat)

dat



## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier
mod1 <- gamm(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat[dat$releaseyear >= 1981,],correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod1)
summary(mod1$gam)
summary(mod1$lme)

#Inspect model object
mod1

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod1$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod1$lme, resid=T, pch=19, rug=F, se=F, pages=1)

gam.check(mod1)

## evaluate models including additional covariates in S-R model ---------------------------
cor(dat$ReproductiveFemales, dat$Ovig_female_CO, use="p") # not correlated

mod2 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod2)
summary(mod2$gam)
summary(mod2$lme)

#Inspect model object
mod2

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod2$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod2)

# 
# 

## cod
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.22
dat$PCod_RA3
dat$ReproductiveFemales
mod3 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(PCod_RA3, k=4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod3)
summary(mod3$gam)
summary(mod3$lme)

#Inspect model object
mod3

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod3$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod3$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod3) # almost the same

## FHS lag 2

cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.17

mod4 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod4)
summary(mod4$gam)
summary(mod4$lme)

#Inspect model object
mod4

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod4$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod4$lme, resid=T, pch=19, rug=F, se=F, pages=1)


## NE wind
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$NE.wind, use="p") # -0.29

mod5 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod5)
summary(mod5$gam)
summary(mod5$lme)

#Inspect model object
mod5

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod5$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod5$lme, resid=T, pch=19, rug=F, se=F, pages=1)


AICc(mod4, mod5) # adding wind not helpful, even makes model worse!



## SE wind
cor(dat$FHS_lag2, dat$SE.wind, use="p") # 0.11
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # 0.17

mod6 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SE.wind, k=4),
             data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod6)
summary(mod6$gam)
summary(mod6$lme)

#Inspect model object
mod6

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod6$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod6$lme, resid=T, pch=19, rug=F, se=F, pages=1)

AICc(mod5, mod6) # adding wind not helpful!


## AO
cor(dat$FHS_lag2, dat$AO_RA3, use="p") #0.67.....highly correlated!
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # 0.55...ALSO highly correlated

mod7 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod7)
summary(mod7$gam)
summary(mod7$lme)

#Inspect model object
mod7

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod7$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod7$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5, mod7) # adding AO not helpful!


## and bottom temp rolling average
plot(FHS_lag2~NBT_3RA,data = dat[dat$releaseyear >= 1981,],pch = 16)

cor(dat$FHS_lag2, dat$NBT_3RA, use="p") # -0.45.....correlated!
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p") #-0.23....weakly correlated
cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") #-0.23....weakly correlated

mod8 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod8)
summary(mod8$gam)
summary(mod8$lme)

#Inspect model object
mod8

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod8$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod8$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5, mod8) # adding bottom temp only marginally better


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")            #-0.38...correlated
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p") #0.16

mod9 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
            data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod9)
summary(mod9$gam)
summary(mod9$lme)

#Inspect model object
mod9

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod9$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod9$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5, mod9)


## refit with linear FHS effect
## No need for correlation calculations as repeat of prior model with modded term

mod10 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
             data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod10)
summary(mod10$gam)
summary(mod10$lme)

#Inspect model object
mod10

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod10$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod10$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5, mod10) # not surprisingly - not any better


## PDO
mod11 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PDO_RA3, k=4),
             data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod11)
summary(mod11$gam)
summary(mod11$lme)

#Inspect model object
mod11

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod11$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod11$lme, resid=T, pch=19, rug=F, se=F, pages=1)


MuMIn::AICc(mod5, mod11) # quite an improvement to include PDO

## Add Pacific cod combined with FHS
cor(dat$FHS_lag2, dat$PCod_RA3, use="p")            #-0.38...correlated
cor(dat$ReproductiveFemales, dat$PCod_RA3, use="p") #0.16
mod12 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PCod_RA3, k=4),
              data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod12)
summary(mod12$gam)
summary(mod12$lme)

#Inspect model object
mod12

#plot
dev.new()
par(mfrow=c(3,1))

plot(mod12$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod12$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5, mod12)

## refit with linear FHS effect ONLY
## No need for correlation calculations as repeat of prior model with modded term

mod13 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod13)
summary(mod13$gam)
summary(mod13$lme)

#Inspect model object
mod13

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod13$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod13$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod13) # not surprisingly - not any better

###################### Wrap up ##########################################

MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12,mod13) ## Model 7 is best performing for Era 1 by AICc

names(dat)
cor.dat <- dat %>%
  select(-era,-releaseyear)

c<-cor(cor.dat, use="p")
c
write.csv(c,"output/Variable_correlations_AICc_Era1.csv")
