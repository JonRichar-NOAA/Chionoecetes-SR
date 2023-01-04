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

## Create Log(R) dataset
logR<-log(dat$Lag3_recruits)

## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier
mod1 <- gam(logR ~ s(ReproductiveFemales, k = 4),
            data = dat)

summary(mod1)
plot(mod1, resid=T, pch=19, rug=F, se=F)

dat$era <- as.factor(if_else(dat$releaseyear <= 1995, 1, 2))

# fit mod1 without years for which FHS_RA2 is unavailable

mod2<- gam(logR ~ s(ReproductiveFemales, k = 4),
            data = dat[dat$releaseyear >= 1983,])

summary(mod2)
plot(mod2, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod1, mod2) 



summary(mod2)
plot(mod2, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod2) ## mod 2 better by ~ 2 units

## examine the residuals for mod1
resid1 <- data.frame(year=dat$releaseyear,
                     resid=residuals(mod1))

ggplot(resid1, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod1)~1)
acf(resid(mod1)) # not perfect, not terrible!

gam.check(mod1)

## evaluate models including additional covariates in S-R model ---------------------------
cor(dat$ReproductiveFemales, dat$Ovig_female_CO, use="p") # not correlated

mod3 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat)

summary(mod3)
MuMIn::AICc(mod1, mod3)

plot(mod3, resid=T, pch=1, se=F, rug=F, pages=1)
# 
# 

## cod
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.22

mod4 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat)

summary(mod4)
MuMIn::AICc(mod1, mod4) # almost the same

## FHS lag 2

cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.17

mod5 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod5)
plot(mod5, resid=T, pch=1, se=F, rug=F, pages=1)

## examine the residuals for mod5

resid5 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                     resid=residuals(mod5))

ggplot(resid5, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod5)~1)
acf(resid(mod5)) # again, not great-great

## FHS rolling average

cor(dat$ReproductiveFemales, dat$FHS_RA2, use="p") # 0.20

mod6 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_RA2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod6)
plot(mod6, resid=T, pch=1, se=F, rug=F, pages=1)


## retain FHS, and see if it's possible to improve with AO/wind/other physical covariates?? ----------
## Model 7 is best performing for Era 1 by AICc

cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.17

mod7 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat)

summary(mod7)

plot(mod7, resid=T, pch=1, se=F, rug=F, pages=1) # FHS effect is negative/linear = plausible!

## NE wind
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$NE.wind, use="p") # -0.29

mod8 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat)

summary(mod8)
plot(mod8, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod6, mod8) # adding wind not helpful, even makes model worse!
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated

## examine the residuals for mod15
resid8 <- data.frame(year=dat$releaseyear[4:17],
                      resid=residuals(mod8, type="response"))

ggplot(resid8, aes(year, resid)) +
  geom_line() +
  geom_point() +
  geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod8)~1)
acf(resid(mod8)) #pretty good!

## SE wind
cor(dat$FHS_lag2, dat$SE.wind, use="p") # 0.11
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # 0.17

mod9 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SE.wind, k=4),
             data = dat)

summary(mod9)
plot(mod9, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod6, mod9) # adding wind not helpful!


## AO
cor(dat$FHS_lag2, dat$AO_RA3, use="p") #0.67.....highly correlated!
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # 0.55...ALSO highly correlated

mod10 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
            data = dat)

summary(mod10)
plot(mod10, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod6, mod10) # adding AO not helpful!


## and bottom temp rolling average
cor(dat$FHS_lag2, dat$NBT_3RA, use="p") # -0.45.....correlated!
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p") #-0.23....weakly correlated

mod11 <- gam(logR ~ s(Reproductiv-eFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
            data = dat)

summary(mod11)

plot(mod11, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod6, mod11) # adding bottom temp only marginally better


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")            #-0.38...correlated
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p") #0.16

mod12 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
            data = dat)

summary(mod12)

plot(mod12, resid=T, pch=1, se=F, rug=F, pages=1) 

## shifts FHS into improbable shape!

## refit with linear FHS effect
## No need for correlation calculations as repeat of prior model with modded term
mod13 <- gam(logR ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
             data = dat)

summary(mod13)

plot(mod13, resid=T, pch=1, se=F, rug=F, pages=1) 
MuMIn::AICc(mod6, mod13) # not surprisingly - not any better


## PDO
mod14 <- gam(logR ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PDO_RA3, k=4),
             data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod14)

plot(mod14, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod6, mod14) # quite an improvement to include PDO

## examine the residuals for mod15
resid14 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod15, type="response"))

ggplot(resid14, aes(year, resid14)) +
  geom_line() +
  geom_point() +
  geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod14)~1)
acf(resid(mod14)) #pretty good!

MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14) ## Model 7 is best performing for Era 1 by AICc

names(dat)
cor.dat <- dat %>%
  select(-era,-releaseyear)

c<-cor(cor.dat, use="p")
c
write.csv(c,"output/Variable_correlations_AICc_Era1.csv")
