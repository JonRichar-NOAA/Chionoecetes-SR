library(stats)
library(MASS)
library(nlme)
library(lmtest)
library(mgcv)
library(ncdf4)
library(chron)
library(lattice)
library(nlstools)
library(MuMIn)
library(tidyverse)
library(corrplot)
#?nlstools
########## package citation
citation("MuMIn")
########## Import data 

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

# load NE wind
wind <- read.csv("./data/mean.tanner.wind.csv", row.names = 1)
names(wind) <- c("releaseyear", "NE.wind")
dat
dat <- left_join(dat, wind)

hist(dat$logRS)

cor.dat <- dat %>%
  select(-era)

cor(cor.dat, use="p")


## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier
?gam()
?gls()
mod1 <- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat)

summary(mod1)
plot(mod1, resid=T, pch=19, rug=F, se=F)

# fit mod1 without years for which FHS_RA2 is unavailable
mod1b<- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat[dat$releaseyear >= 1983,])

summary(mod1b)


##break down by era
dat$era <- as.factor(if_else(dat$releaseyear <= 2000, 1, 2))

mod2 <- gam(logRS ~ s(ReproductiveFemales, k = 4, by = era),
            data = dat)

summary(mod2)
plot(mod2, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod2) ## mod 2 better by ~ 2 units

## plot the raw data by era

ggplot(dat, aes(ReproductiveFemales, logRS)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4)) +
  facet_wrap(~era)

## ok, so there is no basis for splitting by era
## there are simply fewer females in era 2, so the
## full range of the relationship is not present for analysis!
## should make our work easier!

## examine the residuals for mod1
resid1 <- data.frame(year=dat$releaseyear,
                     resid=residuals(mod1))

ggplot(resid1, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod1)~1)
acf(resid(mod1)) # not perfect, not terrible!

gam.check(mod1)

## examine the residuals for mod1b
resid1b <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod1b))

ggplot(resid1b, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod1b)~1)
acf(resid(mod1b)) # again, not great-great

## evaluate models including additional covariates in S-R model ---------------------------
cor(dat$ReproductiveFemales, dat$Ovig_female_CO, use="p") # 0.31
mod3 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat)

summary(mod3)
MuMIn::AICc(mod1, mod3)

plot(mod3, resid=T, pch=1, se=F, rug=F, pages=1)

# this model is an improvement, but I am suspicious of the dome-shaped CO effect - 
# seems likely to be spurious

## cod--RA3
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.06
mod4 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat)

summary(mod4)
plot(mod4, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod1, mod4) # almost the same

## cod--RA2

cor(dat$ReproductiveFemales, dat$PCod_RA2) # -0.06
mod4a <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA2, k=4),
            data = dat)

summary(mod4a)
plot(mod4a, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod4, mod4a) # almost the same

## cod--RA2
names(dat)
cor(dat$ReproductiveFemales, dat$Pcod_lag1) # -0.06
mod4b <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$Pcod_lag1, k=4),
             data = dat)

summary(mod4b)
plot(mod4b, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod4, mod4a, mod4b) # almost the same


## FHS
cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.001
mod5a <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
            data = dat[dat$releaseyear >= 1983,])

summary(mod5a)
plot(mod5a, resid=T, pch=1, se=F, rug=F, pages=1)

# substitute FHS rolling average
mod5b <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_RA2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod5b)
plot(mod5b, resid=T, pch=1, se=F, rug=F, pages=1)


MuMIn::AICc(mod1b, mod5a, mod5b) # better w/ FHS - by a lot!



## retain FHS, and see if it's possible to improve with AO/wind/other physical covariates?? ----------

# refit 5c with all available data (i.e., not limited to years available for FHS_RA2)
mod5c <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat)

summary(mod5c)

plot(mod5c, resid=T, pch=1, se=F, rug=F, pages=1) # FHS effect is negative/linear = plausible!

## NE wind

mod6 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat)

summary(mod6)
plot(mod6, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod5a, mod6) # adding wind not helpful!
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated


## AO
cor(dat$FHS_lag2, dat$AO_RA3, use="p")
mod7 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
            data = dat)

summary(mod7)

plot(mod7, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod5a, mod7) # adding AO not helpful!


## and bottom temp
cor(dat$FHS_lag2, dat$NBT_3RA, use="p")
mod8 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
            data = dat)

summary(mod8)

plot(mod8, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod5a, mod8) # adding bottom temp only marginally better


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")
mod9 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
            data = dat)

summary(mod9)

plot(mod9, resid=T, pch=1, se=F, rug=F, pages=1) 

## shifts FHS into improbable shape!

## refit with linear FHS effect
mod9b <- gam(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
            data = dat)

summary(mod9b)

plot(mod9b, resid=T, pch=1, se=F, rug=F, pages=1) 
MuMIn::AICc(mod5a, mod9b) # not surprisingly - not any better

## drop FHS and consider only S-R and SST
mod9c <- gam(logRS ~ s(ReproductiveFemales, k=4)  + s(SST_May_July, k=4),
            data = dat)

summary(mod9c)

plot(mod9c, resid=T, pch=1, se=F, rug=F, pages=1) 


## PDO
mod10 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PDO_RA3, k=4),
data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod10)

plot(mod10, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod5a, mod10) # quite an improvement to include PDO

## Drop FHS 
mod10a <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(PDO_RA3, k=4),
             data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod10a)

plot(mod10a, resid=T, pch=1, se=F, rug=F, pages=1)


MuMIn::AICc(mod10, mod10a) # quite an improvement to include PDO

## examine the residuals for mod1
resid10 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod10, type="response"))

ggplot(resid10, aes(year, resid)) +
  geom_line() +
  geom_point() +
geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod10)~1)
acf(resid(mod10)) #pretty good!

MuMIn::AICc(mod1,mod1b,mod2,mod3,mod4,mod5a,mod5b,mod5c,mod6,mod7,mod8,mod9, mod9b, mod9c,mod10,mod10a)    #model 10 improves on next best (model 9) by ~5 AIC points.

## Summaries for all models for quick comparison
summary(mod1)
summary(mod1b)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5a)
summary(mod5b)
summary(mod5c)
summary(mod6)
summary(mod7)
summary(mod8)
summary(mod9)
summary(mod9b)
summary(mod9c)
summary(mod10)
summary(mod10a)
