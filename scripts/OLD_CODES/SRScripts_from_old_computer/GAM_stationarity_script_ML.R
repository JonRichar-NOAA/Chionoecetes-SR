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

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

# load NE wind
wind <- read.csv("./data/mean.tanner.wind.csv", row.names = 1)
names(wind) <- c("releaseyear", "NE.wind")

dat <- left_join(dat, wind)

hist(dat$logRS)

cor.dat <- dat %>%
  select(-era)

cor(cor.dat, use="p")


## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier
mod1 <- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat)

summary(mod1)
plot(mod1, resid=T, pch=19, rug=F, se=F)

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

## evaluate models including additional covariates in S-R model ---------------------------
cor(dat$ReproductiveFemales, dat$Ovig_female_CO, use="p") # 0.31
mod3 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat)

summary(mod3)
MuMIn::AICc(mod1, mod3)

plot(mod3, resid=T, pch=1, se=F, rug=F, pages=1)
# this model is an improvement, but I am suspicious of the dome-shaped CO effect - 
# seems likely to be spurious

## cod
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.06
mod4 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat)

summary(mod4)
MuMIn::AICc(mod1, mod4) # almost the same

## FHS
cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.001
mod5a <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
            data = dat[dat$releaseyear >= 1983,])

summary(mod5a)

mod5b <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_RA2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod5b)

# fit mod1 without years for which FHS_RA2 is unavailable
mod1restricted <- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat[dat$releaseyear >= 1983,])

summary(mod1restricted)

MuMIn::AICc(mod1restricted, mod5a, mod5b) # better w/ FHS - by a lot!

## examine the residuals for mod1
resid5a <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                     resid=residuals(mod5))

ggplot(resid5a, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod5a)~1)
acf(resid(mod5a)) # again, not great-great

## retain FHS, and see if it's possible to improve with AO/wind/other physical covariates?? ----------

# refit 5a with all available data (i.e., not limited to years available for FHS_RA2)
mod5a <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat)

summary(mod5a)

plot(mod5a, resid=T, pch=1, se=F, rug=F, pages=1) # FHS effect is negative/linear = plausible!

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


## PDO
mod10 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PDO_RA3, k=4),
data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod10)

plot(mod10, resid=T, pch=1, se=F, rug=F, pages=1) 
MuMIn::AICc(mod5a, mod10) # quite an improvement to include PDO

## examine the residuals for mod1
resid10 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod10, type="response"))

ggplot(resid10, aes(year, resid)) +
  geom_line() +
  geom_point() +
geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod10)~1)
acf(resid(mod10)) #pretty good!