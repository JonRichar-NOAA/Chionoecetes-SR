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
  subset(Era_AICc==2)->dat1

dat1
dat1 %>%
  select(releaseyear:NE.wind) ->dat

head(dat)

dat
cor.dat <- dat %>%
  select(-Era_AICc)

cor(cor.dat, use="p")


## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier

mod1 <- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat)

summary(mod1)
plot(mod1, resid=T, pch=19, rug=F, se=F)

dat$era <- as.factor(if_else(dat$releaseyear <= 1995, 1, 2))


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

mod2 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat)

summary(mod2)
MuMIn::AICc(mod1, mod2)

plot(mod2, resid=T, pch=1, se=F, rug=F, pages=1)

## Re-run ovigerous female CO model but with FHS added 
cor(dat$Ovig_female_CO,dat$FHS_lag2, use="p")
cor(dat$ReproductiveFemales,dat$FHS_lag2, use="p")

mod3 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4)+ s(Ovig_female_CO, k = 4),
             data = dat[dat$releaseyear >= 1983,])
summary(mod3)


plot(mod3, resid=T, pch=1, se=F, rug=F, pages=1)

## Re-run ovigerous female CO model but with FHS and Pcod added 
cor(dat$Ovig_female_CO,dat$FHS_lag2, use="p")
cor(dat$PCod_RA3,dat$FHS_lag2, use="p")
cor(dat$Ovig_female_CO,PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$FHS_lag2, use="p")

mod4 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4)+ s(Ovig_female_CO, k = 4)
             + s(dat$PCod_RA3, k=4), data = dat[dat$releaseyear >= 1983,])
summary(mod4)
plot(mod4, resid=T, pch=1, se=F, rug=F, pages=1)

## Re-run ovigerous female CO model but withPcod added 
cor(dat$PCod_RA3,dat$FHS_lag2, use="p")
cor(dat$Ovig_female_CO,PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")

mod5 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4)
             + s(dat$PCod_RA3, k=4), data = dat[dat$releaseyear >= 1983,])
summary(mod5)
plot(mod5, resid=T, pch=1, se=F, rug=F, pages=1)


## Compare AICc values for models
MuMIn::AICc(mod1, mod2,mod3,mod4,mod5,mod3c)


# this model is an improvement, but I am suspicious of the dome-shaped CO effect - 
# seems likely to be spurious

## cod
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.06

mod6 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat)

summary(mod6)

plot(mod6, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod1, mod4) # almost the same

## FHS
cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.001

mod7 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod7)

## Now use FHS RA
cor(dat$ReproductiveFemales, dat$FHS_RA2, use="p") # -0.12

mod8 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_RA2, k=4),
             data = dat[dat$releaseyear >= 1983,])

summary(mod8)

plot(mod8, resid=T, pch=1, se=F, rug=F, pages=1)

## examine the residuals for mod7
resid8 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod8))

ggplot(resid7, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod8)~1)
acf(resid(mod8)) # again, not great-great

## retain FHS, and see if it's possible to improve with AO/wind/other physical covariates?? ----------

# refit 5a with all available data (i.e., not limited to years available for FHS_RA2)
cor(dat$ReproductiveFemales,dat$FHS_lag2, use="p")

mod9 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat)

summary(mod9)

plot(mod9, resid=T, pch=1, se=F, rug=F, pages=1) # FHS effect is negative/linear = plausible!

## NE wind
cor(dat$ReproductiveFemales, dat$NE.wind, use="p") # wind and female bairdi not correlated
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated

mod10 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat)

summary(mod10)
plot(mod10, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod1, mod10) # adding wind not helpful!


## SE wind
cor(dat$FHS_lag2, dat$SE.wind, use="p") # SE wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # SE wind and correlated

mod11 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SE.wind, k=4),
             data = dat)

summary(mod11)
plot(mod11, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod1, mod11) # adding wind not helpful!



## AO
cor(dat$FHS_lag2, dat$AO_RA3, use="p") # Not correlated
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # Not correlated

mod12 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
            data = dat)

summary(mod12)

plot(mod12, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod12) # adding AO very helpful!


## and bottom temp
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p") # Not correlated
cor(dat$FHS_lag2, dat$NBT_3RA, use="p")            # Not correlated

mod13 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
            data = dat)

summary(mod13)

plot(mod13, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod13) # adding bottom temp only marginally better


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p")

mod14 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
            data = dat)

summary(mod14)

plot(mod14, resid=T, pch=1, se=F, rug=F, pages=1) 

## shifts FHS into improbable shape!

## refit with linear FHS effect
## No need to run correlations here as same as above except FHS now linear

mod15 <- gam(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
             data = dat)

summary(mod15)

plot(mod15, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod14, mod15) # not surprisingly - not any better


## PDO
cor(dat$FHS_lag2, dat$PDO_RA3, use="p")   #somewhat correlated
cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # not correlated

mod16 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3) + s(PDO_RA3, k=4),
             data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod16)

plot(mod16, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod16) # quite an improvement to include PDO

## examine the residuals for mod1
resid16 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod16, type="response"))

ggplot(resid16, aes(year, resid)) +
  geom_line() +
  geom_point() +
  geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod16)~1)
acf(resid(mod16)) #pretty good!

####################### Re-run FHS inclusive models substituting Pcod ############
## cod

cor(dat$ReproductiveFemales, dat$PCod_RA3) # not correlated

mod17 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat)

summary(mod17)

plot(mod17, resid=T, pch=1, se=F, rug=F, pages=1)

MuMIn::AICc(mod1, mod17) # almost the same


## retain Pcod, and see if it's possible to improve with AO/wind/other physical covariates?? ----------

## NE wind
cor(dat$ReproductiveFemales, dat$NE.wind) # not correlated
cor(dat$PCod_RA3, dat$NE.wind) # not correlated

mod18 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(NE.wind, k=4),
            data = dat)

summary(mod18)
plot(mod18, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod1, mod18) # adding wind not helpful!

cor(dat$PCod_RA3, dat$NE.wind, use="p") # wind and FHS not correlated

## SE wind
cor(dat$PCod_RA3, dat$SE.wind, use="p") # wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # correlated

mod19 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(SE.wind, k=4),
             data = dat)

summary(mod19)
plot(mod19, resid=T, pch=1, se=F, rug=F, pages=1) 

AICc(mod1, mod19) # adding wind not helpful!



## AO
cor(dat$PCod_RA3, dat$AO_RA3, use="p")               # Not correlated
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p")    # Not correlated

mod20 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(AO_RA3, k=4),
            data = dat)

summary(mod20)

plot(mod20, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod20) # adding AO not helpful!


## and bottom temp
cor(dat$PCod_RA3, dat$NBT_3RA, use="p")                 # Not correlated
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p")      # Not correlated

mod21 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(NBT_3RA, k=4),
            data = dat)

summary(mod21)

plot(mod21, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1, mod21) # adding bottom temp only marginally better


## sst
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p") # Not correlated
cor(dat$PCod_RA3, dat$SST_May_July, use="p")            # STRONGLY negatively correlated

mod22 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(SST_May_July, k=4),
            data = dat)

summary(mod22)

plot(mod22, resid=T, pch=1, se=F, rug=F, pages=1) 

## shifts FHS into improbable shape!

## refit with linear Pcod effect
mod23 <- gam(logRS ~ s(ReproductiveFemales, k=4) + PCod_RA3 + s(SST_May_July, k=4),
             data = dat)

summary(mod23)

plot(mod23, resid=T, pch=1, se=F, rug=F, pages=1) 
MuMIn::AICc(mod1, mod16, mod23) # not surprisingly - not any better


## PDO

cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # Not correlated
cor(dat$PCod_RA3, dat$PDO_RA3, use="p")            # STRONGLY negatively correlated

mod24 <- gam(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(PDO_RA3, k=4),
             data = dat[dat$releaseyear >= 1983,], na.action = "na.fail")

summary(mod24)

plot(mod24, resid=T, pch=1, se=F, rug=F, pages=1) 

MuMIn::AICc(mod1,mod16, mod24) # quite an improvement to include PDO

## examine the residuals for mod24
resid24 <- data.frame(year=dat$releaseyear[dat$releaseyear >= 1983],
                      resid=residuals(mod24, type="response"))

ggplot(resid23, aes(year, resid)) +
  geom_line() +
  geom_point() +
  geom_smooth(method="gam", formula = y ~ s(x, k=4))

dwtest(resid(mod24)~1)
acf(resid(mod24)) #


MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6, mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,
            mod16,mod17,mod18,mod19,mod20,mod21, mod22,mod23,mod24) #mods 7 and 9 are repeats, Model 12 is best   

names(dat)
cor.dat <- dat %>%
  select(-era,-releaseyear)

cor(cor.dat, use="p")

write.csv(cor.dat,"output/Variable_correlations_AICc_Era2.csv")
