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

mod1 <- gamm(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat, correlation=corAR1())

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


## evaluate models including additional covariates in S-R model ---------------------------
cor(dat$ReproductiveFemales, dat$Ovig_female_CO, use="p") # 0.31

mod2 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4),
            data = dat, correlation=corAR1())

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

MuMIn::AICc(mod2$gam)
## Re-run ovigerous female CO model but with FHS added 
cor(dat$Ovig_female_CO,dat$FHS_lag2, use="p")
cor(dat$ReproductiveFemales,dat$FHS_lag2, use="p")

mod3 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4)+ s(Ovig_female_CO, k = 4),
             data = dat, correlation=corAR1())
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


## Re-run ovigerous female CO model but with FHS and Pcod added 
cor(dat$Ovig_female_CO,dat$FHS_lag2, use="p")
cor(dat$PCod_RA3,dat$FHS_lag2, use="p")
cor(dat$Ovig_female_CO,PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$FHS_lag2, use="p")

mod4 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4)+ s(Ovig_female_CO, k = 4)
             + s(dat$PCod_RA3, k=4), data = dat, correlation=corAR1())
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


## Re-run ovigerous female CO model but withPcod added 
cor(dat$PCod_RA3,dat$FHS_lag2, use="p")
cor(dat$Ovig_female_CO,PCod_RA3, use="p")
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")

mod5 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(Ovig_female_CO, k = 4)
             + s(dat$PCod_RA3, k=4), data = dat, correlation=corAR1())
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


## Compare AICc values for models
MuMIn::AICc(mod1, mod2,mod3,mod4,mod5,mod3c)


## cod
cor(dat$ReproductiveFemales,dat$PCod_RA3, use="p")
cor(dat$ReproductiveFemales, dat$PCod_RA3) # -0.06

mod6 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat, correlation=corAR1())

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


MuMIn::AICc(mod1, mod6) # almost the same

## FHS
cor(dat$ReproductiveFemales, dat$FHS_lag2, use="p") # 0.001

mod7 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4),
             data = dat, correlation=corAR1())

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


MuMIn::AICc(mod1, mod7) # almost the same

## Now use FHS RA
cor(dat$ReproductiveFemales, dat$FHS_RA2, use="p") # -0.12

mod8 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_RA2, k=4),
             data = dat, correlation=corAR1())

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

## NE wind 


cor(dat$ReproductiveFemales, dat$NE.wind, use="p") # wind and female bairdi not correlated
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated

mod9 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat, correlation=corAR1())

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


AICc(mod1, mod9) # adding wind not helpful!


## SE wind
cor(dat$FHS_lag2, dat$SE.wind, use="p") # SE wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # SE wind and correlated

mod10 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SE.wind, k=4),
             data = dat, correlation=corAR1())

#Model summary
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

AICc(mod1, mod10) # adding wind not helpful!

## AO # This is currently best performing model
cor(dat$FHS_lag2, dat$AO_RA3, use="p") # Not correlated
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # Not correlated

mod11 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
            data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1, mod11) # adding AO changes nothing


## and bottom temp
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p") # Not correlated
cor(dat$FHS_lag2, dat$NBT_3RA, use="p")            # Not correlated

mod12 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
            data = dat, correlation=corAR1())

#Model summaries
summary(mod12)
summary(mod12$gam)
summary(mod12$lme)

#Inspect model object
mod12

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod12$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod12$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod12) # Bah


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p")

mod13 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
            data = dat, correlation=corAR1())

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


## refit with linear FHS effect
## No need to run correlations here as same as above except FHS now linear

mod14 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
             data = dat, correlation=corAR1())

#Model summaries
summary(mod14)
summary(mod14$gam)
summary(mod14$lme)

#Inspect model object
mod14

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod14$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod14$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod14, mod15) # not surprisingly - not any better


## PDO
cor(dat$FHS_lag2, dat$PDO_RA3, use="p")   #somewhat correlated
cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # not correlated

mod15 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(PDO_RA3, k=4),
             data = dat, correlation=corAR1())

#Model summaries
summary(mod15)
summary(mod15$gam)
summary(mod15$lme)

#Inspect model object
mod15

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod15$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod15$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod15) # No improvement


##########################################################################################
################################# try linear terms #######################################
## refit with linear FHS effect ONLY
## No need for correlation calculations as repeat of prior model with modded term

mod16 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod16)
summary(mod16$gam)
summary(mod16$lme)

#Inspect model object
mod16

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod16$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod16$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod16) # not surprisingly - not any better

## evaluate models including additional covariates in S-R model ---------------------------


mod17 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + Ovig_female_CO,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod17)
summary(mod17$gam)
summary(mod17$lme)

#Inspect model object
mod17

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod17$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod17$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod17)

## cod

mod18 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + PCod_RA3,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod18)
summary(mod18$gam)
summary(mod18$lme)

#Inspect model object
mod18

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod18$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod18$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod18) # almost the same


## NE wind


mod19 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + NE.wind,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod19)
summary(mod19$gam)
summary(mod19$lme)

#Inspect model object
mod19

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod19$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod19$lme, resid=T, pch=19, rug=F, se=F, pages=1)


AICc(mod1, mod19) # adding wind not helpful, even makes model worse!



## SE wind
mod20 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SE.wind, k=4),
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod20)
summary(mod20$gam)
summary(mod20$lme)

#Inspect model object
mod20

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod20$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod20$lme, resid=T, pch=19, rug=F, se=F, pages=1)

AICc(mod1, mod20) # adding wind not helpful!


## AO

mod21 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + AO_RA3,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod21)

summary(mod21$gam)
summary(mod21$lme)

#Inspect model object
mod21

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod21$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod21$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod21) # adding AO not helpful!


## and bottom temp rolling average

mod22 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + NBT_3RA,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod22)
summary(mod22$gam)
summary(mod22$lme)

#Inspect model object
mod22

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod22$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod22$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod22) 


## sst

mod23 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + SST_May_July,
              data = dat[dat$releaseyear >= 1981,], correlation=corAR1(), na.action=na.omit)

#Model summaries
summary(mod23)
summary(mod23$gam)
summary(mod23$lme)

#Inspect model object
mod23

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod23$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod23$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod23)


## PDO
mod24 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3,
              data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod24)
summary(mod24$gam)
summary(mod24$lme)

#Inspect model object
mod24

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod24$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod24$lme, resid=T, pch=19, rug=F, se=F, pages=1)


MuMIn::AICc(mod1, mod24) # 

## Add Pacific cod combined with FHS

mod25 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PCod_RA3,
              data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod25)
summary(mod25$gam)
summary(mod25$lme)

#Inspect model object
mod25

#plot
dev.new()
par(mfrow=c(3,1))

plot(mod25$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod25$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod25)

## PDO
mod26 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3 + AO_RA3,
              data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod26)
summary(mod26$gam)
summary(mod26$lme)

#Inspect model object
mod26

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod26$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod26$lme, resid=T, pch=19, rug=F, se=F, pages=1)


MuMIn::AICc(mod1, mod26) # 

####### TRY NEW 3 YEAR MIN NBT VARIABLE ##########################
## Add Pacific cod combined with FHS

mod27 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + NBT_3yr_minTemp,
              data = dat[dat$releaseyear >= 1981,], na.action = na.omit, correlation=corAR1())

#Model summaries
summary(mod27)
summary(mod27$gam)
summary(mod27$lme)

#Inspect model object
mod27

#plot
dev.new()
par(mfrow=c(3,1))

plot(mod27$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod27$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod27)
##### Wrap up 
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6, mod7,mod8,mod10,mod11,mod12,mod13,mod14,mod15,
            mod16,mod17,mod18,mod19,mod20,mod21, mod22,mod23,mod24,mod25,mod26) #  

names(dat)
cor.dat <- dat %>%
  select(-era,-releaseyear)

c<-cor(cor.dat, use="p")
c
write.csv(c,"output/Variable_correlations_AICc_Era2.csv")
