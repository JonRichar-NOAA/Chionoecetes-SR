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
## NOTE: Previously, model 7 was repeated as model 9 to start this sequence...this has been removed 
## and in Excel sheet model 10 is reported as model 9, model 11 as 10, etc.

cor(dat$ReproductiveFemales, dat$NE.wind, use="p") # wind and female bairdi not correlated
cor(dat$FHS_lag2, dat$NE.wind, use="p") # wind and FHS not correlated

mod10 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NE.wind, k=4),
            data = dat, correlation=corAR1())

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


AICc(mod1, mod10) # adding wind not helpful!


## SE wind
cor(dat$FHS_lag2, dat$SE.wind, use="p") # SE wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # SE wind and correlated

mod11 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SE.wind, k=4),
             data = dat, correlation=corAR1())

#Model summarie
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

AICc(mod1, mod11) # adding wind not helpful!

## AO # This is currently best performing model
cor(dat$FHS_lag2, dat$AO_RA3, use="p") # Not correlated
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # Not correlated

mod12 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(AO_RA3, k=4),
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

MuMIn::AICc(mod1, mod12) # adding AO changes nothing


## and bottom temp
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p") # Not correlated
cor(dat$FHS_lag2, dat$NBT_3RA, use="p")            # Not correlated

mod13 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(NBT_3RA, k=4),
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

MuMIn::AICc(mod1, mod13) # Bah


## sst
cor(dat$FHS_lag2, dat$SST_May_July, use="p")
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p")

mod14 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(SST_May_July, k=4),
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



## shifts FHS into improbable shape!

## refit with linear FHS effect
## No need to run correlations here as same as above except FHS now linear

mod15 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(SST_May_July, k=4),
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

gam.check(mod15) 

MuMIn::AICc(mod1, mod14, mod15) # not surprisingly - not any better


## PDO
cor(dat$FHS_lag2, dat$PDO_RA3, use="p")   #somewhat correlated
cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # not correlated

mod16 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=4) + s(PDO_RA3, k=4),
             data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1, mod16) # No improvement


####################### Re-run FHS inclusive models substituting Pcod ############
## cod

cor(dat$ReproductiveFemales, dat$PCod_RA3) # not correlated

mod17 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4),
            data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1, mod17) # almost the same


## retain Pcod, and see if it's possible to improve with AO/wind/other physical covariates?? ----------

## NE wind
cor(dat$ReproductiveFemales, dat$NE.wind) # not correlated
cor(dat$PCod_RA3, dat$NE.wind) # not correlated

mod18 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(NE.wind, k=4),
            data = dat, correlation=corAR1())

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

AICc(mod1, mod18) # adding wind not helpful!

cor(dat$PCod_RA3, dat$NE.wind, use="p") # wind and FHS not correlated

## SE wind
cor(dat$PCod_RA3, dat$SE.wind, use="p") # wind and FHS not correlated
cor(dat$ReproductiveFemales, dat$SE.wind, use="p") # correlated

mod19 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(SE.wind, k=4),
             data = dat, correlation=corAR1())

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

AICc(mod1, mod19) # adding wind not helpful!



## AO
cor(dat$PCod_RA3, dat$AO_RA3, use="p")               # Not correlated
cor(dat$ReproductiveFemales, dat$AO_RA3, use="p")    # Not correlated

mod20 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(AO_RA3, k=4),
            data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1, mod20) # adding AO not helpful!


## and bottom temp
cor(dat$PCod_RA3, dat$NBT_3RA, use="p")                 # Not correlated
cor(dat$ReproductiveFemales, dat$NBT_3RA, use="p")      # Not correlated

mod21 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(NBT_3RA, k=4),
            data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1, mod21) # adding bottom temp only marginally better


## sst
cor(dat$ReproductiveFemales, dat$SST_May_July, use="p") # Not correlated
cor(dat$PCod_RA3, dat$SST_May_July, use="p")            # STRONGLY negatively correlated

mod22 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(SST_May_July, k=4),
            data = dat, correlation=corAR1())

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


## refit with linear Pcod effect
mod23 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + PCod_RA3 + s(SST_May_July, k=4),
             data = dat, correlation=corAR1())

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


MuMIn::AICc(mod1, mod16, mod23) # not surprisingly - not any better


## PDO

cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # Not correlated
cor(dat$PCod_RA3, dat$PDO_RA3, use="p")            # STRONGLY negatively correlated

mod24 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(dat$PCod_RA3, k=4) + s(PDO_RA3, k=4),
             data = dat, correlation=corAR1())

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


MuMIn::AICc(mod1,mod16, mod24) # quite an improvement to include PDO



##### AO only #########################################################################################
## AO # This is currently best performing model

cor(dat$ReproductiveFemales, dat$AO_RA3, use="p") # Not correlated

mod25 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(AO_RA3, k=4),
             data = dat, correlation=corAR1())

#Model summaries
summary(mod25)
summary(mod25$gam)
summary(mod25$lme)

#Inspect model object
mod25

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod25$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod25$lme, resid=T, pch=19, rug=F, se=F, pages=1)


MuMIn::AICc(mod1, mod12, mod25) # 

## PDO only

cor(dat$ReproductiveFemales, dat$PDO_RA3, use="p") # Not correlated

mod26 <- gamm(logRS ~ s(ReproductiveFemales, k=4) +  s(PDO_RA3, k=4),
              data = dat, correlation=corAR1())

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

MuMIn::AICc(mod1,mod16, mod26) # quite an improvement to include PDO

## refit with linear FHS effect and AO
## No need to run correlations here as same as above except FHS now linear

mod27 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(AO_RA3, k=4),
              data = dat, correlation=corAR1())

#Model summaries
summary(mod27)
summary(mod27$gam)
summary(mod27$lme)

#Inspect model object
mod27

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod27$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod27$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod27) # Good!

## refit with linear FHS effect and AO
## No need to run correlations here as same as above except FHS now linear

mod28 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + s(PDO_RA3, k=4),
              data = dat, correlation=corAR1())

#Model summaries
summary(mod28)
summary(mod28$gam)
summary(mod28$lme)

#Inspect model object
mod28

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod28$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod28$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod28)


## refit with linear FHS effect ONLY

mod29 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2,
              data = dat, correlation=corAR1())

#Model summaries
summary(mod29)
summary(mod29$gam)
summary(mod29$lme)

#Inspect model object
mod29

#plot
dev.new()
par(mfrow=c(2,1))

plot(mod29$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod29$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod1, mod29)

##### Wrap up 
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6, mod7,mod8,mod10,mod11,mod12,mod13,mod14,mod15,
            mod16,mod17,mod18,mod19,mod20,mod21, mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29) #mods 7 and 9 are repeats, Model 12 is best , model 29 is worse b 0.3 AICc units  

names(dat)
cor.dat <- dat %>%
  select(-era,-releaseyear)

c<-cor(cor.dat, use="p")
c
write.csv(c,"output/Variable_correlations_AICc_Era2.csv")
