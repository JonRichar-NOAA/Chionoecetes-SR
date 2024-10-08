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

dat <- read.csv("./data/EBS_Crab_and_envar_data_full_extent_for_analysis_reducedTime_matchLAG4.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

####################################################### Run models#########################################################

###################################################### S-R effect only ##################################################
############################################ GAMM ##################################################################
mod1 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3),
            data = dat, correlation=corAR1())
#Model summaries
summary(mod1)
summary(mod1$gam)
summary(mod1$lme)

#Inspect model object
mod1

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod1$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod1$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### S-R effect only ##################################################
############################################ GAM ##################################################################
mod1a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3),
             data = dat)
#Model summaries
summary(mod1a)
summary(mod1a$gam)
summary(mod1a$lme)

#Inspect model object
mod1a

#plot,
dev.new()
par(mfrow=c(1,1))

plot(mod1a, resid=T, pch=19, rug=F, se=F)
#plot(mod1$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Add ovigerous female opilio for competitive effect amongst juveniles ###

################################################### Linear ######################################
mod2 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + Ovig_female_CO,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod2)
summary(mod2$gam)
summary(mod2$lme)

#Inspect model object
mod2

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod2$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod2)
################################################### Non-linear #########################################
mod2a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(Ovig_female_CO,k=3),
             data = dat)
#Model summaries
summary(mod2a)
summary(mod2a$gam)
summary(mod2a$lme)

#Inspect model object
mod2a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod2a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod2a)
############################################################################################################################
################################################## Predation models ########################################################

###################################################### Add Pcod lag1 only ##################################################
# NOTE: All cod models have dome shape in GAM models

################################################### linear ############################################
mod3 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + Pcod_lag1,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod3)
summary(mod3$gam)
summary(mod3$lme)

#Inspect model object
mod3

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod3$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod3$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod3)
################################################### GAM Non-linear ############################################
mod3a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(Pcod_lag1,k=3),
             data = dat)
#Model summaries
summary(mod3a)

#Inspect model object
mod3a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod3a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod3$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod3a)

################################################### linear ############################################
mod3b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(Pcod_lag1,k=3),
             data = dat, correlation=corAR1())
#Model summaries
summary(mod3b)
summary(mod3b$gam)
summary(mod3b$lme)

#Inspect model object
mod3b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod3b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod3$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod3b)
###################################################### Add PCod_RA2 only ##################################################


############################## linear #############################################
mod4 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PCod_RA2,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod4)
summary(mod4$gam)
summary(mod4$lme)

#Inspect model object
mod4

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod4$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod4$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod4)
################################### GAM non-linear ##############################################
mod4a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PCod_RA2,k=3),
              data = dat)
#Model summaries
summary(mod4a)

#Inspect model object
mod4a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod4a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod4$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod4a)
############################## GAMM-nonlinear #############################################
mod4b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PCod_RA2,k=3),
             data = dat, correlation=corAR1())
#Model summaries
summary(mod4b)
summary(mod4b$gam)
summary(mod4b$lme)

#Inspect model object
mod4b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod4b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod4$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod4b)
###################################################### Add PCod_RA3 only ##################################################
###################################################### linear ###############################################
mod5 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PCod_RA3,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod5)
summary(mod5$gam)
summary(mod5$lme)

#Inspect model object
mod5

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod5$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod5$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5)
###################################################### GAM-nonlinear ###############################################
mod5a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PCod_RA3,k=3),
             data = dat)
#Model summaries
summary(mod5a)

#Inspect model object
mod5a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod5a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(
MuMIn::AICc(mod5a)
###################################################### GAMM nonlinear ###############################################
mod5b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PCod_RA3,k=3),
             data = dat, correlation=corAR1())
#Model summaries
summary(mod5b)
summary(mod5b$gam)
summary(mod5b$lme)

#Inspect model object
mod5b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod5b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod5$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod5b)

###########################################################################################################################
AICc(mod3,mod3a,mod3b,mod4,mod4a,mod4b,mod5,mod5a,mod5b)
###################################################### Add FHS lag 2 only ##################################################
###################################################### linear #######################################################
mod6 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod6)
summary(mod6$gam)
summary(mod6$lme)

#Inspect model object
mod6

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod6$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod6$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod6)
###################################### non-linear ###################################################################################
mod6a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3),
             data = dat)
#Model summaries
summary(mod6a)


#Inspect model object
mod6a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod6a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod6$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Add FHS RA2 only ###################################################
##################################################### linear ###########################################################
mod7 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_RA2,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod7)
summary(mod7$gam)
summary(mod7$lme)

#Inspect model object
mod7

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod7$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod7$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### non-linear ###########################################################
mod7a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_RA2,k=3),
             data = dat)
#Model summaries
summary(mod7a)

#Inspect model object
mod7a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod7a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod7$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### Combine Pcod_lag1 and FHS lag 2 #####################################
##################################################### linear ##########################################################
mod8 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + Pcod_lag1 + FHS_lag2,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod8)
summary(mod8$gam)
summary(mod8$lme)

#Inspect model object
mod8

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod8$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod8$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### non-linear ##########################################################
mod8a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(Pcod_lag1, k =3) + s(FHS_lag2,k=3),
             data = dat)
#Model summaries
summary(mod8a)

#Inspect model object
mod8a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod8a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod8$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod8a)


##################################################### GAMM nonlinear ##########################################################
mod8b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(Pcod_lag1,k=3) + s(FHS_lag2,k=3),
             data = dat, correlation=corAR1())
#Model summaries
summary(mod8b)
summary(mod8b$gam)
summary(mod8b$lme)

#Inspect model object
mod8b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod8b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#

#######################################################################################################################
MuMIn::AICc(mod8,mod8a,mod8b)

##################################################### Combine PCod_RA2 and FHS_RA2 #####################################
##################################################### linear #######################################################
mod9 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PCod_RA2 + FHS_RA2,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod9)
summary(mod9$gam)
summary(mod9$lme)

#Inspect model object
mod9

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod9$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod9$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### non-linear #######################################################
mod9a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PCod_RA2, k=3) + s(FHS_RA2,k=3),
             data = dat)
#Model summaries
summary(mod9a)

#Inspect model object
mod9a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod9a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod9$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod9a)
########################################################################################################################
##################################################### Environmental factors ############################################

##################################################### NBT 3 yr rolling average #########################################
##################################################### linear ###########################################################
mod10 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + NBT_3RA,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod10)
summary(mod10$gam)
summary(mod10$lme)

#Inspect model object
mod10

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod10$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod10$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### noninear ###########################################################
mod10a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(NBT_3RA,k=3),
              data = dat)
#Model summaries
summary(mod10a)

#Inspect model object
mod10a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod10a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod10$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod10a)
##################################################### NBT 3 yr min temperature #########################################
##################################################### linear ###########################################################
mod11 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + NBT_3yr_minTemp,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod11)
summary(mod11$gam)
summary(mod11$lme)

#Inspect model object
mod11

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod11$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod11$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear ###########################################################
mod11a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(NBT_3yr_minTemp,k=3),
              data = dat)
#Model summaries
summary(mod11a)

#Inspect model object
mod11a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod11a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod11$lme, resid=T, pch=19, rug=F, se=F, pages=1)


##################################################### AO_RA2 #########################################
##################################################### linear ##########################################
mod12 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + AO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod12)
summary(mod12$gam)
summary(mod12$lme)

#Inspect model object
mod12

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod12$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod12$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### non-linear ##########################################
mod12a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(AO_RA2,k=3),
              data = dat)
#Model summaries
summary(mod12a)

#Inspect model object
mod12a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod12a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod12$lme, resid=T, pch=19, rug=F, se=F, pages=1)


##################################################### AO_RA3 #########################################
##################################################### linear #########################################
mod13 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + AO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod13)
summary(mod13$gam)
summary(mod13$lme)

#Inspect model object
mod13

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod13$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod13$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear #########################################
mod13a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(AO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod13a)

#Inspect model object
mod13a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod13a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod13$lme, resid=T, pch=19, rug=F, se=F, pages=1)


##################################################### PDO_RA2 #########################################
##################################################### linear ##########################################
mod14 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PDO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod14)
summary(mod14$gam)
summary(mod14$lme)

#Inspect model object
mod14

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod14$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod14$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### non-linear ##########################################
mod14a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PDO_RA2,k=3),
              data = dat)
#Model summaries
summary(mod14a)

#Inspect model object
mod14a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod14a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod14$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod14a)
##################################################### PDO_RA3 #########################################
##################################################### linear ##########################################
mod15 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PDO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod15)
summary(mod15$gam)
summary(mod15$lme)

#Inspect model object
mod15

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod15$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod15$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod15)
##################################################### nonlinear ##########################################
mod15a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PDO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod15a)

#Inspect model object
mod15a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod15a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod15$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod15a)

##################################################### SST_May_July #########################################
##################################################### linear ###############################################

mod16 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + SST_May_July,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod16)
summary(mod16$gam)
summary(mod16$lme)

#Inspect model object
mod16

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod16$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod16$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear ###############################################

mod16a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(SST_May_July,k=3),
              data = dat)
#Model summaries
summary(mod16a)

#Inspect model object
mod16a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod16a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod16$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### SE.wind #########################################
##################################################### linear ##########################################

mod17 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + SE.wind,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod17)
summary(mod17$gam)
summary(mod17$lme)

#Inspect model object
mod17

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod17$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod17$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear ##########################################

mod17a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(SE.wind,k=3),
              data = dat)
#Model summaries
summary(mod17a)

#Inspect model object
mod17a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod17a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod17$lme, resid=T, pch=19, rug=F, se=F, pages=1)


##################################################### NW.wind #########################################
##################################################### linear ##########################################
mod18 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + NW.wind,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod18)
summary(mod18$gam)
summary(mod18$lme)

#Inspect model object
mod18

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod18$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod18$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear ##########################################
mod18a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(NW.wind,k=3),
              data = dat)
#Model summaries
summary(mod18a)

#Inspect model object
mod18a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod18a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod18$lme, resid=T, pch=19, rug=F, se=F, pages=1)

############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)

MuMIn::AICc(mod1a,mod2a,mod3a,mod4a, mod5a,mod6a,mod7a,mod8a,mod9a,mod10a,mod11a,mod12a,mod13a,mod14a,mod15a,mod16a,mod17a,mod18a)

###############################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + competition/environmental factor######################


###################################################### Combine FHS lag 2 and ovigerous female opilio  #########################################################
###################################################### linear #################################################################################################
mod19 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO,
             data = dat, correlation=corAR1())
#Model summaries
summary(mod19)
summary(mod19$gam)
summary(mod19$lme)

#Inspect model object
mod19

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod19$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod19$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #################################################################################################
mod19a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) +s(FHS_lag2, k=3) + s(Ovig_female_CO,k=3),
              data = dat)
#Model summaries
summary(mod19a)


#Inspect model object
mod19a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod19a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod19$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod19a)
###################################################### Combine FHS lag 2, ovigerous female opilio + PDO_RA3 #####################################################
###################################################### linear ###################################################################################################
mod20 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + PDO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod20)
summary(mod20$gam)
summary(mod20$lme)

#Inspect model object
mod20

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod20$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod20$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear ###################################################################################################
mod20a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO,k=3) + s(PDO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod20a)

#Inspect model object
mod20a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod20a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod20$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, ovigerous female opilio + AO_RA3 #####################################################
##################################################### linear #################################################################################################

mod21 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + AO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod21)
summary(mod21$gam)
summary(mod21$lme)

#Inspect model object
mod21

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod21$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod21$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear #################################################################################################

mod21a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO, k=3) + s(AO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod21a)


#Inspect model object
mod21a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod21a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod21$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
###################################################### linear #################################################################################
mod22 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + PDO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod22)
summary(mod22$gam)
summary(mod22$lme)

#Inspect model object
mod22

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod22$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod22$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #################################################################################
mod22a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO, k=3) + s(PDO_RA2,k=3),
              data = dat)
#Model summaries
summary(mod22a)


#Inspect model object
mod22a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod22a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod22$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, AO_RA2  #####################################################
###################################################### linear #########################################################################
mod23 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + AO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod23)
summary(mod23$gam)
summary(mod23$lme)

#Inspect model object
mod23

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod23$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod23$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #########################################################################
mod23a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO, k=3) + s(AO_RA2,k=3),
              data = dat)
#Model summaries
summary(mod23a)

#Inspect model object
mod23a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod23a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod23$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, ovigerous female opilio + SST_May_July ################################################
###################################################### linear ###########################################################################################

mod24 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + SST_May_July,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod24)
summary(mod24$gam)
summary(mod24$lme)

#Inspect model object
mod24

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod24$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod24$lme, resid=T, pch=19, rug=F, se=F, pages=1)

#################################################### nonlinear ###########################################################################################

mod24a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO,k=3) + s(SST_May_July,k=3),
              data = dat)
#Model summaries
summary(mod24a)

#Inspect model object
mod24a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod24a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod24$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, ovigerous female opilio + SE wind ################################################
##################################################### linear #############################################################################################
mod25 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + SE.wind,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod25)
summary(mod25$gam)
summary(mod25$lme)

#Inspect model object
mod25

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod25$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod25$lme, resid=T, pch=19, rug=F, se=F, pages=1)

##################################################### nonlinear #############################################################################################
mod25a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO,k=3) + s(SE.wind,k=3),
              data = dat)
#Model summaries
summary(mod25a)

#Inspect model object
mod25a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod25a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod25$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, ovigerous female opilio + NE wind ################################################
###################################################### linear ##############################################################################################
mod26 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + Ovig_female_CO + NW.wind,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod26)
summary(mod26$gam)
summary(mod26$lme)

#Inspect model object
mod26

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod26$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod26$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear ##############################################################################################
mod26a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(Ovig_female_CO,k=3) + s(NW.wind,k=3),
              data = dat)
#Model summaries
summary(mod26a)

#Inspect model object
mod26a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod26a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod26$lme, resid=T, pch=19, rug=F, se=F, pages=1)

############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26)
MuMIn::AICc(mod1a,mod2a,mod3a,mod4a, mod5a,mod6a,mod7a,mod8a,mod9a,mod10a,mod11a,mod12a,mod13a,mod14a,mod15a,mod16a,mod17a,mod18a,mod19a,mod20a,mod21a,mod22a,mod23a,mod24a,mod25a,mod26a)


###############################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + environmental factor##################################

###################################################### Combine FHS lag 2 + PDO_RA3 #####################################################
###################################################### linear ##########################################################################
mod27 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod27)
summary(mod27$gam)
summary(mod27$lme)

#Inspect model object
mod27

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod27$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod27$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear ##########################################################################
mod27a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(PDO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod27a)

#Inspect model object
mod27a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod27a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod27$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### GAMM nonlinear PDO ##########################################################################
mod27b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + s(PDO_RA3,k=3),
              data = dat, correlation=corAR1())
#Model summaries
summary(mod27b)
summary(mod27b$gam)
summary(mod27b$lme)

#Inspect model object
mod27b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod27b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod27$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2 + AO_RA3 #####################################################
###################################################### linear #############################################################################
mod28 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + AO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod28)
summary(mod28$gam)
summary(mod28$lme)

#Inspect model object
mod28

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod28$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod28$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #############################################################################
mod28a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(AO_RA3,k=3),
              data = dat)
#Model summaries
summary(mod28a)

#Inspect model object
mod28a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod28a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod28$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2 + PDO_RA2 #####################################################
###################################################### linear ##########################################################################

mod29 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod29)
summary(mod29$gam)
summary(mod29$lme)

#Inspect model object
mod29

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod29$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod29)
#plot(mod29$lme, resid=T, pch=19, rug=F, se=F, pages=1)
################################ non-linear ######################################################
mod29a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2, k=3) + s(PDO_RA2,k=3),
              data = dat)

summary(mod29a)
#summary(mod29a)$dev.expl

dev.new()
par(mfrow=c(2,2))

plot(mod29a, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod29a)

################################ non-linear ######################################################
mod29b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2, k=3) + s(PDO_RA2,k=3),
               data = dat, correlation=corAR1())
#Model summaries
summary(mod29b)
summary(mod29b$gam)
summary(mod29b$lme)

#Inspect model object
mod29b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod29b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod29b)

###################################################### Combine FHS lag 2, AO_RA2  #####################################################
###################################################### linear #########################################################################
mod30 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + AO_RA2,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod30)
summary(mod30$gam)
summary(mod30$lme)

#Inspect model object
mod30

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod30$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod30$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #########################################################################
mod30a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(AO_RA2,k=3),
              data = dat)
#Model summaries
summary(mod30a)

#Inspect model object
mod30a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod30a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod30$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2 + SST_May_July ################################################
###################################################### linear ##########################################################################
mod31 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + SST_May_July,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod31)
summary(mod31$gam)
summary(mod31$lme)

#Inspect model object
mod31

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod31$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod31$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear ##########################################################################
mod31a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(SST_May_July,k=3),
              data = dat)
#Model summaries
summary(mod31a)

#Inspect model object
mod31a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod31a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod31$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2 + SE wind ################################################
###################################################### linear ######################################################################
mod32 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + SE.wind,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod32)
summary(mod32$gam)
summary(mod32$lme)

#Inspect model object
mod32

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod32$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod32$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear ######################################################################
mod32a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(SE.wind,k=3),
              data = dat)
#Model summaries
summary(mod32a)


#Inspect model object
mod32a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod32a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod32$lme, resid=T, pch=19, rug=F, se=F, pages=1)

############################################# wrap up #########################################################################################################
#MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32)

#######################################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + multiple environmental factors ##################################

###################################################### Combine FHS lag 2, AO_RA3 + PDO_RA3 #####################################################
###################################################### linear #################################################################################
mod33 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA3 + AO_RA3,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod33)
summary(mod33$gam)
summary(mod33$lme)

#Inspect model object
mod33

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod33$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod33$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### nonlinear #################################################################################
mod33a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2,k=3) + s(PDO_RA3,k=3) + s(AO_RA3,k=3),
              data = dat, correlation=corAR1())
#Model summaries
summary(mod33a)

#Inspect model object
mod33a

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod33a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod33$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod34 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA2 + AO_RA2,
              data = dat, correlation=corAR1())

#Model summaries
summary(mod34)
summary(mod34$gam)
summary(mod34$lme)

#Inspect model object
mod34

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod34$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod34$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod34)


###################################################### nonlinear #####################################################
mod34a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2, k=3) + s(PDO_RA2, k = 4) + s(AO_RA2,k=3),
              data = dat)

#Model summaries
summary(mod34a)

#Inspect model object
mod34a

#plot,
dev.new()
par(mfrow=c(2,2))

plot(mod34a, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod34$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod34a)
############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34)

MuMIn::AICc(mod1a,mod2a,mod3a,mod4a, mod5a,mod6a,mod7a,mod8a,mod9a,mod10a,mod11a,mod12a,mod13a,mod14a,mod15a,mod16a,mod17a,mod18a,mod19a,mod20a,mod21a,mod22a,mod23a,
            mod24a,mod25a,mod26a,mod27a,mod28a,mod29a,mod30a,mod31a,mod32a,mod33a,mod34a)


AICc_linear<-MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34)

AICc_nonlinear<-MuMIn::AICc(mod1a,mod2a,mod3a,mod4a, mod5a,mod6a,mod7a,mod8a,mod9a,mod10a,mod11a,mod12a,mod13a,mod14a,mod15a,mod16a,mod17a,mod18a,mod19a,mod20a,
                            mod21a,mod22a,mod23a,mod24a,mod25a,mod26a,mod27a,mod28a,mod29a,mod30a,mod31a,mod32a,mod33a,mod34a)

write.csv(AICc_linear, "output/Mixed_GAM_linear_model_AICc_values_SC3_SC4_k3_reducedto_match_Lag3_Lag4.csv")
#write.csv(AICc_nonlinear, "output/Mixed_GAM_linear_model_AICc_values_SC3_SC4_nonlinear_effects.csv")

################################################################################################################################################
###################################################### BELOW ARE EXPLORATORY MODELS NOT INTENDED FOR REPORTING/PUBLISHING ######################
###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod_test <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA2 + AO_RA2+ SST_May_July,
                 data = dat, correlation=corAR1())
#Model summaries
summary(mod_test)
summary(mod_test$gam)
summary(mod_test$lme)

#Inspect model object
mod_test

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod_test$gam, resid=T, pch=19, rug=F, se=F, pages=1)
#plot(mod34$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod_test)

###################################################### Combine FHS lag 2 + PDO_RA2 + Pcod lag 1#####################################################
###################################################### linear ##########################################################################

mod35 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + FHS_lag2 + PDO_RA2+Pcod_lag1,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod35)
summary(mod35$gam)
summary(mod35$lme)

#Inspect model object
mod35

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod35$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod35)
#plot(mod29$lme, resid=T, pch=19, rug=F, se=F, pages=1)
################################ non-linear ######################################################
mod35a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2, k=3) + s(PDO_RA2,k=3) + s(Pcod_lag1,k=3),
              data = dat)

summary(mod35a)
summary(mod35a)$dev.expl

dev.new()
par(mfrow=c(2,2))

plot(mod35a, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod35a)

################################ nonlinear ######################################################
mod35b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(FHS_lag2, k=3) + s(PDO_RA2,k=3) +s(Pcod_lag1, k=3),
               data = dat, correlation=corAR1())
#Model summaries
summary(mod35b)
summary(mod35b$gam)
summary(mod35b$lme)

#Inspect model object
mod35b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod35b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod35b)

MuMIn::AICc(mod35,mod35a,mod35b)

################################################ SWAP PCOD FOR FHS ###########################################################
###################################################### linear ##########################################################################

mod36 <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + PDO_RA2 + Pcod_lag1,
              data = dat, correlation=corAR1())
#Model summaries
summary(mod36)
summary(mod36$gam)
summary(mod36$lme)

#Inspect model object
mod36

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod36$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod36)
################################ non-linear ######################################################
mod36a <- gam(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PDO_RA2,k=3) + s(Pcod_lag1,k=3),
              data = dat)

summary(mod36a)
summary(mod36a)$dev.expl

dev.new()
par(mfrow=c(2,2))

plot(mod36a, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod36a)

################################ nonlinear ######################################################
mod36b <- gamm(SC3_SC4_logRS ~ s(SC3_SC4_ReproductiveFemales, k=3) + s(PDO_RA2,k=3) +s(Pcod_lag1, k=3),
               data = dat, correlation=corAR1())
#Model summaries
summary(mod36b)
summary(mod36b$gam)
summary(mod36b$lme)

#Inspect model object
mod35b

#plot,
dev.new()
par(mfrow=c(2,1))

plot(mod36b$gam, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod36b)

MuMIn::AICc(mod36,mod36a,mod36b)
MuMIn::AICc(mod35,mod35a,mod35b)
