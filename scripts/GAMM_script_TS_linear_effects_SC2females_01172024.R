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

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis_reducedTime_updated.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed, and is reduced to match shortest dataseries (FHS)

head(dat)

#
###################################################### Run models#########################################################

###################################################### S-R effect only ##################################################
mod1 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4),
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



###################################################### Add ovigerous female opilio for competitive effect amongst juveniles ###

mod2 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + Ovig_female_CO,
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



############################################################################################################################
################################################## Predation models ########################################################

###################################################### Add Pcod lag1 only ##################################################
mod3 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + Pcod_lag1,
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


###################################################### Add PCod_RA2 only ##################################################
mod4 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + PCod_RA2,
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


###################################################### Add PCod_RA3 only ##################################################
mod5 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + PCod_RA3,
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

###################################################### Add FHS lag 2 only ##################################################
mod6 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2,
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


###################################################### Add FHS RA2 only ###################################################
mod7 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_RA2,
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

##################################################### Combine Pcod_lag1 and FHS lag 2 #####################################

mod8 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + Pcod_lag1 + FHS_lag2,
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

##################################################### Combine PCod_RA2 and FHS_RA2 #####################################

mod9 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + PCod_RA2 + FHS_RA2,
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

########################################################################################################################
##################################################### Environmental factors ############################################

##################################################### NBT 3 yr rolling average #########################################

mod10 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + NBT_3RA,
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

##################################################### NBT 3 yr min temperature #########################################

mod11 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + NBT_3yr_minTemp,
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



##################################################### AO_RA2 #########################################

mod12 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + AO_RA2,
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


##################################################### AO_RA3 #########################################

mod13 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + AO_RA3,
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



##################################################### PDO_RA2 #########################################

mod14 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + PDO_RA2,
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



##################################################### PDO_RA3 #########################################

mod15 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + PDO_RA3,
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

##################################################### SST_May_July #########################################

mod16 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + SST_May_July,
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

##################################################### SE.wind #########################################

mod17 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + SE.wind,
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


##################################################### NW.wind #########################################

mod18 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + NW.wind,
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


############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18)


###############################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + competition/environmental factor######################


###################################################### Combine FHS lag 2 and ovigerous female opilio  #########################################################
mod19 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO,
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

###################################################### Combine FHS lag 2, ovigerous female opilio + PDO_RA3 #####################################################
mod20 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + PDO_RA3,
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

###################################################### Combine FHS lag 2, ovigerous female opilio + AO_RA3 #####################################################
mod21 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + AO_RA3,
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
###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod22 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + PDO_RA2,
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

###################################################### Combine FHS lag 2, AO_RA2  #####################################################
mod23 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + AO_RA2,
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
###################################################### Combine FHS lag 2, ovigerous female opilio + SST_May_July ################################################
mod24 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + SST_May_July,
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

###################################################### Combine FHS lag 2, ovigerous female opilio + SE wind ################################################
mod25 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + SE.wind,
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

###################################################### Combine FHS lag 2, ovigerous female opilio + NE wind ################################################
mod26 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + Ovig_female_CO + NW.wind,
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

############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26)


###############################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + environmental factor##################################

###################################################### Combine FHS lag 2 + PDO_RA3 #####################################################
mod27 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3,
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

###################################################### Combine FHS lag 2 + AO_RA3 #####################################################
mod28 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + AO_RA3,
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

###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod29 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA2,
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
#plot(mod29$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###################################################### Combine FHS lag 2, AO_RA2  #####################################################
mod30 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + AO_RA2,
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

###################################################### Combine FHS lag 2 + SST_May_July ################################################
mod31 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + SST_May_July,
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

###################################################### Combine FHS lag 2 + SE wind ################################################
mod32 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + SE.wind,
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

############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32)

#######################################################################################################################################################################
##################################################### Combine multiple variables in combined pred-prey + multiple environmental factors ##################################

###################################################### Combine FHS lag 2, AO_RA3 + PDO_RA3 #####################################################
mod33 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3 + AO_RA3,
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

###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod34 <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA2 + AO_RA2,
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



############################################# wrap up #########################################################################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34)
AICc_linear<-MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34)

write.csv(AICc_linear, "output/Mixed_GAM_linear_model_AICc_values_SC2.csv")


###################################################### Combine FHS lag 2, AO_RA2 + PDO_RA2 #####################################################
mod_test <- gamm(SC2logRS ~ s(SC2ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA2 + AO_RA2+ SST_May_July,
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
