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

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

#
###################################################### Run models#########################################################

mod1 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3)+ s(PDO_RA3, k=4)+ s(AO_RA3, k=4),
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

#Now run with era interaction

#########################################################################################################################
mod2 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2*Era_AICc + PDO_RA3*Era_AICc + AO_RA3*Era_AICc,
            data = dat, correlation=corAR1())        # Error message: Error in MEestimate(lmeSt, grps) : 
                                                     # Singularity in backsolve at level 0, block 1
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

########################################################################################################################
mod3 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(FHS_lag2, k=3)+ s(PDO_RA3, k=4, by =Era_AICc),
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

#########################################################################################################################
mod4 <- gamm(logRS ~ s(ReproductiveFemales, k=4) +  s(AO_RA3, k=4, by = Era_AICc) + s(FHS_lag2, k=3),
             data = dat, correlation=corAR1())        # Error message: Error in MEestimate(lmeSt, grps) : Singularity in backsolve at level 0, block 1
                                                      # This error is likely due to correlation between PDO and AO -- see correlation analysis at end of code
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

MuMIn::AICc(mod1, mod2,mod3,mod4) 

n<-length(dat)
dat
cor.dat<-cbind(dat$ReproductiveFemales[4:n],dat$FHS_lag2[4:n],dat$PDO_RA3[4:n],dat$AO_RA3[4:n])
cor(cor.dat) # PDO and AO strongly negatively correlated

###################################################################################################################################################################################
###################################################################################################################################################################################

########################################################################################################################
mod5 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(AO_RA3, k=4, by = Era_AICc) + s(PDO_RA3),
             data = dat, correlation=corAR1())         

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

########################################################################################################################
mod6 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + s(AO_RA3, k=4, by = Era_AICc) ,
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

#########################################################################################################################
mod7 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3 + s(AO_RA3, k=4, by = Era_AICc),
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

#########################################################################################################################
mod8 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3 + AO_RA3:Era_AICc,
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

#########################################################################################################################
mod9 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3:Era_AICc + AO_RA3,
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

#########################################################################################################################
mod10 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2 + PDO_RA3:Era_AICc + AO_RA3:Era_AICc,
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

plot(mod10)
plot(mod10$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod10$lme, resid=T, pch=19, rug=F, se=F, pages=1)

###### Neither of following codes currently works ###################
plotGAMM(mod10, smooth.cov, groupCovs = NULL, orderedAsFactor = F,
         rawOrFitted = F, plotCI = T, grouping = NULL)

r1<-resid(mod10)
r1
#plot(1r,pch=16)
r.fit1<-loess(r1~dat$releaseyear)
plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals")
lines(dat$releaseyear,fitted(r.fit1),col=4)

## examine the residuals for mod1
resid1 <- data.frame(year=dat$releaseyear,
                     resid=residuals(mod1))

#plot(resid1$resid~resid1$year)


ggplot(resid1, aes(year, resid)) +
  geom_line() +
  geom_point()

dwtest(resid(mod1)~1)
acf(resid(mod1)) # not perfect, not terrible!
#########################################################################################################################
mod11 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + FHS_lag2,
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

#########################################################################################################################
mod12 <- gamm(logRS ~ s(ReproductiveFemales, k=4) +s(PDO_RA3, k=4, by = Era_AICc),
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


#########################################################################################################################
mod13 <- gamm(logRS ~ s(ReproductiveFemales, k=4) + PDO_RA3:Era_AICc,
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

############################################# wrap up ##########################################################
MuMIn::AICc(mod1,mod2,mod3,mod4, mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)

