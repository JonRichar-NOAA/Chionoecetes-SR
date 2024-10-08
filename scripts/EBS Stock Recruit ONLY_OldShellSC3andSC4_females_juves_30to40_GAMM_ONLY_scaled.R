##################### Clear environment ############################################################
rm(list=ls())
##################### Load packages  ###############################################################

##### Note for Valkyrie: above should read: setwd('C:/Users/Jon Richar/Desktop/Project/Datasets for Analysis')
##### Note: for Paladin, above should read: setwd('C:/Documents and Settings/Jon/Desktop/Project/Datasets for analysis')
 
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
library(tidyverse)
#?nlstools
########## Citation for packages

citation("mgcv")
########## Import data and define variables####

Recruits<-read.csv("./data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("./data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)

E166Recruits<-read.csv("./data/cb_e166_pop_juvenile.csv")
E166Spawners<-read.csv("./data/cb_e166_pop_sc_cs.csv")
W166Recruits<-read.csv("./data/cb_w166_pop_juvenile.csv")
W166Spawners<-read.csv("./data/cb_w166_pop_sc_cs.csv")
#######################################################################################################
####################################### EDA.norm function #############################################

eda.norm <- function(x, ...)
{
# Examine distribution of data and check against normal distribution
# x is a vector of data. Additional graphics parameters can be supplied
# The function creates a histogram with an empirical density estimate, 
# a boxplot, a normal q-q plot, and a plot of the empirical cumulative
# density function with the corresponding normal cdf. 
# In addition, the function returns the results from the 
# Shapiro-Wilks test of normality
#
# Written by Franz Mueter. Last modified February 24, 2006
#

par(mfrow=c(2,2))
if(sum(is.na(x)) > 0)
warning("NA's were removed before plotting")

x <- x[!is.na(x)]
hist(x, main = "Histogram and non-\nparametric density estimate", prob = T)
iqd <- summary(x)[5] - summary(x)[2]
lines(density(x, width = 2 * iqd))
boxplot(x, main = "Boxplot", ...)
qqnorm(x)
qqline(x)
plot.ecdf(x, main="Empirical and normal cdf")
LIM <- par("usr")
y <- seq(LIM[1],LIM[2],length=100)
lines(y, pnorm(y, mean(x), sqrt(var(x))))
shapiro.test(x)
}



#########################################################################################################
########################### Create recruitment time series for analysis #################################
#########################################################################################################
names(Recruits)

########################### EBS ####################################################################################
rec_30to40<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO40 + Recruits$NUM_FEMALE_30TO40),(E166Recruits$NUM_MALE_30TO40 + E166Recruits$NUM_FEMALE_30TO40),(W166Recruits$NUM_MALE_30TO40 + W166Recruits$NUM_FEMALE_30TO40))))
rec_40to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_40TO50 + Recruits$NUM_FEMALE_40TO50),(E166Recruits$NUM_MALE_40TO50 + E166Recruits$NUM_FEMALE_40TO50),(W166Recruits$NUM_MALE_40TO50 + W166Recruits$NUM_FEMALE_40TO50))))
rec_50to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_50TO60 + Recruits$NUM_FEMALE_50TO60),(E166Recruits$NUM_MALE_50TO60 + E166Recruits$NUM_FEMALE_50TO60),(W166Recruits$NUM_MALE_50TO60 + W166Recruits$NUM_FEMALE_50TO60))))

rec_30to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO50 + Recruits$NUM_FEMALE_30TO50),(E166Recruits$NUM_MALE_30TO50 + E166Recruits$NUM_FEMALE_30TO50),(W166Recruits$NUM_MALE_30TO50 + W166Recruits$NUM_FEMALE_30TO50))))
rec_30to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO60 + Recruits$NUM_FEMALE_30TO60),(E166Recruits$NUM_MALE_30TO50 + E166Recruits$NUM_FEMALE_30TO50),(W166Recruits$NUM_MALE_30TO50 + W166Recruits$NUM_FEMALE_30TO50))))

rec_40to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_40TO50 + Recruits$NUM_FEMALE_40TO50+Recruits$NUM_MALE_50TO60 + Recruits$NUM_FEMALE_50TO60),
(E166Recruits$NUM_MALE_40TO50 + E166Recruits$NUM_FEMALE_40TO50 + E166Recruits$NUM_MALE_50TO60 + E166Recruits$NUM_FEMALE_50TO60),
(W166Recruits$NUM_MALE_40TO50 + W166Recruits$NUM_FEMALE_40TO50+W166Recruits$NUM_MALE_50TO60 + W166Recruits$NUM_FEMALE_50TO60))))

colnames(rec_30to40)<-c("Year", "EBS_Abun_30to40","E166_Abun_30to40","W166_Abun_30to40")
colnames(rec_40to50)<-c("Year", "EBS_Abun_40to50","E166_Abun_40to50","W166_Abun_40to50")
colnames(rec_50to60)<-c("Year", "EBS_Abun_50to60","E166_Abun_50to60","W166_Abun_50to60")

colnames(rec_30to50)<-c("Year", "EBS_Abun_30to50","E166_Abun_30to50","W166_Abun_30to50")
colnames(rec_30to60)<-c("Year", "EBS_Abun_30to60","E166_Abun_30to60","W166_Abun_30to60")
colnames(rec_40to60)<-c("Year", "EBS_Abun_40to60","E166_Abun_40to60","W166_Abun_40to60")

rec_40to60
cor(rec_30to50$EBS_Abun_30to50,rec_30to40$EBS_Abun_30to40)
test<-lm(rec_30to40$EBS_Abun_30to40~rec_30to50$EBS_Abun_30to50)
summary(test)
#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################
names(Spawners)

sp_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_OVIGEROUS))))
sp_sc2<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC2))))
sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))
sp_sc4<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC4))))
sp_os<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4))))

colnames(sp_ovig)<-c("Year", "EBS_ovigfem_abun")
colnames(sp_sc2)<-c("Year", "EBS_sc2fem_abun")
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")
colnames(sp_sc4)<-c("Year", "EBS_sc4fem_abun")
colnames(sp_os)<-c("Year", "EBS_OSfem_abun")

e166_sp_ovig<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_OVIGEROUS))))
e166_sp_sc2<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC2))))
e166_sp_sc3<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3))))
e166_sp_sc4<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC4))))
e166_sp_os<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4))))

colnames(e166_sp_ovig)<-c("Year", "E166_ovigfem_abun")
colnames(e166_sp_sc2)<-c("Year", "E166_sc2fem_abun")
colnames(e166_sp_sc3)<-c("Year", "E166_sc3fem_abun")
colnames(e166_sp_sc4)<-c("Year", "E166_sc4fem_abun")
colnames(e166_sp_os)<-c("Year", "E166_OSfem_abun")

w166_sp_ovig<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_OVIGEROUS))))
w166_sp_sc2<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC2))))
w166_sp_sc3<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3))))
w166_sp_sc4<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC4))))
w166_sp_os<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4))))

colnames(w166_sp_ovig)<-c("Year", "W166_ovigfem_abun")
colnames(w166_sp_sc2)<-c("Year", "W166_sc2fem_abun")
colnames(w166_sp_sc3)<-c("Year", "W166_sc3fem_abun")
colnames(w166_sp_sc4)<-c("Year", "W166_sc4fem_abun")
colnames(w166_sp_os)<-c("Year", "W166_OSfem_abun")

Spawners_SC0<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC0,E166Spawners$NUM_FEMALE_SC0,W166Spawners$NUM_FEMALE_SC0)))
Spawners_SC1<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC1,E166Spawners$NUM_FEMALE_SC1,W166Spawners$NUM_FEMALE_SC1)))
Spawners_SC2<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2)))
Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
Spawners_SC4<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC4,E166Spawners$NUM_FEMALE_SC4,W166Spawners$NUM_FEMALE_SC4)))
Spawners_SC5<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC5,E166Spawners$NUM_FEMALE_SC5,W166Spawners$NUM_FEMALE_SC5)))
Spawners_OS<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4),(E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4),(W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4))))

Spawners_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_OVIGEROUS,E166Spawners$NUM_FEMALE_OVIGEROUS,W166Spawners$NUM_FEMALE_OVIGEROUS)))


colnames(Spawners_SC2)<-c("Year","EBS_SC2","E166_SC2","W166_SC2")
Spawners_SC2

colnames(Spawners_SC3)<-c("Year","EBS_SC3","E166_SC3","W166_SC3")
Spawners_SC3

colnames(Spawners_ovig)<-c("Year","EBS_ovig","E166_ovig","W166_ovig")
Spawners_ovig

colnames(Spawners_OS)<-c("Year","EBS_os","E166_os","W166_os")
Spawners_OS

Spawners_analysis<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,
                                                 Spawners$NUM_FEMALE_SC0,E166Spawners$NUM_FEMALE_SC0,W166Spawners$NUM_FEMALE_SC0,
                                                 Spawners$NUM_FEMALE_SC1,E166Spawners$NUM_FEMALE_SC1,W166Spawners$NUM_FEMALE_SC1,
                                                 Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2,
                                                 Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3,
                                                 Spawners$NUM_FEMALE_SC4,E166Spawners$NUM_FEMALE_SC4,W166Spawners$NUM_FEMALE_SC4,
                                                 Spawners$NUM_FEMALE_SC5,E166Spawners$NUM_FEMALE_SC5,W166Spawners$NUM_FEMALE_SC5,
                                                 (Spawners$NUM_FEMALE_SC0+Spawners$NUM_FEMALE_SC1+Spawners$NUM_FEMALE_SC2+Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4+Spawners$NUM_FEMALE_SC5),
                                                 (E166Spawners$NUM_FEMALE_SC0+E166Spawners$NUM_FEMALE_SC1+E166Spawners$NUM_FEMALE_SC2+E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4+E166Spawners$NUM_FEMALE_SC5),
                                                 (W166Spawners$NUM_FEMALE_SC0+W166Spawners$NUM_FEMALE_SC1+W166Spawners$NUM_FEMALE_SC2+W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4+W166Spawners$NUM_FEMALE_SC5),
                                                 Spawners$NUM_FEMALE_OVIGEROUS,E166Spawners$NUM_FEMALE_OVIGEROUS,W166Spawners$NUM_FEMALE_OVIGEROUS,
                                                 Spawners_OS$EBS_os,Spawners_OS$E166_os, Spawners_OS$W166_os)))

colnames(Spawners_analysis)<-c("Year","EBS_SC0","E166_SC0","W166_SC0",
                               "EBS_SC1","E166_SC1","W166_SC1",
                               "EBS_SC2","E166_SC2","W166_SC2",
                               "EBS_SC3","E166_SC3","W166_SC3",
                               "EBS_SC4","E166_SC4","W166_SC4",
                               "EBS_SC5","E166_SC5","W166_SC5",
                               "EBS_SC_TOTAL","E166_SC_TOTAL","W166_SC_TOTAL",
                               "EBS_ovig","E166_ovig","W166_ovig",
                               "EBS_os","E166_os","W166_os")

#Spawners_analysis
#write.csv(Spawners_analysis,"./data/Female_Tanner_Crab_Series_for_analysis.csv")
#write.csv(rec_30to40,"./data/JuvenileTannerAbun_30to40mmBin.csv")


#########################################################################################################################
########################################### Lag=2 #######################################################################
par(mfrow=c(2,2))

#########################################Define data for lag=2###########################################################
par(mfrow=c(2,2))

rec_30to40
R<-rec_30to40$EBS_Abun_30to40[4:45]
S<-Spawners_OS$EBS_os[4:45]
log.R<-log(R)
log.S<-log(S)


xi<-S
xyear<-Spawners_OS$Year[4:45]
yi<-R
yyear<-rec_30to40$Year[4:45]
lag <- 2				
n <- length(yi)-1			########Note change from xi and addtion of -1 to allow for different length in yi as compared to xi
xi.k <- xi[1:(n-lag)]       # Select reproductive female estimates values
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
S<-xi.k 
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k/0.3
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 2 juvenile recruitment", pch=16,col=4)

xi.k
max(xi.k)

logRS<-log(R/S)
dat<-as.data.frame(cbind(as.matrix(R),as.matrix(S)))
##############################################EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)
################################## Fit GAMM mode #######################################
mod1 <- gamm(log(R/S) ~ s(S, k=3),
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
#plot(mod1$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod1)
lag2aicc<-MuMIn::AICc(mod1)


############################### create simulated data ####
#set.seed(100) 
newdat <- data.frame(S=seq(from=min(S), to=max(S), length.out=100)) 
lag<-as.matrix(rep("Lag 2 yr",times=nrow(as.matrix(newdat))))
predictions = predict(mod1$gam, newdata=newdat, se.fit = TRUE)

# Consolidating new data and predictions
dat2 = as.data.frame(cbind(lag, newdat, predictions))
dat2<- within(dat2, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

view(dat2)

lag_point<-as.matrix(rep("Lag 2 yr",times=nrow(as.matrix(R))))
logRS_2<-log(R/S)
p_dat2<-as.data.frame(cbind(lag_point,S,logRS_2))


colnames(p_dat2)<-c("lag","S","logRS")
view(p_dat2)
#########################################################################################################################
############################################### Drop data point #5  #####################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25)
R2<-R[-5]
R2
R
S2<-S[-5]

dat<-as.data.frame(cbind(as.matrix(R2),as.matrix(S2)))
############################################## GAMM model ###################################################
mod2 <- gamm(log(R2/S2) ~ s(S2, k=3),
              correlation=corAR1())
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

########################################################################################################
######################################### Lag=3 ########################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
R<-rec_30to40$EBS_Abun_30to40[4:45]
S<-Spawners_OS$EBS_os[4:45]
#release.year<-c(1978:2005)

xi<-S
xi
xyear<-Spawners_ovig$Year[4:45]
xyear
yi<-R
yyear<-rec_30to40$Year[4:45]
lag <- 3
n <- length(yi)-1			########Note change from xi and addtion of -1 to allow for different length in yi as compared to xi
xi.k <- xi[1:(n-lag)]       # Select reproductive female estimates values
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
S<-xi.k 
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later

R<-yi.k/0.3
R
yi.k
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 2 juvenile recruitment", pch=16,col=4)

cor(R,S)
plot(R~S)
plot((R/S)~S)
xi.k
max(xi.k)


logRS<-log(R/S)
dat<-as.data.frame(cbind(as.matrix(R),as.matrix(S)))
##############################################EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
dev.new()
eda.norm(log.RS)


################################## Run GAMM model ############################################
mod3 <- gamm(log(R/S) ~ s(S, k=3),
              correlation=corAR1())
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
#plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)

MuMIn::AICc(mod3)
lag3aicc<-MuMIn::AICc(mod3)

################################ create simulated data ##################################################
#set.seed(100) 
newdat <- data.frame(S=seq(from=min(S), to=max(S), length.out=100)) 
lag<-as.matrix(rep("Lag 3 yr",times=nrow(as.matrix(newdat))))
predictions = predict(mod3$gam, newdata=newdat, se.fit = TRUE)

# Consolidating new data and predictions
dat3 = as.data.frame(cbind(lag, newdat, predictions))
dat3<- within(dat3, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

lag_point<-as.matrix(rep("Lag 3 yr",times=nrow(as.matrix(R))))
logRS_3<-log(R/S)
p_dat3<-as.data.frame(cbind(lag_point,S,logRS_3))


colnames(p_dat3)<-c("lag","S","logRS")
view(p_dat3)
#########################################################################################################################
############################################### Drop data point #5  #####################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25)
R2<-R[-5]
R2
R
S2<-S[-5]

############################################## GAMM model ###################################################
mod4 <- gamm(log(R2/S2) ~ s(S2, k=3),
             correlation=corAR1())
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
#plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)
MuMIn::AICc(mod4)



##############################################################################################################################
################################################# Lag=4 ######################################################################
par(mfrow=c(1,1))
R<-rec_30to40$EBS_Abun_30to40[4:45]
S<-Spawners_OS$EBS_os[4:45]

xi<-S
xyear<-Spawners_OS$Year[4:45]
yi<-R
yyear<-rec_30to40$Year[4:45]
lag <- 4
n <- length(yi)			##########note change from xi as previously employed
xi.k <- xi[1:(n-lag)]       # Select reproductive female ests
S<-xi.k
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k/0.3 
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 4 juvenile recruitment",pch=16,col=4)

logRS<-log(R/S)
dat<-as.data.frame(cbind(as.matrix(R),as.matrix(S)))
############################################## EDA.norm analyses of log(R/S)#############################
dev.new()
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)

################################## Run GAMM model ############################################
mod5 <- gamm(log(R/S) ~ s(S, k=3),
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
lag4aicc<-MuMIn::AICc(mod5)

# create simulated data ####
#set.seed(100) 

newdat <- data.frame(S=seq(from=min(S), to=max(S), length.out=100)) 
lag<-as.matrix(rep("Lag 4 yr",times=nrow(as.matrix(newdat))))
predictions = predict(mod5$gam, newdata=newdat, se.fit = TRUE)

# Consolidating new data and predictions
dat4 = as.data.frame(cbind(lag, newdat, predictions))
dat4<- within(dat4, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

lag_point<-as.matrix(rep("Lag 4 yr",times=nrow(as.matrix(R))))
logRS_4<-log(R/S)
p_dat4<-as.data.frame(cbind(lag_point,S,logRS_4))


colnames(p_dat4)<-c("lag","S","logRS")
view(p_dat4)

#########################################################################################################################
############################################### Drop data point #5  #####################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25)
R2<-R[-5]
R2
R
S2<-S[-5]

############################################## GAMM model ###################################################
mod6 <- gamm(log(R2/S2) ~ s(S2, k=3),
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

############################################# AICc for SR models ##########################################

lag2aicc
lag3aicc
lag4aicc
MuMIn::AICc(mod1,mod3,mod5)

dat<-as.data.frame(rbind(dat2,dat3,dat4))

view(dat2)

p_dat<-as.data.frame(rbind(p_dat2,p_dat3,p_dat4))
str(p_dat2)
str(p_dat)

view(p_dat)
#################### TEST PLOT #################################################
view(dat)

view(logRS)

ggplot(dat, aes(S,  fit)) +
  geom_smooth() +
  geom_point(data=p_dat,aes(x=as.numeric(S), y= as.numeric(logRS))) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=S, fill = "band"), alpha = 0.3)+
  facet_wrap(~factor(lag, c("Lag 2 yr","Lag 3 yr","Lag 4 yr")), scales = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank()) +
  theme_bw()+
  theme(legend.position="none")+
  ylab("Ln(R/S)") 
## NOTE: GEOM_POINT currently adds all data points to each of the plots, rather than only the appropriate data set
##############################################################################

p<-ggplot(dat, aes(S,  fit)) +
  geom_smooth() +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=S, fill = "band"), alpha = 0.3)+
  facet_wrap(~factor(lag, c("Lag 2 yr","Lag 3 yr","Lag 4 yr")), scales = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank()) +
  ylab("Ln(R/S)") 
p
