
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
#?nlstools)
library(MuMIn)
library(tidyverse)
library(corrplot)
library(voxel)

########## Import data and define variables####

Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)

E166Recruits<-read.csv("data/cb_e166_pop_juvenile.csv")
E166Spawners<-read.csv("data/cb_e166_pop_sc_cs.csv")
W166Recruits<-read.csv("data/cb_w166_pop_juvenile.csv")
W166Spawners<-read.csv("data/cb_w166_pop_sc_cs.csv")


#########################################################################################################
########################### Create recruitment time series for analysis #################################
#########################################################################################################
names(Recruits)

########################### EBS ####################################################################################

rec_30to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO50 + Recruits$NUM_FEMALE_30TO50),(E166Recruits$NUM_MALE_30TO50 + E166Recruits$NUM_FEMALE_30TO50),(W166Recruits$NUM_MALE_30TO50 + W166Recruits$NUM_FEMALE_30TO50))))
colnames(rec_30to50)<-c("Year", "EBS_Abun_30to50","E166_Abun_30to50","W166_Abun_30to50")


#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################
names(Spawners)

sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")

e166_sp_sc3<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3))))
colnames(e166_sp_sc3)<-c("Year", "E166_sc3fem_abun")

w166_sp_sc3<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3))))
colnames(w166_sp_sc3)<-c("Year", "W166_sc3fem_abun")

Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
colnames(Spawners_SC3)<-c("Year", "EBS_SC3","E166_SC3","W166_SC3")

########################################################################################################################
########################################## ACF analysis of series by stanza ############################################
########################################################################################################################

########################################################################################################################
########################################## ACF for juvenile indices ####################################################
rec_30to50$Year[4:45] #check years
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

R<-rec_30to50$EBS_Abun_30to50[4:45]


########################################## Divide into two stanzas based on AICc)#######################################
R1<-R[1:21]  #1978-1998
R2<-R[22:42] #1999-2019



acf(R1,main="Juvenile index, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R2,main="Juvenile index, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

########################################################################################################################
########################################## Now do ACF for female indices ###############################################
########################################################################################################################

########################################################################################################################
########################################## Shell condition 3 ###########################################################
Spawners_SC3$Year
S<-Spawners_SC3$EBS_SC3[4:45]

dev.new()
acf(S,main="SC3 females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title

S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019



acf(S1,main="SC3 females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="SC3 females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)



########################################################################################################################
###################################### Plot females using customized graphic ###########################################
y.range<-range(0,4e+08)
par(mfrow=c(2,1),cex.main=1.5, cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels

names(Spawners_SC3)
EBS.fem<-Spawners_SC3$EBS_SC3

plot(EBS.fem,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)


axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4,labels=c("0","100","200","300","400"),tck=0.03)

box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main=" Reproductive female abundance by year")



#######################################################################################################################
###################################### Plot juveniles using customized graphic  #######################################
y.range<-range(0,3e+08)

plot(EBS.juv,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)

#abline(h=0)

axis(1, at=1:45,lab=Recruits$SURVEY_YEAR,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)

box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main=" Juvenile abundance by year")

#######################################################################################################################
########################################### STOCK RECRUIT MODELS  #####################################################
#######################################################################################################################


#######################################################################################################################
########################################### SHELL CONDITION 3 FEMALES  ################################################
#######################################################################################################################

########################################### Lag=2 #####################################################################

########################################### Recreate data series for actual analysis ##################################
par(mfrow=c(3,1))
n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_SC3$EBS_SC3[4:n]
S
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_SC3$Year[4:n]
recYear
R<-R[-7]
S<-S[-5]
recYear<-recYear[-7]
spYear<-spYear[-5]
recYear
spYear



xi<-S
xyear<-spYear
yi<-R
yyear<-recYear
lag <- 2				
n <- length(yi)-1			########Note change from xi and addtion of -1 to allow for different length in yi as compared to xi
xi.k <- xi[1:(n-lag)]       # Select reproductive female estimates values
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
S<-xi.k 
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k


########################################## Now do ACF for log(R/S) ######################################################

logR<-log(R)

dev.new()
par(mfrow=c(2,1))

logR_1<-logR[1:20]  #1978-1998
logR_2<-logR[21:38] #1999-2019



acf(logR_1,main="Log(R), lag 2-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logR_2,main="Log(R), lag 2-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

####################################Fit GAMM models with 1st and 2nd order autocorrelation#############################


mod1<-gamm(log(R)~s(S),correlation=corAR1(), na.action=na.omit)

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

# Calculate AICc value
MuMIn::AICc(mod1)      #
AIC_lag2<-MuMIn::AICc(mod1)


#######################################################################################################################
######################################### Lag=3 #######################################################################
#par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_SC3$EBS_SC3[4:n]

R
S
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_SC3$Year[4:n]
release.year<-c(1978:20016)

R<-R[-8]
S<-S[-5]
recYear<-recYear[-8]
spYear<-spYear[-5]


xi<-S
xyear<-spYear
yi<-R
yyear<-recYear
lag <- 3
n <- length(yi)			##########note change from xi as previously employed
xi.k <- xi[1:(n-lag)]       # Select reproductive female ests
S<-xi.k
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k 
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k

########################################## Now do ACF for log(R/S) ######################################################

logR<-log(R)

dev.new()
par(mfrow=c(2,1))

logR_1<-logR[1:20]  #1978-1998
logR_2<-logR[21:38] #1999-2019



acf(logR_1,main="Log(R), lag 3-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logR_2,main="Log(R), lag 3-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

#################################Fit GAMM model with 1st order autocorrelation###########################################

mod2<-gamm(log(R)~s(S),correlation=corAR1(), na.action=na.omit)

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


# Calculate AICc value

MuMIn::AICc(mod2)      #
AIC_lag3<-MuMIn::AICc(mod2)



####################################################################################################################
################################################# Lag=4 ############################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_SC3$EBS_SC3[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_SC3$Year[4:n]

xi<-S
xyear<-spYear[-5]
yi<-R
yyear<-recYear[-9]
lag <- 4
n <- length(yi)			##########note change from xi as previously employed
xi.k <- xi[1:(n-lag)]       # Select reproductive female ests
S<-xi.k
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k 
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k

########################################## Now do ACF for log(R/S) #################################################

logR<-log(R)

dev.new()
par(mfrow=c(2,1))

logR_1<-logR[1:20]  #1978-1998
logR_2<-logR[21:38] #1999-2019



acf(logR_1,main="Log(R), lag 4-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logR_2,main="Log(R), lag 4-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
####################################################################################################################
############################################## Fit GAMM model with 1st order autocorrelation #######################

mod3<-gamm(log(R)~s(S),correlation=corAR1(), na.action=na.omit)

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


# Calculate AICc value
MuMIn::AICc(mod3)      #
AIC_lag4<-MuMIn::AICc(mod3)

########################################### AICc values for all models ############################################
AIC_lag2
AIC_lag3
AIC_lag4
