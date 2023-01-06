
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


########################################## Divide into two stanzas based on AICc)##############################################
R1<-R[1:21]  #1978-1998
R2<-R[22:42] #1999-2019



acf(R1,main="Juvenile index, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R2,main="Juvenile index, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

###############################################################################################################################
########################################## Now do ACF for female indices ######################################################
###############################################################################################################################

###############################################################################################################################
########################################## Shell condition 3 #################################################################
Spawners_SC3$Year
S<-Spawners_SC3$EBS_SC3[4:45]

dev.new()
acf(S,main="SC3 females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title

S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019



acf(S1,main="SC3 females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="SC3 females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)



########################################################################################################################

for(i in names(Spawners_SC3)) {		# cycle through all names in data frame
	x <- Spawners_SC3$Year			# convert to numeric variable
	y <- Spawners_SC3[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
}


########################################################################################################################
###################################### Plot females using customized graphic ###########################################
y.range<-range(0,4e+08)
par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
names(Spawners_SC3)
W166.fem<-Spawners_SC3$W166_SC3
EBS.fem<-Spawners_SC3$EBS_SC3
E166.fem<-Spawners_SC3$E166_SC3
plot(EBS.fem)

plot(EBS.fem,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)
lines(W166.fem,type="b",pch=18,col=1)
lines(E166.fem,type="b",pch=22,col=2)

axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4,labels=c("0","100","200","300","400"),tck=0.03)
legend(35,4e+08,c("Total EBS", "Western area", "Eastern area"),cex=1.25,col=c(4,1,2),pch=c(16,18,22))
#abline(h=0)
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
#title(main=" Reproductive female abundance by year")



########################################################################################################################
###################################### Plot juveniles using customized graphic  #######################################
par(mfrow=c(1,1),cex.lab=1.4,cex.axis=1.4,cex=1.5) #configure axis labels
y.range<-range(0,9e+08)

plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)
lines(W166.juv,type="b",pch=18,col=1,lty=2)
lines(E166.juv,type="b",pch=22,col=1,lty=3)
#abline(h=0)

axis(1, at=1:31,lab=Recruits$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("Total EBS", "Western area", "Eastern area"),cex=1.25,col=c(1,1,1),pch=c(2,18,22))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
#title(main=" Juvenile abundance by year")



####################################################################################################################
###################################### Plot EBS juve abundance for creation of recruitment bins ####################
n<-nrow(Recruits)
par(mfrow=c(1,1),cex.lab=1.25,cex.axis=1.25,cex=1.25)
plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)
axis(1, at=1:n,lab=Recruits$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1),pch=c(16,18))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()
mean<-mean(EBS.juv)
abline(h=mean,col=1)
dev<-sd(EBS.juv)
dev
abline(h=mean+.65*dev)
abline(h=mean-.65*dev)
abline(h=0)

title(xlab="Release Year")
title(ylab="Juvenile abundance (Millions)")
#########################################################################################################################
###################################### Plot juv abundance separately ####################################################

Year<-Recruits$Year
plot(R~Year, type="b",ylab="Estimated abundance of EBS juvenile bairdi", xlab="Year", main = "Estimated abundance of EBS juvenile bairdi",pch=16)
abline(lm(R ~ Year), col=4)
fit <- loess(R ~ Year)
lines(x, fitted(fit), col=2)

#########################################################################################################################################
########################################### STOCK RECRUIT MODELS  #######################################################################
#######################################################################################################################################


##############################################################################################################################################
########################################### SHELL CONDITION 3 FEMALES  #######################################################################
##############################################################################################################################################

########################################### Lag=2 ###########################################################################################

########################################### Recreate data series for actual analysis ########################################################
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

logRS<-log(R/S)

dev.new()
par(mfrow=c(2,1))

logRS_1<-logRS[1:20]  #1978-1998
logRS_2<-logRS[21:38] #1999-2019



acf(logRS_1,main="Log(R/S), lag 2-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logRS_2,main="Log(R/S), lag 2-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
####################################Fit gls models with 1st and 2nd order autocorrelation######################
#par(mfrow=c(1,1))
mod1<-gamm(log(R/S)~s(S),correlation=corAR1(), na.action=na.omit)

summary(mod1)

summary(mod1$gam)
summary(mod1$lme)

mod1

plot(mod1$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod1$lme, resid=T, pch=19, rug=F, se=F, pages=1)


#plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main= "Lag 2 S-R model",pch=16)

MuMIn::AICc(mod1)      #
AIC_lag2<-MuMIn::AICc(mod1)


########################################################################################################
######################################### Lag=3 ########################################################
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

logRS<-log(R/S)

dev.new()
par(mfrow=c(2,1))

logRS_1<-logRS[1:20]  #1978-1998
logRS_2<-logRS[21:38] #1999-2019



acf(logRS_1,main="Log(R/S), lag 3-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logRS_2,main="Log(R/S), lag 3-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

#################################Fit gls model with 1st order autocorrelation#######################################################
mod2<-gamm(log(R/S)~s(S),correlation=corAR1(), na.action=na.omit)
summary(mod2)
summary(mod2$gam)
summary(mod2$lme)

mod2

plot(mod2$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod2$lme, resid=T, pch=19, rug=F, se=F, pages=1)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main= "Lag 3 S-R model",pch=16)

MuMIn::AICc(mod2)      #
AIC_lag3<-MuMIn::AICc(mod2)



##############################################################################################################################
################################################# Lag=4 ######################################################################

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

########################################## Now do ACF for log(R/S) ######################################################

logRS<-log(R/S)

dev.new()
par(mfrow=c(2,1))

logRS_1<-logRS[1:20]  #1978-1998
logRS_2<-logRS[21:38] #1999-2019



acf(logRS_1,main="Log(R/S), lag 4-yr, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(logRS_2,main="Log(R/S), lag 4-yr, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
####################################################################################################################################
############################################## Fit GAMM model with 1st order autocorrelation ########################################

mod3<-gamm(log(R/S)~s(S),correlation=corAR1(), na.action=na.omit)

summary(mod3)
summary(mod3$gam)
summary(mod3$lme)

mod3

plot(mod3$gam, resid=T, pch=19, rug=F, se=F, pages=1)
plot(mod3$lme, resid=T, pch=19, rug=F, se=F, pages=1)

plot(mod3, resid=T, pch=19, rug=F, se=F, pages=1)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main= "Lag 4 S-R model",pch=16)

MuMIn::AICc(mod3)      #
AIC_lag4<-MuMIn::AICc(mod3)

########################################### AICc values for all models ######################################################
AIC_lag2
AIC_lag3
AIC_lag4
