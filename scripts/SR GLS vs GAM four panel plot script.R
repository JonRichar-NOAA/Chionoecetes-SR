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


cor.dat <- dat %>%
  select(-era)

cor(cor.dat, use="p")


## fit models to explain recruitment variability -------------

## begin with spawner-recruit models ----------------
## with and without 2000 breakpoint suggested by earlier
mod1 <- gam(logRS ~ s(ReproductiveFemales, k = 4),
            data = dat)

summary(mod1)

dev.new()
par(mfrow=c(2,1),cex.main=1, cex.lab=1,cex.axis=1,cex=1) #configure axis labels

plot(mod1, resid=T, pch=19, rug=F, se=F, ylab="Log(recruits/spawners)",xlab = "SC3 female abundance")


r1<-resid(mod1)
#plot(1r,pch=16)
r.fit1<-loess(r1~dat$releaseyear)
plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals")
lines(dat$releaseyear,fitted(r.fit1),col=4)

######################################### Linear GLS ##################################################

Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)


E166Recruits<-read.csv("data/cb_e166_pop_juvenile.csv")
E166Spawners<-read.csv("data/cb_e166_pop_sc_cs.csv")
W166Recruits<-read.csv("data/cb_w166_pop_juvenile.csv")
W166Spawners<-read.csv("data/cb_w166_pop_sc_cs.csv")
#######################################################################################################



#########################################################################################################
################################### Define baseline Ricker model ########################################
#Ricker <- function(S,a,b){
#a*S*exp(-b*S)
#}
#START <- list(a=1, b=0.01)

#########################################################################################################
############################ Define linearized Ricker model--note log(R/S)= #############################

Ricker <- function(S,a,b){
  a+b*S
}
START<-list(a=1,b=0.1)


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

Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
Spawners_SC3
colnames(Spawners_SC3)<-c("Year", "EBS_SC3","E166_SC3","W166_SC3")
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
release.year<-c(1978:2016)

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


#################################Fit gls model with 1st order autocorrelation#######################################################

Fit4<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit4)
#par(mfrow=c(1,1))
par(mfrow=c(2,2),cex.main=1, cex.lab=1,cex.axis=1,cex=1) #configure axis labels

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="Log(R/S)",main="Lag 3 S-R model - GLS",pch=16)
abline(Fit4)

#plot(Fit3)
coef(Fit4)

r4<-resid(Fit4)
#plot(r3,pch=16)
r.fit4<-loess(r4~xyear.k,span=0.5)
plot(r.fit4,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit4),col=4)

MuMIn::AICc(Fit4)     # 131.4126
AIC_lag3<-MuMIn::AICc(Fit4)

plot(mod1, resid=T, pch=16, rug=F, se=F, ylab="Log(R/S)",xlab = "SC3 female abundance",main="Lag 3 S-R model - GAM")


r1<-resid(mod1)
#plot(1r,pch=16)
r.fit1<-loess(r1~dat$releaseyear,span=0.5)
plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(dat$releaseyear,fitted(r.fit1),col=4)

##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()
dev.new()
par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
r6<-resid(Fit4)
xyear<-xyear.k

### Raw abundances
#plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)
#axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
#axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
#box()
#title(xlab="EBS female abundance (Millions)")
#title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")


### Log-survival
plot(log(R/S)~S,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit4)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Female abundance (Millions)")
title(ylab="EBS log-survival")
title(main="SC3 females")


### Residual trend

#plot(r6,pch=16)
#r.fit3<-loess(r6~xyea.kr,span=0.30)
#plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
#lines(xyear.k,fitted(r.fit3),col=1)

r.fit3<-loess(r6~xyear.k,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear.k,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()
title(xlab="Hatch year")
title(ylab="EBS S-R residuals")
