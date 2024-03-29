
###### par(mfrow=c(1,1),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels

######set work directory and attach relevant libraries#####
#setwd('C:/Users/jon.richar/Work/GitRepos/Chionoecetes-SR/data')
#setwd('C:/Documents and Settings/Jon/Desktop/Project/Datasets for analysis')
 


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
###########################Define experimental inverse model#############################################
################################# Note: Model disabled ##################################################
#Inverse <- function(S,a){
#a+(1/S)
#}
#START<-list(a=4)

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
colnames(w166_sp_os)<-c("Year", "W166_sc4fem_abun")

Spawners_SC2<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2)))
Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
Spawners_SC4<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC4,E166Spawners$NUM_FEMALE_SC4,W166Spawners$NUM_FEMALE_SC4)))
Spawners_OS<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4),(E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4),(W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4))))

Spawners_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_OVIGEROUS,E166Spawners$NUM_FEMALE_OVIGEROUS,W166Spawners$NUM_FEMALE_OVIGEROUS)))


colnames(Spawners_SC2)<-c("Year","EBS_SC2","E166_SC2","W166_SC2")
Spawners_SC2

colnames(Spawners_SC3)<-c("Year","EBS_SC3","E166_SC3","W166_SC3")
Spawners_SC3

colnames(Spawners_ovig)<-c("Year","EBS_ovig","E166_ovig","W166_ovig")
Spawners_ovig

colnames(Spawners_OS)<-c("Year","EBS_os","E166_os","W166_os")
Spawners_ovig

Spawners_analysis<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3,
                                                 Spawners$NUM_FEMALE_OVIGEROUS,E166Spawners$NUM_FEMALE_OVIGEROUS,W166Spawners$NUM_FEMALE_OVIGEROUS,Spawners_OS$EBS_os,Spawners_OS$E166_os, Spawners_OS$W166_os)))

colnames(Spawners_analysis)<-c("Year","EBS_SC2","E166_SC2","W166_SC2","EBS_SC3","E166_SC3","W166_SC3","EBS_ovig","E166_ovig","W166_ovig","EBS_os","E166_os","W166_os")

Spawners_analysis
write.csv(Spawners_analysis,"./data/Female_Tanner_Crab_Series_for_analysis.csv")
#########################################################################################################
##################################### Analyses using UNEDITTED DATA######################################
######################################Plot of recruits vs spawners#######################################


#########################################################################################################
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)
#Graphically evaluate reproductive female data
for(i in names(Spawners_SC3)) {		# cycle through all names in data frame
	x <- Spawners_SC3$Year			# convert to numeric variable
	y <- Spawners_SC3[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
	abline(lm(y ~ x), col=4)
	# Fit a LOWESS (or LOESS) smooth trend and add to plot:
	fit <- loess(y ~ x)	
	lines(x, fitted(fit), col=2)	# Add fitted line
}
################################ Juveniles 30 to 40 #######################################################################
dev.new()
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

#Graphically evaluate juvenile data
for(i in names(rec_30to40)) {		# cycle through all names in data frame
	x <- rec_30to40$Year			# convert to numeric variable
	y <- rec_30to40[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
	abline(lm(y ~ x), col=4)
	# Fit a LOWESS (or LOESS) smooth trend and add to plot:
	fit <- loess(y ~ x)	
	lines(x, fitted(fit), col=2)	# Add fitted line
}

################################ Juveniles 40 to 50 #######################################################################
dev.new()
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

#Graphically evaluate juvenile data
for(i in names(rec_40to50)) {		# cycle through all names in data frame
	x <- rec_40to50$Year			# convert to numeric variable
	y <- rec_40to50[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
	abline(lm(y ~ x), col=4)
	# Fit a LOWESS (or LOESS) smooth trend and add to plot:
	fit <- loess(y ~ x)	
	lines(x, fitted(fit), col=2)	# Add fitted line
}

################################ Juveniles 50 to 60 #######################################################################
dev.new()
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

#Graphically evaluate juvenile data
for(i in names(rec_50to60)) {		# cycle through all names in data frame
	x <- rec_50to60$Year			# convert to numeric variable
	y <- rec_50to60[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
	abline(lm(y ~ x), col=4)
	# Fit a LOWESS (or LOESS) smooth trend and add to plot:
	fit <- loess(y ~ x)	
	lines(x, fitted(fit), col=2)	# Add fitted line
}

################################ Juveniles 30 to 50 #######################################################################
dev.new()
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

#Graphically evaluate juvenile data
for(i in names(rec_30to50)) {		# cycle through all names in data frame
	x <- rec_30to50$Year			# convert to numeric variable
	y <- rec_30to50[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
	abline(lm(y ~ x), col=4)
	# Fit a LOWESS (or LOESS) smooth trend and add to plot:
	fit <- loess(y ~ x)	
	lines(x, fitted(fit), col=2)	# Add fitted line
}

########################################################################################################################
####################################### EDA.norm analyses of R #########################################################
par(mfrow=c(1,1),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)
R<-rec_30to50$EBS_Abun_30to50[4:45]
S<-Spawners_ovig$EBS_ovig
log.R<-log(R)
log.S<-log(S)

eda.norm(R)
eda.norm(log.R)



########################################################################################################################
########################################## ACF analysis of series by stanza ############################################
########################################################################################################################

########################################################################################################################
########################################## ACF for juvenile indices ####################################################
rec_30to50$Year[4:45]
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

R<-rec_30to50$EBS_Abun_30to50[4:45]
S<-Spawners_ovig$EBS_ovig

########################################## Full Timeseries #############################################################

acf(R,main="Juvenile index, 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title
#?acf()

########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
R1<-R[1:21]  #1978-1998
R2<-R[22:42] #1999-2019
R3<-R[1:31]  #1978-2008

acf(R3,main="Juvenile index, 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R1,main="Juvenile index, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R2,main="Juvenile index, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

###############################################################################################################################
########################################## Now do ACF for female indices ######################################################
Spawners_ovig$Year
S<-Spawners_ovig$EBS_ovig[4:45]
########################################## Full Timeseries #############################################################

acf(S,main="SC2 females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title
#?acf()

########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008

acf(S3,main="SC2 females 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S1,main="SC2 females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="SC2 females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)


########################################################################################################################

for(i in names(Spawners_ovig)) {		# cycle through all names in data frame
	x <- Spawners_ovig$Year			# convert to numeric variable
	y <- Spawners_ovig[,i]
	# Plot time series for variable i using both points and lines. 
	# Suppress default axis labels and add variable name as title:
	plot(x, y, type = "b", xlab="", ylab="", main = i,pch=16)
	# Fit linear time trend and add line to plot:
}

########################################################################################################################
####################################### Plot all females  ######################################################


y.range<-range(0,4e+08)
par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
names(Spawners_ovig)
Spawners_SC2
W166.fem<-Spawners_ovig$W166_ovig
EBS.fem<-Spawners_ovig$EBS_ovig
E166.fem<-Spawners_ovig$E166_ovig
plot(EBS.fem)
W166.fem
plot(EBS.fem,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)
lines(W166.fem,type="b",pch=16,col=1)
lines(E166.fem,type="b",pch=22,col=2)

axis(1, at=1:45,lab=Spawners_SC2$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4)
legend(8,4e+08,c("EBS female estimate", "Western group female estimate", "Eastern group female estimate"),cex=1,col=c(4,1,2),pch=c(16,20,22))
abline(h=0)
box()

title(xlab="Year")
#title(ylab="Reproductive female abundance")
title(main=" Reproductive female abundance by year")





#########################################################################################################################
########################################### Lag=2 #######################################################################
par(mfrow=c(2,2))

#########################################Define data for lag=2####################http://127.0.0.1:44653/graphics/5095cdb0-250e-4a89-b3e1-bea9e1861ea9.png
par(mfrow=c(2,2))
R<-rec_30to50$EBS_Abun_30to50[4:45]
S<-Spawners_ovig$EBS_ovig
log.R<-log(R)
log.S<-log(S)


xi<-S
xyear<-Spawners$Year
yi<-R
yyear<-Recruits$Year
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
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 2 juvenile recruitment", pch=16,col=4)

###################################EDA.norm analysis of log(R/S)########################
log.RS<-log(R/S)
eda.norm(log.RS)

########################################Fit Ricker model#########################################

Rick.fit<-nls(log(R/S) ~ Ricker(S,a,b),start=START, na.action=na.omit)
summary(Rick.fit)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 2 log(Recruits/Repfem)",pch=16)
abline(Rick.fit)

#plot(Rick.fit)
#identify(S,log(R/S))

r0<-resid(Rick.fit)
#plot(r0,pch=16)
r.fit2<-loess(r0~xyear.k)
plot(r.fit2,pch=16,xlab="Hatch year", ylab="Lag 2 S/R residuals")
#lines(r.fit2)
lines(xyear.k,fitted(r.fit2))


####################################Fit gls models with 1st and 2nd order autocorrelation######################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
 summary(Fit3)
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 2 log(Recruits/Repfem)",pch=16)
abline(Fit3)
plot(Fit3)

r1<-resid(Fit3)
#plot(1r,pch=16)
r.fit3<-loess(r1~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 2 S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)


Fit4<-gls(log(R/S)~S,correlation=corARMA(p=2))
summary(Fit4)

########################################################################################################
######################################### Lag=3 ########################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
R<-rec_30to50$EBS_Abun_30to50[4:45]
S<-Spawners_ovig$EBS_ovig
#release.year<-c(1978:2005)

xi<-S
xyear<-Spawners_ovig$Year
xyear
yi<-R
yyear<-rec_30to50$Year
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
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 3 juvenile recruitment", main= "EBS lag 3 S/R comparison",pch=16,col=4,cex=1.25)

cor(R,S)
##############################################EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)

##############################################Fit Ricker model##########################################

Rick.fit2<-nls(log(R/S) ~ Ricker(S,a,b),start=START)
summary(Rick.fit2)
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS lag 3 Ricker model",pch=16)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS lag 3 S/R",pch=16)

abline(Rick.fit2)


r2<-resid(Rick.fit2)
r.fit2<-loess(r2~xyear.k)
plot(r.fit2,xlab="Hatch year",ylab="EBS lag 3 S/R residuals",pch=16)
lines(xyear.k,fitted(r.fit2))

par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
plot(r.fit2,xlab="Hatch year",ylab="EBS lag 3 S/R residuals",pch=16)
lines(xyear.k,fitted(r.fit2))
par(mfrow=c(2,2))

#lines(r.fit2)
#plot(r2,pch=16)
#identify(S,log(R/S))
#plot(Rick.fit2)

#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS group lag 3 Ricker model with ovigerous females",pch=16)
abline(Fit3)
#plot(Fit3)

par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k,span=0.25)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)


Fit4<-gls(log(R/S)~S,correlation=corARMA(p=2))
summary(Fit4)



Lm.Fit<-lm(log(R/S)~S, na.action=na.omit)
dwtest(Lm.Fit)


#########################################################################################################################
############################################### Drop data point #5  #####################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25)
R2<-R[-5]
R2
R
S2<-S[-5]

plot(R2~S2,xlab= "Reproductive female abundance", ylab="Juvenile abundance",main= "EBS juveniles vs reproductive female abundance",pch=16)
Fit3<-gls(log(R2/S2)~S2,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
plot(log(R2/S2)~S2, xlab="Reproductive female abundance",ylab="LN(Recruits/Repfem)",main="EBS log-survival vs. abundance of reproductive females",pch=16)
abline(Fit3)
#plot(Fit3)


par(mfrow=c(1,1))
r6<-resid(Fit3)
xyear<-xyear.k[-5]
#plot(r6,pch=16)
r.fit3<-loess(r6~xyear,span=0.45)
plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S/R residuals", main=" EBS lag 3 S/R residuals vs. hatch year")
lines(xyear,fitted(r.fit3),col=4)

cor(R2,S2)

##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()

par(mfrow=c(1,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
plot(R2~S2,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
box()

title(xlab="EBS female abundance (Millions)")
title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")


plot(log(R2/S2)~S2,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Female abundance (Millions)")
title(ylab="EBS log-survival")

r6<-resid(Fit3)
xyear<-xyear.k[-5]
#plot(r6,pch=16)
r.fit3<-loess(r6~xyear,span=0.45)
plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
lines(xyear,fitted(r.fit3),col=1)

r.fit3<-loess(r6~xyear,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005),labels=c("1980","1985","1990","1995","2000","2005"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")

##############################################################################################################################
################################################# Lag=4 ######################################################################
par(mfrow=c(1,1))
R<-Juv_abun
S<-Repfem_abun
release.year<-c(1978:2004)

xi<-S
xyear<-Spawners$Year
yi<-R
yyear<-Recruits$Year
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
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 4 juvenile recruitment",pch=16,col=4)


############################################## EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)

##################################################################################################################################
################################################### Fit Ricker model##############################################################

Rick.fit3<-gnls(log(R/S) ~ Ricker(S,a,b),start=START)
summary(Rick.fit3)
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 4 log(Recruits/Repfem)",pch=16)
abline(Rick.fit2)


r4<-resid(Rick.fit3)
r.fit4<-loess(r4~xyear.k)
plot(r.fit3,xlab="Hatch year",ylab="Lag 4 S/R residuals",pch=16)
lines(xyear.k,fitted(r.fit4))

#lines(r.fit2)
#plot(r2,pch=16)
#identify(S,log(R/S))
#plot(Rick.fit2)

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

Fit5<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit5)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 4 log(Recruits/Repfem)",pch=16)
abline(Fit5)
#plot(Fit3)

par(mfrow=c(1,1))
r5<-resid(Fit3)
#plot(r5,pch=16)
r.fit5<-loess(r5~xyear.k)
plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 4 S/R residuals")
lines(xyear.k,fitted(r.fit5),col=4)


Fit6<-gls(log(R/S)~S,correlation=corARMA(p=2))
summary(Fit6)





