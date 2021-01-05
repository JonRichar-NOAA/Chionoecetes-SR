

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
library(MuMIn)
library(tidyverse)
library(corrplot)
#?nlstools

########## Import data and define variables####

Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)


#########################################################################################################
########################### Create recruitment time series for analysis #################################
#########################################################################################################

rec_30to40<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO40 + Recruits$NUM_FEMALE_30TO40))))

rec_40to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_40TO50 + Recruits$NUM_FEMALE_40TO50))))

rec_50to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_50TO60 + Recruits$NUM_FEMALE_50TO60))))

rec_30to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO50 + Recruits$NUM_FEMALE_30TO50))))

rec_30to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO60 + Recruits$NUM_FEMALE_30TO60))))

rec_40to60<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_40TO50 + Recruits$NUM_FEMALE_40TO50+Recruits$NUM_MALE_50TO60 + Recruits$NUM_FEMALE_50TO60))))

colnames(rec_30to40)<-c("Year", "EBS_Abun_30to40")
colnames(rec_40to50)<-c("Year", "EBS_Abun_40to50")
colnames(rec_50to60)<-c("Year", "EBS_Abun_50to60")

colnames(rec_30to50)<-c("Year", "EBS_Abun_30to50")
colnames(rec_30to60)<-c("Year", "EBS_Abun_30to60")
colnames(rec_40to60)<-c("Year", "EBS_Abun_40to60")

#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################

sp_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_OVIGEROUS))))

sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))

colnames(sp_ovig)<-c("Year", "EBS_ovigfem_abun")
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")


Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3)))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3")

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
S<-Spawners_SC3$EBS_SC3
log.R<-log(R)
log.S<-log(S)
eda.norm(R)
eda.norm(log.R)


########################################################################################################################
########################################## ACF analysis of series by stanza ############################################
########################################################################################################################

########################################################################################################################
########################################## ACF for juvenile indices ####################################################
rec_30to50$Year[4:45] #check years
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

R<-rec_30to50$EBS_Abun_30to50[4:45]

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
Spawners_SC3$Year
S<-Spawners_SC3$EBS_SC3[4:45]
########################################## Full Timeseries #############################################################

acf(S,main="SC3 females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title
#?acf()

########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008

acf(S3,main="SC3 females 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
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
####################################### Plot all females  ######################################################

dev.new()
y.range<-range(0,4e+08)
par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
names(Spawners_SC3)
Spawners_SC3

EBS.fem<-Spawners_SC3$EBS_SC3

plot(EBS.fem)

plot(EBS.fem,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE,cex=1.25)


axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4)
#legend(8,4e+08,c("EBS female estimate"),cex=1.25,col=c(4),pch=c(16))
abline(h=0)
box()

title(xlab="Year")
#title(ylab="Reproductive female abundance")
title(main=" Reproductive female abundance by year")
########################################################################################################################
###################################### Repeat above customized for MEPS ###########################################
y.range<-range(0,4e+08)
par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
names(Spawners_SC3)

EBS.fem<-Spawners_SC3$EBS_SC3

plot(EBS.fem)

plot(EBS.fem,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4,labels=c("0","100","200","300","400"),tck=0.03)
#legend(8,4e+08,c("Total EBS"),cex=1.25,col=c(1),pch=c(2))
#abline(h=0)
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main=" Reproductive female abundance by year")

#########################################################################################################################
################################## Plot all juveniles together ##########################################################
names(rec_30to50)
n<-nrow(rec_30to50)

EBS.juv<-rec_30to50$EBS_Abun_30to50


par(mfrow=c(1,3))

par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.35,cex=1.35) #configure axis labels
y.range<-range(0,9e+08)

plot(EBS.juv,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)


axis(1, at=1:n,lab=rec_30to50$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("Total EBS"),cex=1,col=c(4),pch=c(22))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
#title(main=" Juvenile abundance by year")

#?axis()

########################################################################################################################
###################################### Repeat above customized for MEPS ###########################################
par(mfrow=c(1,1),cex.lab=1.4,cex.axis=1.4,cex=1.5) #configure axis labels
y.range<-range(0,4e+08)

plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)

#abline(h=0)

axis(1, at=1:n,lab=rec_30to50$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("Total EBS"),cex=1.25,col=c(1),pch=c(2))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main=" Juvenile abundance by year")



########################################################################################################
######################################### Lag=3 ########################################################
par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
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
plot(yi.k~xi.k,xlab="Reproductive female abundance", ylab="Lag 3 juvenile recruitment", main= "EBS lag 3 S/R comparison",pch=16,col=4,cex=1.25)

plot(yi.k~xi.k,xlab="X", ylab="Y", main= "EBS lag 3",pch=16,col=4,cex=1.25)
recYear_lag3<-yyear.k

##############################################EDA.norm analyses of log(R/S)#############################
#eda.norm(R)
log.RS<-log(R/S)
#eda.norm(log.RS)

cor(R,S)
#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
coef(Fit3)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS group lag 3 Ricker model",pch=16)
abline(Fit3)
#plot(Fit3)

par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)

##################################################################################################################################
################################################## Plot above for Publication #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()

par(mfrow=c(1,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
box()

title(xlab="EBS female abundance (Millions)")
title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")


plot(log(R/S)~S2,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Female abundance (Millions)")
title(ylab="EBS log-survival")


r.fit3<-loess(r6~xyear,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")

##############################################################################################################################################
################################################# Stock recruit models with environmental covariates integrated into model ####################
###############################################################################################################################################

######################################## Restrict spawner and recruit timeseries to same extent as envar data #################################

R1<-R[5:length(R)]
S1<-S[5:length(S)]
xyear<-xyear.k[5:length(xyear.k)]
#####################################################################################################################
#################################################     ERSST     #####################################################
#################################################import data set#####################################################
ersst0<-read.csv("data/ERSST_SST_avgs_anoms.csv")
ersst<-ersst0[6:nrow(ersst0),]
#########################################################################################################################################################
######################################## Edit Series for analyses using editted residual series and to match extent of FHS data set######################
May.anom<-ersst$May.anom
June.anom<-ersst$June.anom
July.anom<-ersst$July.anom
April.May.anom<-ersst$April.May.anom
May.June.anom<-ersst$May.June.anom
June.July.anom<-ersst$June.July.anom
April.July.anom<-ersst$April.July.anom
May.July.anom<-ersst$May.July.anom
ersst.year<-ersst$Year
ersst.year


May.July.anom
S1
##############################################################################################################################################
################################################# Add ENVAR database with non-custom ERSST and NCAR-NCEP and other additional datasets #######
################################################# Flathead sole data start with 1982, recruitmrnt series starts with 1983 due to previous and 1982 female outlier so need to standardize data inputs for start with 1982 cohort

envar<-read.csv("data/envars.csv")
envar
names(envar)

#########################################################################################################################
####################################### Danielson wind data from Mike ###################################################
names(envar)
SE.wind<-envar$SE.wind.May.Sep
year<-envar$Year

wind.dat<-cbind(year,SE.wind)
wind.dat

SE.wind.MS<-SE.wind[9:42]
##############################################################################################################################################
#################################################     NBT        #############################################################################
#################################################import data set##############################################################################

NBT<-envar$EBS_mean_NBT

n2<-length(envar$Year)
###############################################define variables#####################################################


ts_mean_NBT<-mean(NBT)
NBT.Anom<-NBT-ts_mean_NBT




###################################Lag=4 for effects during embryonic stage##########################################
xi<-NBT.Anom
xyear<-envar$Year
NBT.lag.4<-xi[8:(n2-4)]
xyear.lag.4<-xyear[8:(n2-4)]



#######################################  Lag=3 for effects in first summer######################################
######################################## Note: below only works for lag=3#######################################
xi<-NBT.Anom
xyear<-envar$Year

NBT.lag.3<-xi[9:n2]
xyear.lag.3<-xyear[9:n2]

#####################################################################################################################################################################################
######################################################### Lag=2 for effects over course of first winter##############################################################################
xi<-NBT.Anom
xyear<-envar$Year
NBT.lag.2<-xi[10:(n2)]
xyear.lag.2<-xyear[10:(n2)]

######################################################################################################################################################################################
##################################### 3 year rolling average of mean NBT for effects during juvenile stages, note: temp in year y= mean((y)+(y-1)+(y-2) ##############################
envar$EBS_NBT_RA3_final_year
xyear<-envar$Year
anom1<-envar$EBS_NBT_RA3_final_year[11:44]
NBT.RA.3.end<-anom1

End.Year<-c(1986:2019)
End.Year<-envar$Year[11:44]
RA.3.end.check<-cbind(End.Year,anom1)
RA.3.end.check

######################################################################################################################################################################################
##################################### 3 year rolling average of mean NBT for effects during juvenile stages, note: temp in year y= mean((y-1)+y+(y+1) ##############################
names(envar)
envar
envar$EBS_NBT_RA3_midyear
anom1<-envar$EBS_NBT_RA3_midyear[10:43]
anom1
n<-length(anom1)


Mid.Year<-c(1984:2017)
Mid.Year
RA.3.mid<-cbind(Mid.Year,anom1)
RA.3.mid


NBT.RA.3a<-anom1 #starts in 1980
NBT.RA.3a


#######################################################################################################################################################################################
##################################### 3 year rolling average of mean NBT for effects during embryonic and early juvenile stages, note: temp in year y= mean((y-1)+y+(y+1)###########

anom1<-envar$EBS_NBT_RA3_midyear[9:42]
anom1
n<-length(anom1)


Mid.Year<-c(1983:2016)
Mid.Year
RA.3.mid<-cbind(Mid.Year,anom1)
RA.3.mid


NBT.RA.3b<-anom1
NBT.RA.3b


########################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during juvenile stages, note: temp in year y= mean((y-1)+y ################################
names(envar)
envar$Year
envar$EBS_NBT_RA2
anom3<-envar$EBS_NBT_RA2[10:43]
anom3
First.Year<-c(1984:2017)
First.Year
NBT.RA.2a<-anom3



NBT.RA.2a

####################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during embryonic and 1st juvenile year, note: temp in year y= mean((y-1)+y)############

anom<-envar$EBS_NBT_RA2[9:42]
anom3
First.Year<-c(1983:2016)
First.Year
NBT.RA.2b<-cbind(First.Year,anom3)



NBT.RA.2b
###########################################################################################################################################################################################
############################################################# Age 3 to 7 Pacific cod abundance ############################################################################################
#######################################################################################################################################################################################

cod<-envar$Age3to7Pcodabun
year<-envar$Year

PC<-cbind(year,cod)
PC

#######################################################################################################################################################################
############################################################ lag 1-predation release year +2 e.g. 1985-2018 ###########################################################################

Pcod_lag1<-cod[11:44]

#########################################################################################################################################################################################
############################################################# Flathead sole model abundance vs. S/R residuals ##########################################################################
################################################################ Import data and define vatiables #######################################################################################
FHS<-envar$FHS_TBM[8:45]
avg.FHS<-mean(FHS)
FHS.anom<-FHS-avg.FHS
FHS.Year<-envar$Year[8:45]


plot(FHS.anom~FHS.Year,ylab="flathead sole age 3+ model biomass anomoly",xlab="model year",pch=16)
abline(h=0)

plot(FHS~FHS.Year,ylab="flathead sole age 3+ model biomass",xlab="model year",pch=16)
plot(FHS.rec~FHS.Year,ylab="flathead sole age 3 recruitment anomoly",xlab="model year",pch=16)
abline(h=0)


###############################################################################################################################################################################
#################################################### Lag=2: Select data for effects during first year after settling #########################################################

FHS.anom
FHS.lag.2<-FHS[3:36]
FHS.lag.2

First.year<-c(1984:2017)
First.year
FHS.a<-cbind(First.year,FHS.lag.2)
FHS.a


FHS.lag.2



############################################## Danielson wind data ##############################################################
SE.wind<-envar$SE.wind.May.Sep[9:42]



####################################################################################################################################################################
############################################## FIT MODELS AND COMPARE ##############################################################################################

############################## set up a data frame with all the data ##############################################
dat <- data.frame(releaseyear = c(1983:2016),
                  recr = R1 ,
                  spawn = S1,
                  logRS = log(R1/S1),
                  cod = Pcod_lag1,
                  fhs = FHS.lag.2,
                  temp = NBT.RA.3.end,
                  SE.wind = SE.wind) # or whatever parsimonious set of variables you want to use! 

dat
dat2<-dat[,2:ncol(dat)]
cor2 <- cor(dat2, use="complete.obs")
corrplot(cor2)

## examine distribution for recruits
## if skewed, may need a family other than Gaussian, e.g. Gamma

hist(dat$recr)   ## Data are skewed, use of family = Gamma is appropriate
hist(log(dat$recr))
hist(dat$recr^0.25)

########################################## now fit the full model ##################################################################################

mod1 <- gam(recr ~ s(spawn, k = 4) + s(cod, k = 4) + s(fhs, k = 4) + s(temp, k = 4) + s(SE.wind, k = 4),
            data = dat, family="Gamma",method = "REML")

summary(mod1)
plot(mod1, pages=1)

## drop lag 1 Pcod due to terrible p-value


mod2 <- gam(recr ~ s(spawn, k = 4) + s(fhs, k = 4) + s(temp, k = 4) + s(SE.wind, k = 4),
            data = dat,family="Gamma",method="REML")

summary(mod2)
plot(mod2, pages=1)
## drop Southeast wind due to terrible p-value


mod3 <- gam(recr ~ s(spawn, k = 4) + s(fhs, k = 4) + s(temp, k = 4),
            data = dat,family="Gamma",method="REML")

summary(mod3)
plot(mod3, pages=1) 
## note positive effect of FHS - not consistent with predation!
AICc(mod1, mod2, mod3)

########################################## use log ratio of R/S as response ######################################################################
hist(dat$logRS)
mod4 <- gam(logRS ~  s(spawn, k = 4) + s(cod, k = 4) + s(fhs, k = 4) + s(temp, k = 4) + s(SE.wind, k = 4),
            data = dat,method = "REML")

summary(mod4)
plot(mod4, pages=1)

## drop wind
mod5 <- gam(logRS ~ s(spawn, k = 4) + s(cod, k = 4) + s(fhs, k = 4) + s(temp, k = 4),
            data = dat,method="REML")

summary(mod5)
plot(mod5, pages=1)

## drop cod

mod6 <- gam(logRS ~ s(spawn, k = 4) + s(fhs, k = 4) + s(temp, k = 4),
            data = dat,method="REML")

summary(mod6)
plot(mod6, pages=1)

## and drop temp

mod7 <- gam(logRS ~ s(spawn, k = 4) + s(fhs, k = 4),
            data = dat,method="REML")

summary(mod7)
plot(mod7, pages=1, resid=T, pch=19)

AICc(mod4, mod5, mod6, mod7) # mod7 is marginaly better than mod6, much better than others

acf(resid(mod7)) # some autocorrelation at lag1

mod8 <- gam(logRS ~ s(log(spawn), k = 4) + s(fhs, k = 4),
            data = dat,method="REML")

summary(mod8)
plot(mod8, pages=1, resid=T, pch=19, shade=T, rug=F)

AICc(mod7, mod8)
acf(resid(mod8)) # same autocorrelation at lag1


# parameterize in GLS for good p-value in the presence of autocorrelation
fit1<-gls(logRS ~ log(spawn) + fhs, correlation = corAR1(), method="ML", data = dat)
summary(fit1)
# phi = 0.48

