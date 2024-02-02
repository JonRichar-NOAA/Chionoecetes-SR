

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

#R<-R[-8]
#S<-S[-5]
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
################################################# EXPORT DATA #############################################################################################
#dat <- data.frame(releaseyear = xyear.k,
#                  Lag3_recruits = R ,
#                  ReproductiveFemales = S,
#                  logRS = log(R/S)) # or whatever parsimonious set of variables you want to use! 

#write.csv(dat,"data/SR_data.csv")
##############################################################################################################################################
################################################# Stock recruit models with environmental covariates integrated into model ####################
###############################################################################################################################################

######################################## Restrict spawner and recruit timeseries to same extent as envar data #################################

#R1<-R[5:length(R)]
#S1<-S[5:length(S)]
#S1
#xyear<-xyear.k[5:length(S)]
####################################### For using as-is #############################################################
R1<-R
S1<-S
S1
xyear<-xyear.k
#####################################################################################################################
#################################################     ERSST     #####################################################
#################################################import data set#####################################################
ersst0<-read.csv("data/ERSST_SST_avgs_anoms.csv")

#for restricted length
#ersst<-ersst0[6:nrow(ersst0),]
#for full length
ersst<-ersst0[1:nrow(ersst0),]
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

#for restricted length
#SE.wind.MS<-SE.wind[9:42]

#for fulluseable length

SE.wind<-SE.wind[4:42]
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


NBT.lag.4<-xi[3:(n2-4)]
xyear.lag.4<-xyear[3:(n2-4)]

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
#envar$EBS_NBT_RA3_final_year
#xyear<-envar$Year
#restricted extent
#anom1<-envar$EBS_NBT_RA3_final_year[11:44]
#NBT.RA.3.end<-anom1

#full useable extent
anom1<-envar$EBS_NBT_RA3_final_year[6:44]
NBT.RA.3.end<-anom1

End.Year<-c(1980:2018)
End.Year<-envar$Year[6:44]
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
############################################################ lag 1-predation release year +2 e.g. 1980-2018 ###########################################################################

#Pcod_lag1<-cod[11:44]
#full extent (1980-2018)
Pcod_lag1<-cod[6:44]


#######################################################################################################################################################################
############################################################ lag 2-predation release year +1 e.g. 1979-2017 ###########################################################################

#Pcod_lag1<-cod[11:44]
#full extent (1980-2018)
Pcod_lag2<-cod[5:43]

########################################################################################################################################################################
############################################################ Data for 2 year Rolling average of Pcod abundance     ########################################################################
##################################################### Predation release year +1 and release year+2 (end years=1980 to 2018) #########################################################
Pcod_RA2<-envar$Age3to7Pcodabun_RA2_end[6:44]

########################################################################################################################################################################
############################################################ Data for 3 year Rolling average of Pcod abundance     ########################################################################
##################################################### Predation release year +1 and release year+2 and recruitment year (end years=1981 to 2019) #########################################################
Pcod_RA3<-envar$Age3to7Pcodabun_RA3_end[7:45]


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
#FHS.lag.2<-FHS[3:36]
#FHS.lag.2

#First.year<-c(1984:2017)
#First.year
#FHS.a<-cbind(First.year,FHS.lag.2)
#FHS.a


FHS.lag.2

# Full extent (1979-2018)
FHS.lag.2<-envar$FHS_TBM[5:43]
avg.FHS<-mean(FHS)
FHS.anom<-FHS-avg.FHS
FHS.Year<-envar$Year[5:43]

First.year<-c(1979:2017)
First.year
FHS.a<-cbind(FHS.Year,FHS.lag.2)
FHS.a


########################################################################################################################################################################
############################################################ Data for 2 year Rolling average of FHS TBM    ########################################################################
##################################################### Predation release year +1 and release year+2 (end years=1980 to 2018) #########################################################
FHS_RA2<-envar$FHS_TBM_RA2_end[6:44]
############################################## Danielson wind data ##############################################################
SE.wind<-envar$SE.wind.May.Sep[9:42]

#full extent

SE.wind<-envar$SE.wind.May.Sep[4:42]
############################################## Opilio data #############################################################################################
opies_edit<-read.csv("data/co_male_female_group_abun_eb_noCI_CVs.csv")
#opies_edit
#colnames(opies_edit)
#opies_edit$SURVEY_YEAR
#ovig_co_fem<-opies_edit$NUM_FEMALE_OVIGEROUS[4:37]

#full extent
opies_edit$SURVEY_YEAR
ovig_co_fem.a<-opies_edit$NUM_FEMALE_OVIGEROUS[1:37]
e<-c("NA","NA")
ovig_co_fem<-c(e,ovig_co_fem.a)
ovig_co_fem
####################################################################################################################################################################
############################## PDO and AO #########################################################################################################################

names(envar)
envar$Year

#pdo_djf<-envar$PDO_djf[9:42]
#AO_jfm<-envar$AO_jfm[9:42]

#Full extent
pdo_djf<-envar$PDO_djf[4:42]
AO_jfm<-envar$AO_jfm[4:42]

####################################################################################################################################################################
############################## 2 year rolling  for PDO and AO for releaae year and following year (end year = 1979-2017) ######################################################################################################
pdo_djf_RA2<-envar$PDO_djf_RA2_end[5:43]
AO_jfm_RA2<-envar$AO_jfm_RA2_end[5:43]

####################################################################################################################################################################
############################## 3 year rolling  for PDO and AO for releaae year and following 2 years (end year = 1980-2018) ######################################################################################################
pdo_djf_RA3<-envar$PDO_djf_RA3_end[5:43]
AO_jfm_RA3<-envar$AO_jfm_RA3_end[5:43]
############################################## EXPORT DATAFRAME ####################################################################################################
############################## set up a data frame with all the data ###############################################################################################
R
S
ovig_co_fem
Pcod_lag1
FHS.lag.2
NBT.RA.3.end
SE.wind
AO_jfm
pdo_djf
May.July.anom
length(R)
length(S)
length(ovig_co_fem)
length(Pcod_lag1)
length(FHS.lag.2)
length(NBT.RA.3.end)
length(SE.wind)
length(AO_jfm)
length(pdo_djf)
length(May.July.anom)
dat <- data.frame(releaseyear = c(1978:2016),
                  Lag3_recruits = R ,
                  ReproductiveFemales = S,
                  logRS = log(R/S),
                  Ovig_female_CO=ovig_co_fem,
                  Pcod_lag1 = Pcod_lag1,
                  PCod_RA2 = Pcod_RA2,
                  PCod_RA3 = Pcod_RA3,
                  FHS_lag2 = FHS.lag.2,
                  FHS_RA2 = FHS_RA2,
                  NBT_3RA = NBT.RA.3.end,
                   SE.wind = SE.wind,
                  AO_jfm = AO_jfm,
                  AO_RA2 = AO_jfm_RA2,
                  AO_RA3 = AO_jfm_RA3,
                  PDO_djf = pdo_djf,
                  PDO_RA2 = pdo_djf_RA2,
                  PDO_RA3 = pdo_djf_RA3,
                  SST_May_July = May.July.anom) # or whatever parsimonious set of variables you want to use! 

write.csv(dat,"data/EBS_Crab_and_envar_data_full_extent.csv")
####################################################################################################################################################################
############################################## FIT MODELS AND COMPARE ##############################################################################################


