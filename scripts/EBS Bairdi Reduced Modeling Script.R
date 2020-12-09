

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
library(MuMIN)
#?nlstools

getwd()
########## Import data and define variables####

Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)


#######################################################################################################
####################################### AICc formula ##################################################
# Copyright 2008 Alexander W Blocker
# You are free to use, modify, and redistribute this code under the terms of the GNU
# General Public License, version 2
#


AICc <- function(x, ... , k=2) {
  UseMethod("AICc");
}

AICc.default <- function(x, ... , k=2) {
  if ( length( list( ... ) ) ) {
    x <- list( x, ... );
    val <- lapply( x, logLik );
    val <- as.data.frame( t( sapply( val, function(el) c( attr( el, 
                                                                "df"), AICc( el, k = k) ) ) ) )
    names(val) <- c("df", "AICc")
    row.names(val) <- as.character(match.call()[-1])
    return(val)
  }
  AICc( logLik( x ) );
}

AICc.logLik <- function(x, ... , k=2) {
  -2 * c(x) + k * attr( x, "df" ) + 
    2 * attr( x, "df" ) * ( attr( x, "df" ) + 1 )/
    ( attr( x, "nobs" ) - attr( x, "df" ) - 1 );
}


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
names(Spawners)

sp_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_OVIGEROUS))))

sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))

colnames(sp_ovig)<-c("Year", "EBS_ovigfem_abun")
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")



Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3)))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3")
Spawners_SC3
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
legend(8,4e+08,c("EBS female estimate"),cex=1.25,col=c(4),pch=c(16))
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
legend(8,4e+08,c("Total EBS"),cex=1.25,col=c(1),pch=c(2))
#abline(h=0)
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
#title(main=" Reproductive female abundance by year")

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
y.range<-range(0,9e+08)

plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)

#abline(h=0)

axis(1, at=1:31,lab=Recruits$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("Total EBS"),cex=1.25,col=c(1),pch=c(2))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
#title(main=" Juvenile abundance by year")



#########################################################################################################################
###################################### Plot juv abundance separately ####################################################

Year<-Recruits$Year
plot(R~Year, type="b",ylab="Estimated abundance of EBS juvenile bairdi", xlab="Year", main = "Estimated abundance of EBS juvenile bairdi",pch=16)
abline(lm(R ~ Year), col=4)
fit <- loess(R ~ Year)
lines(x, fitted(fit), col=2)


#############################################################################################################################################
########################################### Recreate data series for actual analysis ########################################################
n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_SC3$EBS_SC3[4:n]

recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_SC3$Year[4:n]

R<-R[-5]
S<-S[-5]
recYear<-recYear[-5]
spYear<-spYear[-5]


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

R<-R[-5]
S<-S[-5]
recYear<-recYear[-5]
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
cor(R,S)
##############################################EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)


#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS group lag 3 Ricker model",pch=16)
abline(Fit3)
#plot(Fit3)

par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)


Fit4<-gls(log(R/S)~S,correlation=corARMA(p=2))
summary(Fit4)



Lm.Fit<-lm(log(R/S)~S, na.action=na.omit)
dwtest(Lm.Fit)
resid.series<-cbind(xyear.k,r3)
colnames(resid.series)<-c("HatchYear","Lag3SR_residual")
write.csv(resid.series,"data/EBS_Lag3_SR_residuals.csv")
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
par(mfrow=c(2,2))
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
anom1<-envar$EBS_NBT_RA3_final_year[11:45]
NBT.RA.3.end<-anom1

End.Year<-c(1986:2019)
RA.3.end.check<-cbind(End.Year,anom1)


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
############################################################ lag 3-predation release year e.g. 1983-2016 ###########################################################################

Pcod_lag3<-cod[9:42]

plot(Pcod_lag3~envar$Year[9:42])
#######################################################################################################################################################################
############################################################ lag 2-predation release year +1 e.g. 1984-2017 ###########################################################################

Pcod_lag2<-cod[10:43]
plot(Pcod_lag2~envar$Year[10:43])
#######################################################################################################################################################################
############################################################ lag 1-predation release year +2 e.g. 1985-2018 ###########################################################################

Pcod_lag1<-cod[11:44]

#######################################################################################################################################################################
############################################################ 3 year rolling average on end year, lagged for effect beginning hatch year e.g. 1983-2016 ###########################################################################
cod<-envar$Age3to7Pcodabun_RA3_end
Pcod_RA3_end_hatch<-cod[11:44]

plot(Pcod_RA3_end_hatch~envar$Year[11:44])
#######################################################################################################################################################################
############################################################ 3 year rolling average on end year, lagged for effect beginning year after hatch year e.g. 1984-2017 ###########################################################################
cod<-envar$Age3to7Pcodabun_RA3_end
Pcod_RA3_end_1ah<-cod[12:45]

#######################################################################################################################################################################
############################################################ 3 year rolling average on mid, lagged for effect beginning year after hatch year e.g. mid years = 1985-2018 ###########################################################################
cod<-envar$Age3to7Pcodabun_RA3_mid
Pcod_RA3_mid<-cod[11:44]
#########################################################################################################################################################################################
############################################################# Flathead sole model abundance vs. S/R residuals ##########################################################################
################################################################ Import data and define vatiables #######################################################################################
envar<-read.csv("data/envars.csv")
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
#################################################### Lag=3: Select data for effects during year of settling #########################################################
envar$Year[8:45]
FHS.anom
FHS.lag.3<-FHS.anom[2:35]
FHS.lag.3
settle.year<-c(1983:2016)
settle.year
FHS.a<-cbind(settle.year,FHS.lag.3)
FHS.a


FHS.lag.3

###############################################################################################################################################################################
#################################################### Lag=2: Select data for effects during first year after settling #########################################################

FHS.anom
FHS.lag.2<-FHS.anom[3:36]
FHS.lag.2

First.year<-c(1984:2017)
First.year
FHS.a<-cbind(First.year,FHS.lag.2)
FHS.a


FHS.lag.2
FHS.lag.2.ed

###############################################################################################################################################################################
#################################################### Lag=1: Select data for effects during second year after settling #########################################################

FHS.anom
FHS.lag.1<-FHS.anom[4:37]
FHS.lag.1
Second.year<-c(1985:2018)
Second.year
FHS.a<-cbind(Second.year,FHS.lag.1)
FHS.a


FHS.lag.1

###########################################################################################################################################
########################################## Set work directory and attach libraries ########################################################

############################################# PDO winter-release year ##########################################################################################
envar$Year
pdo<-envar$PDO_djf[9:42]
pdo
attach(pdo)

############################################# Summer PDO - release year ###############################################################################
pdos<-envar$PDO_jja[9:42]

pdos

pdos.ed<-pdos
############################################## Arctic oscillation - release year #########################################################################

AO<-envar$AO_jfm[9:42]


AO
AO.ed<-AO

############################################## Danielson wind data ##############################################################
SE.wind<-envar$SE.wind.May.Sep[9:42]
NW.wind<-envar$NW.wind.May.Sep[9:42]



####################################################################################################################################################################
############################################## FIT MODELS AND COMPARE ##############################################################################################

################################### baseline S/R model ############################################################################################################
fit0<-gls(log(R1/S1)~S1,correlation=corAR1(),method="ML")
summary(fit0)

r<-resid(fit0)

r.fit<-loess(r~xyear.k)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit),col=4)

######################################################################################################################################################################
#################################### Add Pcod ########################################################################################################################
######################################################################################################################################################################

################################### Add lag 2 Pcod linear term #######################################################################################################
cor.test(S1,Pcod_lag2) 

fit1<-gls(log(R1/S1)~S1+Pcod_lag2,correlation=corAR1(),method="ML")
summary(fit1)

r<-resid(fit1)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Add lag 2 Pcod linear and nonlinear term ########################################################################################

fit2<-gls(log(R1/S1)~S1+Pcod_lag2+I(Pcod_lag2^2),correlation=corAR1(),method="ML")
summary(fit2)

r<-resid(fit2)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Swap lag 1 for lag 2 Pcod ###################################################################################################################################
################################### Add lag 1 Pcod linear term - INSIGNIFICANT #######################################################################################################
cor.test(S1,Pcod_lag1) 

fit3<-gls(log(R1/S1)~S1+Pcod_lag1,correlation=corAR1(),method="ML")
summary(fit3)

r<-resid(fit3)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Add lag 1 Pcod linear and nonlinear term - BOTH SIGNIFICANT ########################################################################################
#################################### MODEL IS SIGNIFICANT ############################################################################################

fit4<-gls(log(R1/S1)~S1+Pcod_lag1+I(Pcod_lag1^2),correlation=corAR1(),method="ML")
summary(fit4)

r<-resid(fit4)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### lag 1 Pcod nonlinear term only - INSIGNIFICANT ########################################################################################

fit5<-gls(log(R1/S1)~S1+I(Pcod_lag1^2),correlation=corAR1(),method="ML")
summary(fit5)

r<-resid(fit5)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Pcod RA3_end beginning hatch year linear term only ########################################################################################

fit6<-gls(log(R1/S1)~S1+Pcod_RA3_end_hatch,correlation=corAR1(),method="ML")
summary(fit6)

r<-resid(fit6)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Pcod RA3_end beginning hatch year nonlinear term only ########################################################################################

fit6<-gls(log(R1/S1)~S1+I(Pcod_RA3_end_hatch^2),correlation=corAR1(),method="ML")
summary(fit6)

r<-resid(fit6)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Pcod RA3_end beginning hatch year linear and nonlinear terms ########################################################################################

fit7<-gls(log(R1/S1)~S1+I(Pcod_RA3_end_hatch^2)+Pcod_RA3_end_hatch,correlation=corAR1(),method="ML")
summary(fit7)

r<-resid(fit7)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)



################################### Pcod RA3_end beginning one year after hatch year linear term only ########################################################################################

fit8<-gls(log(R1/S1)~S1+Pcod_RA3_end_1ah,correlation=corAR1(),method="ML")
summary(fit8)

r<-resid(fit8)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Pcod RA3_end beginning one year after hatch year nonlinear term only ########################################################################################

fit9<-gls(log(R1/S1)~S1+I(Pcod_RA3_end_1ah^2),correlation=corAR1(),method="ML")
summary(fit9)

r<-resid(fit9)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Pcod RA3_end beginning one year after hatch year linear and nonlinear terms - BOTH TERMS ARE SIGNIFICANT ########################################################################################

#################################### MODEL IS SIGNIFICANT ############################################################################################

fit10<-gls(log(R1/S1)~S1+I(Pcod_RA3_end_1ah^2)+Pcod_RA3_end_1ah,correlation=corAR1(),method="ML")
summary(fit10)

r<-resid(fit10)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)






######################################################################################################################################################################
#################################### Add Flathead sole ###############################################################################################################
######################################################################################################################################################################

################################### Add lag 2 FHS linear term ########################################################################################################
################################### MODEL IS SIGNIFICANT #############################################################################################################
cor.test(S1,FHS.lag.2) 

fit11<-gls(log(R1/S1)~S1+FHS.lag.2,correlation=corAR1(),method="ML")
summary(fit11)

r<-resid(fit11)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Add lag 2 Pcod nonlinear term ###################################################################################################
################################### MODEL IS NOT SIGNIFICANT #################################################################################################################

fit11a<-gls(log(R1/S1)~S1+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit11a)

r<-resid(fit11a)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################### Add lag 2 FHS linear and nonlinear term ########################################################################################
################################### LINEAR TERM SIGNIFICANT, NONLINEAR IS NOT #######################################################################################

fit12<-gls(log(R1/S1)~S1+FHS.lag.2+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit12)

r<-resid(fit12)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

################################################################################################################################################################################
################################### Swap lag 1 for lag 2 FS ###################################################################################################################################
cor.test(S1,FHS.lag.1) 

fit13<-gls(log(R1/S1)~S1+FHS.lag.1,correlation=corAR1(),method="ML")
summary(fit13)

r<-resid(fit13)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")

################################### lag 1 FHS nonlinear term only - INSIGNIFICANT ########################################################################################

fit14<-gls(log(R1/S1)~S1+I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit14)

r<-resid(fit14)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

lines(xyear,fitted(r.fit),col=4)

################################### Add lag 1 FHS linear and nonlinear term - BOTH INSIGNIFICANT ########################################################################################

fit15<-gls(log(R1/S1)~S1+FHS.lag.1+I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit15)

r<-resid(fit15)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)


######################################################################################################################################################################
#################################### Add SST from ERSST###############################################################################################################
######################################################################################################################################################################

################################### Add linear term ########################################################################################################
################################### MODEL IS SIGNIFICANT #############################################################################################################
cor.test(S1,May.July.anom) 

fit16<-gls(log(R1/S1)~S1+May.July.anom,correlation=corAR1(),method="ML")
summary(fit16)

r<-resid(fit16)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

anova(fit1,fit2)
anova(fit3,fit4)


######################################################################################################################################################################
#################################### Add PDOw###############################################################################################################
######################################################################################################################################################################

################################### Add linear term ########################################################################################################
cor.test(S1,pdo) 

fit17<-gls(log(R1/S1)~S1+pdo,correlation=corAR1(),method="ML")
summary(fit17)

r<-resid(fit17)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

######################################################################################################################################################################
#################################### Add AO ################################################################################################################
######################################################################################################################################################################

################################### Add linear term ########################################################################################################
cor.test(S1,AO) 

fit18<-gls(log(R1/S1)~S1+AO,correlation=corAR1(),method="ML")
summary(fit18)

r<-resid(fit18)

r.fit<-loess(r~xyear)
plot(r.fit,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear,fitted(r.fit),col=4)

##################################################################################################################################
anova(fit1,fit2)
anova(fit3,fit4)

anova(fit1,fit2)
anova(fit3,fit4)
