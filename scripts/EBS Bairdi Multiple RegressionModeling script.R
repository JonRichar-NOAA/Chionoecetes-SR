

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
#?nlstools

########## Import data and define variables####

Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)

E166Recruits<-read.csv("data/cb_e166_pop_juvenile.csv")
E166Spawners<-read.csv("data/cb_e166_pop_sc_cs.csv")
W166Recruits<-read.csv("data/cb_w166_pop_juvenile.csv")
W166Spawners<-read.csv("data/cb_w166_pop_sc_cs.csv")

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

#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################
names(Spawners)

sp_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_OVIGEROUS))))
sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))
colnames(sp_ovig)<-c("Year", "EBS_ovigfem_abun")
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")

e166_sp_ovig<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_OVIGEROUS))))
e166_sp_sc3<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3))))
colnames(e166_sp_ovig)<-c("Year", "E166_ovigfem_abun")
colnames(e166_sp_sc3)<-c("Year", "E166_sc3fem_abun")

w166_sp_ovig<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_OVIGEROUS))))
w166_sp_sc3<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3))))
colnames(w166_sp_ovig)<-c("Year", "W166_ovigfem_abun")
colnames(w166_sp_sc3)<-c("Year", "W166_sc3fem_abun")

Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3","E166_SC3","W166_SC3")
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
W166.fem<-Spawners_SC3$W166_SC3
EBS.fem<-Spawners_SC3$EBS_SC3
E166.fem<-Spawners_SC3$E166_SC3
plot(EBS.fem)
W166.fem
plot(EBS.fem,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE,cex=1.25)
lines(W166.fem,type="b",pch=16,col=1,cex=1.25)
lines(E166.fem,type="b",pch=22,col=2,cex=1.25)

axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4)
legend(8,4e+08,c("EBS female estimate", "Western group female estimate", "Eastern group female estimate"),cex=1.25,col=c(4,1,2),pch=c(16,20,22))
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
W166.fem<-Spawners_SC3$W166_SC3
EBS.fem<-Spawners_SC3$EBS_SC3
E166.fem<-Spawners_SC3$E166_SC3
plot(EBS.fem)

plot(EBS.fem,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)
lines(W166.fem,type="b",pch=20,col=1)
lines(E166.fem,type="b",pch=22,col=1)

axis(1, at=1:45,lab=Spawners_SC3$Year,tck=0.03)
axis(2, las=1, at=1e+08*0:4,labels=c("0","100","200","300","400"),tck=0.03)
legend(8,4e+08,c("Total EBS", "Western area", "Eastern area"),cex=1.25,col=c(1,1,1),pch=c(2,20,22))
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
W166.juv<-rec_30to50$W166_Abun_30to50
E166.juv<-rec_30to50$E166_Abun_30to50

par(mfrow=c(1,3))

plot(EBS.juv)
plot(E166.juv)
plot(W166.juv)

par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.35,cex=1.35) #configure axis labels
y.range<-range(0,9e+08)

plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)
lines(W166.juv,type="b",pch=18,col=1)
lines(E166.juv,type="b",pch=22,col=2)
#abline(h=0)

axis(1, at=1:n,lab=rec_30to50$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,9e+08,c("Total EBS", "Western area", "Eastern area"),cex=1,col=c(4,1,2),pch=c(16,18,22))
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
#########################################################################################################################
########################################### Lag=2 #######################################################################
par(mfrow=c(2,2))

plot(R~recYear, ylab="Estimated juvenile abundance", xlab="Year", pch=16)
plot(S~spYear, ylab="Estimated reproductive female crab abundance", xlab="Year",pch=16)




#########################################Define data for lag=2####################
par(mfrow=c(2,2))

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
n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_SC3$EBS_SC3[4:n]

R
S
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_SC3$Year[4:n]
release.year<-c(1978:2005)

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

cor(R,S)
##############################################EDA.norm analyses of log(R/S)#############################
eda.norm(R)
log.RS<-log(R/S)
eda.norm(log.RS)

##############################################Fit Ricker model##########################################

Rick.fit2<-nls(log(R/S) ~ Ricker(S,a,b),start=START)
summary(Rick.fit2)
plot(log(R/S)~S, xlab="Reproductive SC3 female abundance",ylab="lag 3 log(Recruits/Repfem)",main="EBS lag 3 Ricker model",pch=16)
abline(Rick.fit2)

dev.new()
plot(log(R/S)~S, xlab="SC3 ovigerous female abundance",ylab="log(Juveniles/females)",main="EBS lag 3",pch=16)


r2<-resid(Rick.fit2)
r.fit2<-loess(r2~xyear.k,span=0.45)
plot(r.fit2,xlab="Hatch year",ylab="EBS lag 3 S/R residuals",pch=16)
lines(xyear.k,fitted(r.fit2))

plot(r.fit2,xlab="Year",ylab="Y",pch=16)
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

plot(r.fit3,pch=16,xlab="Year", ylab=" Y")
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
#######################################################################################################################################
############################################### Designate residual series to use in modeling##########################################
######################################################################################################################################
r
r<-r3

#########################################################################################################################################
########################################## Analyze S/R residuals against environmental variables ######################################################
#########################################################################################################################################

Year<-c(1981:2019)
resid.series<-cbind(Year,r)  #pairs year of recruitment to index with S/R residual of that pseudocohort
resid.series
plot(resid.series,ylab="Lag 3 S/R residual",xlab="Year",pch=16)



###############################################################################################################################
################################################# NCAR NCEP wind data ########################################################
###############################################################################################################################
wind_dat<-read.csv("data/wind_vector_data.csv")

wind_dat
#####################################################################################################################
#################################################     ERSST     #####################################################
#################################################import data set#####################################################
ersst<-read.csv("data/ERSST_SST_avgs_anoms.csv")
#####################################################################################################################
#################################################     Survey SST     ################################################
#################################################import data set#####################################################
par(mfrow=c(2,2))
SST<-read.csv("EBS_SST.csv")
attach(SST)

###############################################define variables#####################################################

names(SST)
#Year<-SST$Year
SST.avg<-Avg.SST
SST.Anom<-SST.Anom
SST$Year
plot(SST.Anom~SST$Year,pch=16,ylab="SST anomoly",xlab="Year")
abline(h=0)



#######################################  Lag=3 for effects in first summer######################################
######################################## Note: below only works for lag=3#######################################
xi<-SST.Anom
xyear<-SST$Year
yi<-r
yyear<-Year

n <- length(yi)
xi.t<-xi[1:n]
xyear.t<-xyear[1:n]
yi
Year
xi.t
xyear.t


##############################################################################################################################################
#################################################     NBT        #############################################################################
#################################################import data set##############################################################################
par(mfrow=c(2,2))
NBT<-read.csv("NBT_Anomoly.csv")
attach(NBT)

yi<-r
n <- length(yi)
##############################################3#define variables#####################################################

names(NBT)
#Year<-NBT$Year
NBT.mean<-Mean_NBT
NBT.Anom<-NBT_Anom
NBT$Year
plot(NBT.Anom~NBT$Year,pch=16,ylab="NBT anomoly",xlab="Year")
abline(h=0)



###################################Lag=4 for effects during embryonic stage##########################################

xi<-NBT.Anom
xyear<-NBT$Year

lag.4.juvs<-yi[2:n]
yyear.t4<-Year[2:n]


NBT.lag.4<-xi[1:(n-1)]
xyear.lag.4<-xyear[1:(n-1)]

NBT.lag.4
xyear.lag.4
lag.4.juvs
yyear.t4

#######################################  Lag=3 for effects in first summer######################################
######################################## Note: below only works for lag=3#######################################
xi<-NBT.Anom
xyear<-NBT$Year

NBT.lag.3<-xi[1:n]
xyear.lag.3<-xyear[1:n]

NBT.lag.3
xyear.t3

NBT.lag.3

#####################################################################################################################################################################################
######################################################### Lag=2 for effects over course of first winter##############################################################################
names(NBT)
xi<-NBT.Anom
xyear<-NBT$Year
xyear
NBT.lag.2<-xi[2:(n+1)]
NBT.lag.2
xyear.t2<-xyear[2:(n+1)]
xyear.t2


NBT.lag.2

######################################################################################################################################################################################
##################################### 3 year rolling average of NBT anomoly for ffects during juvenile stages, note: temp in year y= mean((y-1)+y+(y+1) ##############################

anom1<-NBT$RA_3[1:29]
anom1
Mid.Year<-c(1979:2007)
Mid.Year
RA.3<-cbind(Mid.Year,anom1)
RA.3
NBT.RA.3a<-anom1[1:28]
NBT.RA.3a
xyear<-Mid.Year[1:28]
xyear


NBT.RA.3a

#######################################################################################################################################################################################
##################################### 3 year rolling average of NBT anomoly for effects during embryonic and early juvenile stages, note: temp in year y= mean((y-1)+y+(y+1)###########

anom2<-NBT$RA_3[1:29]
anom2
Mid.Year<-c(1979:2007)
Mid.Year
RA.3<-cbind(Mid.Year,anom2)
RA.3
NBT.RA.3b<-anom2[1:27]
xyear<-Mid.Year[1:27]



NBT.RA.3b

########################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during juvenile stages, note: temp in year y= mean((y-1)+y+(y+1) ################################

NBT$RA_2
anom3<-NBT$RA_2[1:30]
anom3
First.Year<-c(1978:2007)
First.Year
RA.2<-cbind(First.Year,anom3)
RA.2
NBT.RA.2a<-anom3[1:28]
xyear<-First.Year[1:28]
xyear


RA.2a.temp<-cbind(xyear,NBT.RA.2a)
RA.2a.temp


NBT.RA.2a

####################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during embryonic and 1st juvenile year, note: temp in year y= mean((y-1)+y+(y+1)############

NBT$RA_2
anom<-NBT$RA_2[1:30]
anom
Embryo.Year<-c(1978:2007)
Embryo.Year

NBT.RA.2b<-anom[1:27]
xyear<-Embryo.Year[1:27]
xyear
yi<-r[2:28]
yyear<-Year[2:28]
yi
yyear

RA.2b.temp<-cbind(xyear,NBT.RA.2b)
RA.2b.temp
RA.2b.recruit<-cbind(yyear,yi)
RA.2b.recruit


NBT.RA.2b


##############################################################################################################################################
#################################################   NBT using editted residual series  #############################################################################
#################################################import data set##############################################################################

par(mfrow=c(2,2))
NBT<-read.csv("NBT_Anomaly.csv")
attach(NBT)

##############################################3#define variables#####################################################

names(NBT)
#Year<-NBT$Year
NBT.mean<-Mean_NBT
NBT.Anom<-NBT_Anom
NBT$Year
plot(NBT.Anom~NBT$Year,pch=16,ylab="NBT anomoly",xlab="Year")
abline(h=0)



###################################Lag=4 for effects during embryonic stage##########################################

xi<-NBT.Anom
xyear<-NBT$Year

NBT.lag.4<-xi[1:27]
NBT.lag.4.ed<-NBT.lag.4[-4]

xyear
xyear.t4a<-xyear[1:27]
xyear.t4<-xyear.t4a[-4]


NBT4.series<-cbind(xyear.t4,NBT.4)
NBT4.series


NBT.lag.4.ed

##############################################################################################################
#######################################  Lag=3 for effects in first summer######################################
######################################## Note: below only works for lag=3#######################################
xi<-NBT.Anom
xyear<-NBT$Year



#n <- length(yi)
NBT.lag.3<-xi[1:28]
NBT.lag.3.ed<-NBT.lag.3[-5]

xyear
xyear.t3<-xyear[1:28]
NBTyear.lag.3.ed<-xyear.t3a[-5]


NBT.3
NBT.3a

NBT3.series<-cbind(xyear.t3,NBT.lag.3)
NBT3a.series<-cbind(xyear.t3a,NBT.lag.3a.edit)

NBT3.series
NBT3a.series



NBT.lag.3.ed

###################################################################################################################################################################################
######################################################### Lag=2 for effects over course of first winter##############################################################################
names(NBT)
xi<-NBT.Anom
xyear<-NBT$Year
xi
xyear

NBT.series<-cbind(xyear,xi)
NBT.series


NBT.lag.2<-xi[2:29]
NBT.lag.2.ed<-NBT.lag.2[-6]
xyear.t2<-xyear[2:29]
xyear.t2a<-xyear.t2a[-6]


NBT2.series<-cbind(xyear.t2,NBT.lag.2.edit)
NBT2.series


NBT.lag.2.ed

######################################################################################################################################################################################
##################################### 3 year rolling average of NBT anomoly for effects during juvenile stages, note: temp in year y= mean((y-1)+y+(y+1) ##############################

anom<-NBT$RA_3[1:29]
anom
Mid.Year<-c(1979:2007)
Mid.Year
RA.3<-cbind(Mid.Year,anom)
RA.3

NBT.RA.3a<-anom[1:28]
NBT.RA.3a.ed<-NBT.RA.3a[-5]
xyear.3<-Mid.Year[1:28]
xyear<-xyear.3[-5]
xyear


RA.3.series<-cbind(xyear.3,RA.3)
RA.3a.series<-cbind(xyear.3a,RA.3a)



NBT.RA.3a.ed

#######################################################################################################################################################################################
####################################1 3 year rolling average of NBT anomoly for effects during embryonic and early juvenile stages, note: temp in year y= mean((y-1)+y+(y+1)###########

anom<-NBT$RA_3[1:29]
anom
Mid.Year<-c(1979:2007)
Mid.Year
RA.3<-cbind(Mid.Year,anom)
RA.3
NBT.RA.3b<-anom[1:27]
NBT.RA.3b.ed<-NBT.RA.3b[-4]
xyear.3<-Mid.Year[1:27]
xyear<-xyear.3[-4]
xi
xyear



NBT.RA.3b.ed

########################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during juvenile stages, note: temp in year y= mean((y-1)+y+(y+1) ################################

NBT$RA_2
anom2a<-NBT$RA_2[1:30]
anom2a
First.Year<-c(1978:2007)
First.Year
RA.2<-cbind(First.Year,anom2a)
RA.2
NBT.RA.2a<-anom2a[1:28]
NBT.RA.2a.ed<-NBT.RA.2a[-5]
xyear.2<-First.Year[1:28]
xyear<-xyear.2[-5]
xyear




RA.2a.temp<-cbind(xyear,xi)
RA.2a.temp



NBT.RA.2a.ed

####################################################################################################################################################################################
##################################### 2 year rolling average of NBT anomoly for effects during embryonic and 1st juvenile year, note: temp in year y= mean(y+y+1)############

NBT$RA_2
anom2b<-NBT$RA_2[1:30]
anom2b
Embryo.Year<-c(1978:2007)
Embryo.Year

NBT.RA.2b<-anom2b[1:27]
NBT.RA.2b

NBT.RA.2b.ed<-NBT.RA.2a[-4]
NBT.RA.2b.ed

xyear.2<-Embryo.Year[1:27]
xyear.2
xyear<-xyear.2[-4]
xyear



NBT.RA.2b.ed



##############################################################################################################################################
#################################################     Large males         #############################################################################
#################################################import data set##############################################################################

par(mfrow=c(1,1),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
Males<-read.csv("Large_males_TS.csv")
attach(Males)
names(Males)

#Year<-Males$Year
males<-Male_abun
male.anom<-Male_anom
Males$Year
plot(males~Males$Year,pch=16,ylab="Large male abundance",xlab="Year")
abline(h=0)

yi<-r
n <- length(yi)

#######################################  Lag=3 for predation on settling #######################################
######################################## Note: below only works for lag=3#######################################
xi<-males
xyear<-Males$Year

lag.3.males<-xi[1:n]
xyear.t3<-xyear[1:n]

lag.3.males
xyear.t3


#####################################################################################################################################################################################
######################################################### Lag=2 for predation during first juve year ##############################################################################
names(Males)
xi<-males
xyear<-Males$Year
xyear

lag.2.males<-xi[2:29]
lag.2.males
xyear.t2<-xyear[2:29]
xyear.t2


lag.2.males

#####################################################################################################################################################################################
######################################################### Lag=1 for predation during second juv year ##############################################################################
names(Males)
xi<-males
xyear<-Males$Year
xyear

lag.1.males<-xi[3:30]
lag.1.males
xyear.t1<-xyear[3:30]
xyear.t1

lag1.males<-cbind(xyear.t1,lag.1.males)
lag1.males

######################################################################################################################################################################################
##################################### 3 year rolling average of male abundance for predation through juvenile stages, note: abundance in year y= mean((y-1)+y+(y+1) ##################
Males$RA.3

RA.3.a<-Males$RA.3[1:29]
RA.3.a
Mid.Year<-c(1979:2007)
Mid.Year
RA.3<-cbind(Mid.Year,RA.3.a)
RA.3

RA.3.males<-RA.3.a[2:29]
xyear<-Mid.Year[2:29]

RA.3.b<-cbind(xyear,RA.3.males)
RA.3.b

RA.3.males

########################################################################################################################################################################################
##################################### 2 year rolling average of male abundance for predation throughout juvenile stages, note: abun in year y= mean((y-1)+y+(y+1) ################################
Males$RA.2
anom1<-Males$RA.2[1:30]
anom1

First.Year<-c(1978:2007)
First.Year

xyear<-First.Year[1:28]
xyear



RA.2a.males<-anom1[1:28]

RA.2<-cbind(First.Year,anom1)
RA.2
RA.2a.males.check<-cbind(xyear,RA.2a.males)
RA.2a.males.check


RA.2a.males

#######################################################################################################################################################################################
##################################### 2 year rolling average of male abundance for effects during 2nd juvenile year and beginning of third, note: temp in year y= mean((y-1)+y+(y+1)###

Males$RA.2
anom2<-Males$RA.2[1:30]
anom2
Second.Year<-c(1979:2007)
Second.Year

RA.2bb<-cbind(Second.Year,anom2)
RA.2bb

xi.b<-anom2[2:30]

xyear<-Second.Year
xyear

RA.2b.check<-cbind(xyear,xi.b)
RA.2b.check

RA.2b.males<-xi.b[2:29]
RA.2b.males


RA.2b.recruit<-cbind(yyear,yi)
RA.2b.recruit

###############################################################################################################################################
################################################## Generate editted series for use with editted residual set ##################################

lag.3.males.ed<-lag.3.males[-5]	
lag.2.males.ed<-lag.2.males[-5]	
lag.1.males.ed<-lag.1.males[-5]
RA.3.males.ed<-RA.3.males[-5]
RA.2a.males.ed<-RA.2a.males[-5]	
RA.2b.males.ed<-RA.2b.males[-5] 

#########################################################################################################################################################################################
############################################################# Flathead sole model abundance vs. S/R residuals ##########################################################################
################################################################ Import data and define vatiables #######################################################################################

FHS.TS<-read.csv("FHS_SAFE.csv")
names(FHS.TS)
FHS.TS
attach(FHS.TS)
FHS<-Tot_biomass
FHS.anom<-BM_anom
FHS.Year<-FHS.TS$Year
FHS.rec<-Rec_anom

plot(FHS.anom~FHS.Year,ylab="flathead sole age 3+ model biomass anomoly",xlab="model year",pch=16)
abline(h=0)

plot(FHS~FHS.Year,ylab="flathead sole age 3+ model biomass",xlab="model year",pch=16)
plot(FHS.rec~FHS.Year,ylab="flathead sole age 3 recruitment anomoly",xlab="model year",pch=16)
abline(h=0)



###############################################################################################################################################################################
#################################################### Lag=3: Select data for effects during year of settling #########################################################

FHS.anom
FHS.lag.3<-FHS.anom[2:29]
FHS.lag.3
FHS.lag.3.ed<-FHS.lag.3[-5]
settle.year<-c(1978:2005)
settle.year
FHS.a<-cbind(settle.year,FHS.lag.3)
FHS.a


FHS.lag.3
FHS.lag.3.ed
###############################################################################################################################################################################
#################################################### Lag=2: Select data for effects during first year after settling #########################################################
r
FHS.anom
FHS.lag.2<-FHS.anom[3:30]
FHS.lag.2
FHS.lag.2.ed<-FHS.lag.2[-5]
First.year<-c(1979:2006)
First.year
FHS.a<-cbind(First.year,FHS.lag.2)
FHS.a


FHS.lag.2
FHS.lag.2.ed

###############################################################################################################################################################################
#################################################### Lag=1: Select data for effects during second year after settling #########################################################

FHS.anom
FHS.lag.1<-FHS.anom[4:31]
FHS.lag.1
FHS.lag.1.ed<-FHS.lag.1[-5]
Second.year<-c(1980:2007)
Second.year
FHS.a<-cbind(Second.year,FHS.lag.1)
FHS.a


FHS.lag.1
FHS.lag.1.ed
###########################################################################################################################################
########################################## Set work directory and attach libraries ########################################################

setwd('C:/Users/Jon Richar/Desktop/Project/Datasets for Analysis/BSC.data')#Valkrie		
setwd('C:/Documents and Settings/Jon/Desktop/Project/Datasets for analysis/BSC.data')#Paladin

############################################# PDO ##########################################################################################
pdo<-read.csv("PDO.csv")
pdo
attach(pdo)

pdos<-PDOs[79:106]

pdos

pdos.ed<-pdos[-5]
pdos.ed
############################################## Arctic oscillation #########################################################################


Arctic.osc<-read.csv("Arctic_oscillation.csv")
Arctic.osc
attach(Arctic.osc)
AO<-AOIa[28:55]


AO

AO.ed<-AO[-5]

AO.ed
##################################################################################################################################################################################
######################################################### Factors for models##################################################################################################
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels

################################ Residual series ################
r
r6
################################ SST ############################

April.anom #28 entries 1978-
May.anom	#28 entries 1978-
June.anom	#28 entries 1978-
July.anom	#28 entries 1978-
April.May.anom	#28 entries 1978-
May.June.anom	#28 entries 1978-
June.July.anom	#28 entries 1978-
April.July.anom	#28 entries 1978-
May.July.anom	#28 entries 1978-

############################### NBT ##############################

NBT.lag.4	#27 enties beginning w/ 1978
NBT.lag.3	#28 entries beginning w/1978
NBT.lag.2 	#27 entries beginning w/1979
NBT.RA.3a	#28 entries
NBT.RA.3b	#27 entries
NBT.RA.2a	#28 entries
NBT.RA.2b	#27 entries

############################## Large males ######################

lag.3.males	#28 entries
lag.2.males	#28 entries
lag.1.males #28 entries
RA.3.males #28 entries
RA.2a.males	#28 entries
RA.2b.males #28 entries

############################## Flathead sole #######################

FHS.lag.3	#28 entries
FHS.lag.2	#28 entries
FHS.lag.1	#28 entries

############################## Wind vector components ##############

D.90.May
D.90.June
D.90.July
D.90.May.June
D.90.JJ
D.90.MJ

D.75.May
D.75.June
D.75.July
D.75.May.June
D.75.JJ
D.75.MJ

D.60.May
D.60.June
D.60.July
D.60.May.June
D.60.JJ
D.60.MJ

D.45.May
D.45.June
D.45.July
D.45.May.June
D.45.JJ
D.45.MJ

D.30.May
D.30.June
D.30.July
D.30.May.June
D.30.JJ
D.30.MJ

D.15.May
D.15.June
D.15.July
D.15.May.June
D.15.JJ
D.15.MJ

D.0.May
D.0.June
D.0.July
D.0.May.June
D.0.JJ
D.0.MJ

D.neg.15.May
D.neg.15.June
D.neg.15.July
D.neg.15.May.June
D.neg.15.JJ
D.neg.15.MJ

D.neg.30.May
D.neg.30.June
D.neg.30.July
D.neg.30.May.June
D.neg.30.JJ
D.neg.30.MJ

D.neg.45.May
D.neg.45.June
D.neg.45.July
D.neg.45.May.June
D.neg.45.JJ
D.neg.45.MJ

D.neg.60.May
D.neg.60.June
D.neg.60.July
D.neg.60.May.June
D.neg.60.JJ
D.neg.60.MJ


D.neg.75.May
D.neg.75.June
D.neg.75.July
D.neg.75.May.June
D.neg.75.JJ
D.neg.75.MJ

D.neg.90.May
D.neg.90.June
D.neg.90.July
D.neg.90.May.June
D.neg.90.JJ
D.neg.90.MJ
########################################################################################
############################### potential variables ####################################
May.anom	#28 entries 1978-
June.anom	#28 entries 1978-
May.June.anom	#28 entries 1978-

FHS.lag.3	#28 entries
FHS.lag.2	#28 entries
FHS.lag.1	#28 entries

D.neg.60.May.June
D.neg.90.May.June
D.90.May.June
#############################################################################################################################
#############################################################################################################################

par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
Year<-c(1978:2005)
#############################################################################################################################
################################# Fit multivariate models models ###########################################################
r


############################################################################################################################
###################### Minus 90 degree component and lag 1 FHS with interaction term  #######################################

fit.1aa<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.1aa)
D<-resid(fit.1aa)
plot(D~Year, ylab="EBS reduced model residuals",xlab="Release Year",pch=16)

AIC(fit.1aa)
AICc(fit.1aa)
aic.1<-AIC(fit.1aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.1a<-gls(r~D.neg.90.May.June,correlation=corAR1(),method="ML")
summary(fit.1a)
D<-resid(fit.1a)
plot(D)
AIC(fit.1a)
AICc(fit.1a)
aic.1a<-AIC(fit.1a)

############################################################################################################################
################### Drop minus 90 degree and add lag 1 FHS #################################################################

fit.1b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.1b)
D<-resid(fit.1b)
plot(D)
AIC(fit.1b)
AICc(fit.1b)
aic.1b<-AIC(fit.1b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.1bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.1bb)
D<-resid(fit.1bb)
plot(D)
AIC(fit.1bb)
AICc(fit.1bb)
aic.1b<-AIC(fit.1bb)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.1c<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.1c)
D<-resid(fit.1c)
plot(D)
AIC(fit.1c)
AICc(fit.1c)
aic.1c<-AIC(fit.1c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 90 degree component and lag 2 FHS ####################################

fit.1d<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.1d)
D<-resid(fit.1d)
plot(D)

AIC(fit.1d)
AICc(fit.1d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.1e<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.1e)
D<-resid(fit.1e)
plot(D)

AIC(fit.1e)
AICc(fit.1e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.1ee<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.1ee)
D<-resid(fit.1ee)
plot(D)

AIC(fit.1ee)
AICc(fit.1ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.1e3<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.1e3)
D<-resid(fit.1e3)
plot(D)

AIC(fit.1e3)
AICc(fit.1e3)

###############################################################################################################################
##################### swap lag 2 FHS for lag 1 FHS  #################################################################################

fit.1e4<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.1e4)
D<-resid(fit.1e4)
plot(D)

AIC(fit.1e4)
AICc(fit.1e4)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.1f<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+AO*D.neg.90.May.June,correlation=corAR1(),method="ML")
summary(fit.1f)
D<-resid(fit.1f)
plot(D)

AIC(fit.1f)
AICc(fit.1f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component AND swap in FHS lag 1 ###########################

fit.1ff<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2)+AO+AO*D.neg.90.May.June,correlation=corAR1(),method="ML")
summary(fit.1ff)
D<-resid(fit.1ff)
plot(D)

AIC(fit.1ff)
AICc(fit.1ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.1g<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.1g)
D<-resid(fit.1g)
plot(D)

AIC(fit.1g)
AICc(fit.1g)

lm.1g<-glm(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)

##################################################################################################################################
##################### Add Summer PDO and swap in lag 1 FHS #############################################################################################

fit.1gg<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.1gg)
D<-resid(fit.1gg)
plot(D)

AIC(fit.1gg)
AICc(fit.1gg)

lm.1g<-glm(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 90 degree component ########################################################

fit.1h<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.90.May.June,correlation=corAR1(),method="ML")
summary(fit.1h)
D<-resid(fit.1h)
plot(D)

AIC(fit.1h)
AICc(fit.1h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.1i<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.90.May.June+AO*D.neg.90.May.June,correlation=corAR1(),method="ML")
summary(fit.1i)
D<-resid(fit.1i)
plot(D)
AIC(fit.1i)
AICc(fit.1i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.1j<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.1j)
D<-resid(fit.1j)
plot(D)
AIC(fit.1j)
AICc(fit.1j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.1jj<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.1jj)
D<-resid(fit.1jj)
plot(D)
AIC(fit.1jj)
AICc(fit.1jj)

################################################################################################################################
#####################  Drop all interaction terms and pdos and swap in lag 1 FHS ##############################################################################

fit.1jjj<-gls(r~D.neg.90.May.June+I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.1jjj)
D<-resid(fit.1jjj)
plot(D)
AIC(fit.1jjj)
AICc(fit.1jjj)

################################################################################################################################
#####################  Drop all interaction terms and pdos and neg.90 and swap in lag 1 FHS ##############################################################################

fit.1jjjj<-gls(r~I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.1jjjj)
D<-resid(fit.1jjjj)
plot(D)
AIC(fit.1jjjj)
AICc(fit.1jjjj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.1k<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.1k)
D<-resid(fit.1k)
plot(D)

AIC(fit.1k)
AICc(fit.1k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.1m<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.neg.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.1m)
D<-resid(fit.1m)
plot(D)

AIC(fit.1m)
AICc(fit.1m)

lm.5m<-glm(r~D.neg.90.May.June+I(FHS.lag.1^2)+May.June.anom+AO+pdos+D.neg.90.May.June*I(FHS.lag.1^2)+May.anom+lag.2.males+lag.1.males)
step(lm.5m)

summary(lm.5m)


###########################################################################################################################################
##################### Exchange neg 75 for neg 90 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### Minus 75 degree component and lag 1 FHS with interaction term  #######################################

fit.2aa<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+D.neg.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.2aa)
D<-resid(fit.2aa)
plot(D)

AIC(fit.2aa)
AICc(fit.2aa)
aic.1<-AIC(fit.2aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.2a<-gls(r~D.neg.75.May.June,correlation=corAR1(),method="ML")
summary(fit.2a)
D<-resid(fit.2a)
plot(D)
AIC(fit.2a)
AICc(fit.2a)
aic.1a<-AIC(fit.2a)

############################################################################################################################
################### Drop minus 75 degree and add lag 1 FHS #################################################################

fit.2b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.2b)
D<-resid(fit.2b)
plot(D)
AIC(fit.2b)
AICc(fit.2b)
aic.1b<-AIC(fit.2b)

############################################################################################################################
################### Drop minus 75 degree and add AO #################################################################

fit.2ba<-gls(r~AO,correlation=corAR1(),method="ML")
summary(fit.2ba)
D<-resid(fit.2ba)
plot(D)
AIC(fit.2ba)
AICc(fit.2ba)
aic.1b<-AIC(fit.2ba)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.2bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.2bb)
D<-resid(fit.2bb)
plot(D)
AIC(fit.2bb)
AICc(fit.2bb)
aic.1b<-AIC(fit.2bb)
############################################################################################################################
################### Add minus 75 degree component #########################################################################

fit.2c<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.2c)
D<-resid(fit.2c)
plot(D)
AIC(fit.2c)
AICc(fit.2c)

aic.1c<-AIC(fit.1c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 75 degree component and lag 2 FHS ####################################

fit.2d<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.2d)
D<-resid(fit.2d)
plot(D)

AIC(fit.2d)
AICc(fit.2d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.2e<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.2e)
D<-resid(fit.2e)
plot(D)

AIC(fit.2e)
AICc(fit.2e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.2ee<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+D.neg.75.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.2ee)
D<-resid(fit.2ee)
plot(D)

AIC(fit.2ee)
AICc(fit.2ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.2e3<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+D.neg.75.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.2e3)
D<-resid(fit.2e3)
plot(D)

AIC(fit.2e3)
AICc(fit.2e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 75 degree component ###########################

fit.2f<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+AO*D.neg.75.May.June,correlation=corAR1(),method="ML")
summary(fit.2f)
D<-resid(fit.2f)
plot(D)

AIC(fit.2f)
AICc(fit.2f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 75 degree component AND swap in FHS lag 1 ###########################

fit.2ff<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+D.neg.75.May.June*I(FHS.lag.1^2)+AO+AO*D.neg.75.May.June,correlation=corAR1(),method="ML")
summary(fit.2ff)
D<-resid(fit.2ff)
plot(D)

AIC(fit.2ff)
AICc(fit.2ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.2g<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.2g)
D<-resid(fit.2g)
plot(D)

AIC(fit.2g)
AICc(fit.2g)

lm.2g<-glm(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.2g)

##################################################################################################################################
##################### Add Summer PDO and swap in lag 1 FHS#############################################################################################

fit.2gg<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+D.neg.75.May.June*I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.2gg)
D<-resid(fit.2gg)
plot(D)

AIC(fit.2gg)
AICc(fit.2gg)

lm.2g<-glm(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.2g)
###############################################################################################################################
##################### Add interaction term between Summer PDO and minus 75 degree component ########################################################

fit.2h<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.75.May.June,correlation=corAR1(),method="ML")
summary(fit.2h)
D<-resid(fit.2h)
plot(D)

AIC(fit.2h)
AICc(fit.2h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.2i<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.75.May.June+AO*D.neg.75.May.June,correlation=corAR1(),method="ML")
summary(fit.2i)
D<-resid(fit.2i)
plot(D)
AIC(fit.2i)
AICc(fit.2i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.2j<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.2j)
D<-resid(fit.2j)
plot(D)
AIC(fit.2j)
AICc(fit.2j)

################################################################################################################################
#####################  Drop all interaction terms and swap in FHS lag 1 ##############################################################################

fit.2jj<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.2jj)
D<-resid(fit.2jj)
plot(D)
AIC(fit.2jj)
AICc(fit.2jj)


################################################################################################################################
#####################  Drop all interaction terms AND PDOS and swap in FHS lag 1  ##############################################################################

fit.2jjj<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.2jjj)
D<-resid(fit.2jjj)
plot(D)
AIC(fit.2jjj)
AICc(fit.2jjj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.2k<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.2k)
D<-resid(fit.2k)
plot(D)

AIC(fit.2k)
AICc(fit.2k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.2m<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.neg.75.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.2m)
D<-resid(fit.2m)
plot(D)

AIC(fit.2m)
AICc(fit.2m)

################################################################################################################################
#####################  Run step wise model #########################################################################################

fit.2m<-gls(r~D.neg.75.May.June+I(FHS.lag.1^2)+May.June.anom+AO+pdos+D.neg.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.2m)
D<-resid(fit.2m)
plot(D)

AIC(fit.2m)
AICc(fit.2m)

lm.5m<-glm(r~D.neg.75.May.June+I(FHS.lag.1^2)+May.June.anom+AO+pdos+D.neg.75.May.June*I(FHS.lag.1^2)+May.anom+lag.2.males+lag.1.males)
step(lm.5m)

###########################################################################################################################################
##################### Exchange neg 75 for neg 60 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### Minus 60 degree component and lag 1 FHS with interaction term  #######################################

fit.3aa<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)

AIC(fit.3aa)
AICc(fit.3aa)
aic.1<-AIC(fit.3aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.3a<-gls(r~D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.3a)
D<-resid(fit.3a)
plot(D)
AIC(fit.3a)
AICc(fit.3a)
aic.1a<-AIC(fit.3a)

############################################################################################################################
################### Drop neg 60 and add AO #########################################################################################

fit.3aa<-gls(r~AO,correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)
AIC(fit.3aa)
AICc(fit.3aa)
aic.1a<-AIC(fit.3aa)
############################################################################################################################
################### Drop minus 75 degree and add lag 1 FHS #################################################################

fit.3b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.3b)
D<-resid(fit.3b)
plot(D)
AIC(fit.3b)
AICc(fit.3b)
aic.1b<-AIC(fit.3b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.3bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.3bb)
D<-resid(fit.3bb)
plot(D)
AIC(fit.3bb)
AICc(fit.3bb)
aic.1b<-AIC(fit.3bb)
############################################################################################################################
################### Add minus 60 degree component #########################################################################

fit.3c<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.3c)
D<-resid(fit.3c)
plot(D)
AIC(fit.3c)
AICc(fit.3c)

aic.1c<-AIC(fit.3c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 75 degree component and lag 2 FHS ####################################

fit.3d<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.3d)
D<-resid(fit.3d)
plot(D)

AIC(fit.3d)
AICc(fit.3d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.3e<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.3e)
D<-resid(fit.3e)
plot(D)

AIC(fit.3e)
AICc(fit.3e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.3ee<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.3ee)
D<-resid(fit.3ee)
plot(D)

AIC(fit.3ee)
AICc(fit.3ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.3e3<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.3e3)
D<-resid(fit.3e3)
plot(D)

AIC(fit.3e3)
AICc(fit.3e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 60 degree component ###########################

fit.3f<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+AO*D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.3f)
D<-resid(fit.3f)
plot(D)

AIC(fit.3f)
AICc(fit.3f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 60 degree component AND swap in FHS lag 1 ###########################

fit.3ff<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^2)+AO+AO*D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.3ff)
D<-resid(fit.3ff)
plot(D)

AIC(fit.3ff)
AICc(fit.3ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.3g<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.3g)
D<-resid(fit.3g)
plot(D)

AIC(fit.3g)
AICc(fit.3g)

lm.2g<-glm(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.2g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 60 degree component ########################################################

fit.3h<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.3h)
D<-resid(fit.3h)
plot(D)

AIC(fit.3h)
AICc(fit.3h)

###############################################################################################################################
#################### Reinstate interaction term between minus 60 and AO #######################################################
fit.3i<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.neg.60.May.June+AO*D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.3i)
D<-resid(fit.3i)
plot(D)
AIC(fit.3i)
AICc(fit.3i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.3j<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.3j)
D<-resid(fit.3j)
plot(D)
AIC(fit.3j)
AICc(fit.3j)

################################################################################################################################
##################### AS EXPERIMENT: Drop pdos ##############################################################################

fit.3jj<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.3jj)
D<-resid(fit.3jj)
plot(D)
AIC(fit.3jj)
AICc(fit.3jj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.3k<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.3k)
D<-resid(fit.3k)
plot(D)

AIC(fit.3k)
AICc(fit.3k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.3m<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.neg.60.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.3m)
D<-resid(fit.3m)
plot(D)

AIC(fit.3m)
AICc(fit.3m)



###########################################################################################################################################
##################### Exchange POSITIVE 90 for neg 60 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### Minus 60 degree component and lag 1 FHS with interaction term  #######################################

fit.4aa<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.4aa)
D<-resid(fit.4aa)
plot(D)

AIC(fit.4aa)
AICc(fit.4aa)
aic.1<-AIC(fit.4aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.4a<-gls(r~D.90.May.June,correlation=corAR1(),method="ML")
summary(fit.4a)
D<-resid(fit.4a)
plot(D)
AIC(fit.4a)
AICc(fit.4a)
aic.1a<-AIC(fit.4a)

############################################################################################################################
################### Drop 90 degree and add lag 1 FHS #################################################################

fit.4b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.4b)
D<-resid(fit.4b)
plot(D)
AIC(fit.4b)
AICc(fit.4b)
aic.1b<-AIC(fit.3b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.4bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4bb)
D<-resid(fit.4bb)
plot(D)
AIC(fit.4bb)
AICc(fit.4bb)
aic.1b<-AIC(fit.4bb)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.4c<-gls(r~D.90.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4c)
D<-resid(fit.4c)
plot(D)
AIC(fit.4c)
AICc(fit.4c)

aic.1c<-AIC(fit.4c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between 90 degree component and lag 2 FHS ####################################

fit.4d<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4d)
D<-resid(fit.4d)
plot(D)

AIC(fit.4d)
AICc(fit.4d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.4e<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.4e)
D<-resid(fit.4e)
plot(D)

AIC(fit.4e)
AICc(fit.4e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.4ee<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.4ee)
D<-resid(fit.4ee)
plot(D)

AIC(fit.4ee)
AICc(fit.4ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.4e3<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.4e3)
D<-resid(fit.4e3)
plot(D)

AIC(fit.4e3)
AICc(fit.4e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.4f<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2)+AO+AO*D.90.May.June,correlation=corAR1(),method="ML")
summary(fit.4f)
D<-resid(fit.4f)
plot(D)

AIC(fit.4f)
AICc(fit.4f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and 90 degree component AND swap in FHS lag 1 ###########################

fit.4ff<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2)+AO+AO*D.90.May.June,correlation=corAR1(),method="ML")
summary(fit.4ff)
D<-resid(fit.4ff)
plot(D)

AIC(fit.4ff)
AICc(fit.4ff)

##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.4g<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4g)
D<-resid(fit.4g)
plot(D)

AIC(fit.4g)
AICc(fit.4g)

lm.2g<-glm(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.2g)

##################################################################################################################################
##################### Add Summer PDO and swap in lag 1 FHS #############################################################################################

fit.4gg<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4gg)
D<-resid(fit.4gg)
plot(D)

AIC(fit.4gg)
AICc(fit.4gg)

lm.2g<-glm(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.2g)
###############################################################################################################################
##################### Add interaction term between Summer PDO and 90 degree component ########################################################

fit.4h<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.90.May.June,correlation=corAR1(),method="ML")
summary(fit.4h)
D<-resid(fit.4h)
plot(D)

AIC(fit.4h)
AICc(fit.4h)

###############################################################################################################################
#################### Reinstate interaction term between 90 and AO #######################################################
fit.4i<-gls(r~D.90.May.June+I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.90.May.June+AO*D.90.May.June,correlation=corAR1(),method="ML")
summary(fit.4i)
D<-resid(fit.4i)
plot(D)
AIC(fit.4i)
AICc(fit.4i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.4j<-gls(r~D.90.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4j)
D<-resid(fit.4j)
plot(D)
AIC(fit.4j)
AICc(fit.4j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.4jj<-gls(r~D.90.May.June+I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4jj)
D<-resid(fit.4jj)
plot(D)
AIC(fit.4jj)
AICc(fit.4jj)

###############################################################################################################################
##################### Drop all interaction terms and PDOS and swap in lag 1 FHS #########################################

fit.4jjj<-gls(r~D.90.May.June+I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.4jjj)
D<-resid(fit.4jjj)
plot(D)

AIC(fit.4jjj)
AICc(fit.4jjj)
###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.4k<-gls(r~D.90.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4k)
D<-resid(fit.4k)
plot(D)

AIC(fit.4k)
AICc(fit.4k)

################################################################################################################################
##################### Add interaction term for FHS and May-June 90 #########################################################################################

fit.4m<-gls(r~D.90.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4m)
D<-resid(fit.4m)
plot(D)

AIC(fit.4m)
AICc(fit.4m)

###############################################################################################################################
###################### Add lag 2 males #########################################################################################

fit.4n<-gls(r~D.90.May.June+I(FHS.lag.2^2)+lag.2.males+May.June.anom+AO+pdos+D.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4n)
D<-resid(fit.4n)
plot(D)

AIC(fit.4n)
AICc(fit.4n)
###############################################################################################################################
###################### TEST #########################################################################################

fit.4nn<-gls(r~D.90.May.June+I(FHS.lag.2^2)+lag.2.males+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.4nn)
D<-resid(fit.4nn)
plot(D)

AIC(fit.4nn)
AICc(fit.4nn)

################################################################################################################################
##################### Add RA3a NBT #########################################################################################

fit.4o<-gls(r~D.90.May.June+I(FHS.lag.2^2)+NBT.RA.3a+May.June.anom+AO+pdos+D.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4o)
D<-resid(fit.4o)
plot(D)

AIC(fit.4o)
AICc(fit.4o)

################################################################################################################################
##################### TEST #########################################################################################

fit.4oo<-gls(r~D.90.May.June+I(FHS.lag.1^2)+NBT.RA.3b+AO+pdos+D.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.4oo)
D<-resid(fit.4oo)
plot(D)

AIC(fit.4oo)
AICc(fit.4oo)

################################################################################################################################
##################### Add lag 2 NBT #########################################################################################

fit.4p<-gls(r~D.90.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.90.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.4p)
D<-resid(fit.4p)
plot(D)

AIC(fit.4p)
AICc(fit.4p)


############################################################################################################################
###################### Drop 90 degrees and use POSITIVE 75 degrees ###################################################
##############################################################################################################################

fit.5aa<-gls(r~D.75.May.June+I(FHS.lag.1^2)+D.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5aa)
D<-resid(fit.5aa)
plot(D)

AIC(fit.5aa)
AICc(fit.5aa)
aic.5<-AIC(fit.5aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.5a<-gls(r~D.75.May.June,correlation=corAR1(),method="ML")
summary(fit.5a)
D<-resid(fit.5a)
plot(D)
AIC(fit.5a)
AICc(fit.5a)
aic.5a<-AIC(fit.5a)

############################################################################################################################
################### Drop minus 90 degree and add lag 1 FHS #################################################################

fit.5b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5b)
D<-resid(fit.5b)
plot(D)
AIC(fit.5b)
AICc(fit.5b)
aic.5b<-AIC(fit.5b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.5bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5bb)
D<-resid(fit.5bb)
plot(D)
AIC(fit.5bb)
AICc(fit.5bb)
aic.5bb<-AIC(fit.5bb)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.5c<-gls(r~D.75.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5c)
D<-resid(fit.5c)
plot(D)
AIC(fit.5c)
AICc(fit.5c)
aic.1c<-AIC(fit.5c)


model.5.aic<-c(aic.5,aic.5a,aic.5b,aic.5c)
model.5.aic

##################################################################################################################################
##################### Add interaction term between minus 90 degree component and lag 2 FHS ####################################

fit.5d<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5d)
D<-resid(fit.5d)
plot(D)

AIC(fit.5d)
AICc(fit.5d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.5e<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.5e)
D<-resid(fit.5e)
plot(D)

AIC(fit.5e)
AICc(fit.5e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.5ee<-gls(r~D.75.May.June+I(FHS.lag.1^2)+D.75.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.5ee)
D<-resid(fit.5ee)
plot(D)

AIC(fit.5ee)
AICc(fit.5ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.5e3<-gls(r~D.75.May.June+I(FHS.lag.1^2)+D.75.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5e3)
D<-resid(fit.5e3)
plot(D)

AIC(fit.5e3)
AICc(fit.5e3)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.5f<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+AO*D.75.May.June,correlation=corAR1(),method="ML")
summary(fit.5f)
D<-resid(fit.5f)
plot(D)

AIC(fit.5f)
AICc(fit.5f)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component AND swap in FHS lag 1 ###########################

fit.5ff<-gls(r~D.75.May.June+I(FHS.lag.1^2)+D.75.May.June*I(FHS.lag.1^2)+AO+AO*D.75.May.June,correlation=corAR1(),method="ML")
summary(fit.5ff)
D<-resid(fit.5ff)
plot(D)

AIC(fit.5ff)
AICc(fit.5ff)

##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.5g<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5g)
D<-resid(fit.5g)
plot(D)

AIC(fit.5g)
AICc(fit.5g)

lm.1g<-glm(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 90 degree component ########################################################

fit.5h<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.75.May.June,correlation=corAR1(),method="ML")
summary(fit.5h)
D<-resid(fit.5h)
plot(D)

AIC(fit.5h)
AICc(fit.5h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.5i<-gls(r~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.75.May.June+AO*D.75.May.June,correlation=corAR1(),method="ML")
summary(fit.5i)
D<-resid(fit.5i)
plot(D)
AIC(fit.5i)
AICc(fit.5i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.5j<-gls(r~D.75.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5j)
D<-resid(fit.5j)
plot(D)
AIC(fit.5j)
AICc(fit.5j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.5jj<-gls(r~D.75.May.June+I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5jj)
D<-resid(fit.5jj)
plot(D)
AIC(fit.5jj)
AICc(fit.5jj)


###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.5k<-gls(r~D.75.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5k)
D<-resid(fit.5k)
plot(D)

AIC(fit.5k)
AICc(fit.5k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.5m<-gls(r~D.75.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.75.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5m)
D<-resid(fit.5m)
plot(D)

AIC(fit.5m)
AICc(fit.5m)




############################################################################################################################
###################### Drop 75 degrees and use POSITIVE 60 degrees ###################################################
##############################################################################################################################

fit.5aa<-gls(r~D.60.May.June+I(FHS.lag.1^2)+D.60.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5aa)
D<-resid(fit.5aa)
plot(D)

AIC(fit.5aa)
AICc(fit.5aa)
aic.5<-AIC(fit.5aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.5a<-gls(r~D.60.May.June,correlation=corAR1(),method="ML")
summary(fit.5a)
D<-resid(fit.5a)
plot(D)
AIC(fit.5a)
AICc(fit.5a)
aic.5a<-AIC(fit.5a)

############################################################################################################################
################### Drop minus 90 degree and add lag 1 FHS #################################################################

fit.5b<-gls(r~I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5b)
D<-resid(fit.5b)
plot(D)
AIC(fit.5b)
AICc(fit.5b)
aic.5b<-AIC(fit.5b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.5bb<-gls(r~I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5bb)
D<-resid(fit.5bb)
plot(D)
AIC(fit.5bb)
AICc(fit.5bb)
aic.5bb<-AIC(fit.5bb)

############################################################################################################################
################### Drop lag 1 FHS and add AO #################################################################

fit.5ba<-gls(r~AO,correlation=corAR1(),method="ML")
summary(fit.5ba)
D<-resid(fit.5ba)
plot(D)
AIC(fit.5ba)
AICc(fit.5ba)
aic.5bb<-AIC(fit.5ba)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.5c<-gls(r~D.60.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5c)
D<-resid(fit.5c)
plot(D)
AIC(fit.5c)
AICc(fit.5c)
aic.1c<-AIC(fit.5c)


model.5.aic<-c(aic.5,aic.5a,aic.5b,aic.5c)
model.5.aic

##################################################################################################################################
##################### Add interaction term between minus 90 degree component and lag 2 FHS ####################################

fit.5d<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5d)
D<-resid(fit.5d)
plot(D)

AIC(fit.5d)
AICc(fit.5d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.5e<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO,correlation=corAR1(),method="ML")
summary(fit.5e)
D<-resid(fit.5e)
plot(D)

AIC(fit.5e)
AICc(fit.5e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.5ee<-gls(r~D.60.May.June+I(FHS.lag.1^2)+D.60.May.June*I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.5ee)
D<-resid(fit.5ee)
plot(D)

AIC(fit.5ee)
AICc(fit.5ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.5e3<-gls(r~D.60.May.June+I(FHS.lag.1^2)+D.60.May.June*I(FHS.lag.1^2)+AO+I(FHS.lag.2^2)*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5e3)
D<-resid(fit.5e3)
plot(D)

AIC(fit.5e3)
AICc(fit.5e3)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.5f<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO+AO*D.60.May.June,correlation=corAR1(),method="ML")
summary(fit.5f)
D<-resid(fit.5f)
plot(D)

AIC(fit.5f)
AICc(fit.5f)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component AND swap in FHS lag 1 ###########################

fit.5ff<-gls(r~D.60.May.June+I(FHS.lag.1^2)+D.60.May.June*I(FHS.lag.1^2)+AO+AO*D.60.May.June,correlation=corAR1(),method="ML")
summary(fit.5ff)
D<-resid(fit.5ff)
plot(D)

AIC(fit.5ff)
AICc(fit.5ff)

##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.5g<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5g)
D<-resid(fit.5g)
plot(D)

AIC(fit.5g)
AICc(fit.5g)

lm.1g<-glm(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 90 degree component ########################################################

fit.5h<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.60.May.June,correlation=corAR1(),method="ML")
summary(fit.5h)
D<-resid(fit.5h)
plot(D)

AIC(fit.5h)
AICc(fit.5h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.5i<-gls(r~D.60.May.June+I(FHS.lag.2^2)+D.60.May.June*I(FHS.lag.2^2)+AO+pdos+pdos*D.60.May.June+AO*D.60.May.June,correlation=corAR1(),method="ML")
summary(fit.5i)
D<-resid(fit.5i)
plot(D)
AIC(fit.5i)
AICc(fit.5i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.5j<-gls(r~D.60.May.June+I(FHS.lag.2^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5j)
D<-resid(fit.5j)
plot(D)
AIC(fit.5j)
AICc(fit.5j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.5jj<-gls(r~D.60.May.June+I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5jj)
D<-resid(fit.5jj)
plot(D)
AIC(fit.5jj)
AICc(fit.5jj)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS and drop 60 degree ##############################################################################

fit.5jjj<-gls(r~I(FHS.lag.1^2)+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5jjj)
D<-resid(fit.5jjj)
plot(D)
AIC(fit.5jjj)
AICc(fit.5jjj)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS and drop 60 degree and pdos ##############################################################################

fit.5jjjj<-gls(r~I(FHS.lag.1^2)+AO,correlation=corAR1(),method="ML")
summary(fit.5jjjj)
D<-resid(fit.5jjjj)
plot(D)
AIC(fit.5jjjj)
AICc(fit.5jjjj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.5k<-gls(r~D.60.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos,correlation=corAR1(),method="ML")
summary(fit.5k)
D<-resid(fit.5k)
plot(D)

AIC(fit.5k)
AICc(fit.5k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.5m<-gls(r~D.60.May.June+I(FHS.lag.2^2)+May.June.anom+AO+pdos+D.60.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.5m)
D<-resid(fit.5m)
plot(D)

AIC(fit.5m)
AICc(fit.5m)






############################################################################################################################################
################################ Re-run using editted residual series ######################################################################
##########################################################################################################################################

r6


AO.ed

pdos.ed

FHS.lag.1.ed
FHS.lag.2.ed
FHS.lag.3.ed

lag.3.males.ed	
lag.2.males.ed	
lag.1.males.ed
RA.3.males.ed
RA.2a.males.ed	
RA.2b.males.ed


NBT.lag.4.ed
NBT.lag.3.ed
NBT.lag.2.ed
NBT.RA.3a.ed
NBT.RA.3b.ed
NBT.RA.2a.ed
NBT.RA.2b.ed

May.ed
June.ed
July.ed
April.May.ed
May.June.ed
June.July.ed
April.July.ed
May.July.ed


D.90.May.ed
D.90.June.ed
D.90.July.ed
D.90.May.June.ed
D.90.JJ.ed
D.90.MJ.ed

D.75.May.ed
D.75.June.ed
D.75.July.ed
D.75.May.June.ed
D.75.JJ.ed
D.75.MJ.ed


D.60.May.ed
D.60.June.ed
D.60.July.ed
D.60.May.June.ed
D.60.JJ.ed
D.60.MJ.ed

D.45.May.ed
D.45.June.ed
D.45.July.ed
D.45.May.June.ed
D.45.JJ.ed
D.45.MJ.ed

D.30.May.ed
D.30.June.ed
D.30.July.ed
D.30.May.June.ed
D.30.JJ.ed
D.30.MJ.ed

D.15.May.ed
D.15.June.ed
D.15.July.ed
D.15.May.June.ed
D.15.JJ.ed
D.15.MJ.ed

D.0.May.ed 
D.0.June.ed
D.0.July.ed
D.0.May.June.ed
D.0.JJ.ed
D.0.MJ.ed

D.neg.15.May.ed
D.neg.15.June.ed
D.neg.15.July.ed
D.neg.15.May.June.ed
D.neg.15.JJ.ed
D.neg.15.MJ.ed

D.neg.30.May.ed
D.neg.30.June.ed
D.neg.30.July.ed
D.neg.30.May.June.ed
D.neg.30.JJ.ed
D.neg.30.MJ.ed


D.neg.45.May.ed
D.neg.45.June.ed
D.neg.45.July.ed
D.neg.45.May.June.ed
D.neg.45.JJ.ed
D.neg.45.MJ.ed

D.neg.60.May.ed
D.neg.60.June.ed
D.neg.60.July.ed
D.neg.60.May.June.ed
D.neg.60.JJ.ed
D.neg.60.MJ.ed


D.neg.75.May.ed
D.neg.75.June.ed
D.neg.75.July.ed
D.neg.75.May.June.ed
D.neg.75.JJ.ed
D.neg.75.MJ.ed

D.neg.90.May.ed
D.neg.90.June.ed
D.neg.90.July.ed
D.neg.90.May.June.ed
D.neg.90.JJ.ed
D.neg.90.MJ.ed




Year.ed<-Year[-5]



############################################################################################################################
###################### Minus 90 degree component and lag 1 FHS with interaction term  #######################################

fit.1aa<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.1aa)
D<-resid(fit.1aa)
plot(D)

AIC(fit.1aa)
AICc(fit.1aa)
aic.1<-AIC(fit.1aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.1a<-gls(r6~D.neg.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.1a)
D<-resid(fit.1a)
plot(D)
AIC(fit.1a)
AICc(fit.1a)
aic.1a<-AIC(fit.1a)

############################################################################################################################
################### Drop minus 90 degree and add lag 1 FHS #################################################################

##################### BEST REDUCED MODEL FOR EBS group #####################################################

par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.25,cex=1.25)
fit.1b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.1b)
D<-resid(fit.1b)
plot(D~Year.ed,ylab="Reduced model residuals",xlab="Hatch Year",pch=16)
abline(h=0,lty=2)
AIC(fit.1b)
AICc(fit.1b)
aic.1b<-AIC(fit.1b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.1bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.1bb)
D<-resid(fit.1bb)
plot(D)
AIC(fit.1bb)
AICc(fit.1bb)
aic.1b<-AIC(fit.1bb)

############################################################################################################################
################### Drop lag 1 FHS and add AO #################################################################

fit.1ba<-gls(r6~AO.ed,correlation=corAR1(),method="ML")
summary(fit.1ba)
D<-resid(fit.1ba)
plot(D)
AIC(fit.1ba)
AICc(fit.1ba)
aic.1b<-AIC(fit.1ba)

############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.1c<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.1c)
D<-resid(fit.1c)
plot(D)
AIC(fit.1c)
AICc(fit.1c)
aic.1c<-AIC(fit.1c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

############################################################################################################################
################### Add minus 90 degree component WITH lag 1 FHS AS EXPERIMENT--NOT RECORDED #########################################################################

fit.1cc<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.1cc)
D<-resid(fit.1cc)
plot(D)
AIC(fit.1cc)
AICc(fit.1cc)
aic.1cc<-AIC(fit.1cc)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 90 degree component and lag 2 FHS ####################################

fit.1d<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.1d)
D<-resid(fit.1d)
plot(D)

AIC(fit.1d)
AICc(fit.1d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.1e<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.1e)
D<-resid(fit.1e)
plot(D)

AIC(fit.1e)
AICc(fit.1e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.1ee<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.1ee)
D<-resid(fit.1ee)
plot(D)

AIC(fit.1ee)
AICc(fit.1ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.1e3<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.1e3)
D<-resid(fit.1e3)
plot(D)

AIC(fit.1e3)
AICc(fit.1e3)

###############################################################################################################################
##################### swap lag 2 FHS for lag 1 FHS  #################################################################################

fit.1e4<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.1e4)
D<-resid(fit.1e4)
plot(D)

AIC(fit.1e4)
AICc(fit.1e4)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.1f<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.neg.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.1f)
D<-resid(fit.1f)
plot(D)

AIC(fit.1f)
AICc(fit.1f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component AND swap in FHS lag 1 ###########################

fit.1ff<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.neg.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.1ff)
D<-resid(fit.1ff)
plot(D)

AIC(fit.1ff)
AICc(fit.1ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.1g<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.1g)
D<-resid(fit.1g)
plot(D)

AIC(fit.1g)
AICc(fit.1g)

lm.1g<-glm(r6~D.neg.90.May.June+I(FHS.lag.2^2)+D.neg.90.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 90 degree component ########################################################

fit.1h<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.1h)
D<-resid(fit.1h)
plot(D)

AIC(fit.1h)
AICc(fit.1h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.1i<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.90.May.June.ed+AO.ed*D.neg.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.1i)
D<-resid(fit.1i)
plot(D)
AIC(fit.1i)
AICc(fit.1i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.1j<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.1j)
D<-resid(fit.1j)
plot(D)
AIC(fit.1j)
AICc(fit.1j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

###################### NOTE: BEST EBS FULL MODEL ####################################################
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.25,cex=1.25) #configure axis labels
fit.1jj<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.1jj)
D<-resid(fit.1jj)
plot(D~Year.ed, ylab="Full model residuals",xlab="Hatch Year",pch=16)
abline(h=0,lty=2)
AIC(fit.1jj)
AICc(fit.1jj)


###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.1k<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.1k)
D<-resid(fit.1k)
plot(D)

AIC(fit.1k)
AICc(fit.1k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.1m<-gls(r6~D.neg.90.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.neg.90.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.1m)
D<-resid(fit.1m)
plot(D)

AIC(fit.1m)
AICc(fit.1m)


###########################################################################################################################################
##################### Exchange neg 75 for neg 90 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### Minus 75 degree component and lag 1 FHS with interaction term  #######################################

fit.2aa<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.2aa)
D<-resid(fit.2aa)
plot(D)

AIC(fit.2aa)
AICc(fit.2aa)
aic.1<-AIC(fit.2aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.2a<-gls(r6~D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2a)
D<-resid(fit.2a)
plot(D)
AIC(fit.2a)
AICc(fit.2a)
aic.1a<-AIC(fit.2a)

############################################################################################################################
################### Drop minus 75 degree and add lag 1 FHS #################################################################

fit.2b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.2b)
D<-resid(fit.2b)
plot(D)
AIC(fit.2b)
AICc(fit.2b)
aic.1b<-AIC(fit.2b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.2bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.2bb)
D<-resid(fit.2bb)
plot(D)
AIC(fit.2bb)
AICc(fit.2bb)
aic.1b<-AIC(fit.2bb)
############################################################################################################################
################### Add minus 75 degree component #########################################################################

fit.2c<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.2c)
D<-resid(fit.2c)
plot(D)
AIC(fit.2c)
AICc(fit.2c)

aic.1c<-AIC(fit.1c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 75 degree component and lag 2 FHS ####################################

fit.2d<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.2d)
D<-resid(fit.2d)
plot(D)

AIC(fit.2d)
AICc(fit.2d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.2e<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.2e)
D<-resid(fit.2e)
plot(D)

AIC(fit.2e)
AICc(fit.2e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.2ee<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.2ee)
D<-resid(fit.2ee)
plot(D)

AIC(fit.2ee)
AICc(fit.2ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.2e3<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.2e3)
D<-resid(fit.2e3)
plot(D)

AIC(fit.2e3)
AICc(fit.2e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 75 degree component ###########################

fit.2f<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2f)
D<-resid(fit.2f)
plot(D)

AIC(fit.2f)
AICc(fit.2f)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 75 degree component AND swap in FHS lag 1 ###########################

fit.2ff<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2ff)
D<-resid(fit.2ff)
plot(D)

AIC(fit.2ff)
AICc(fit.2ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.2g<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2g)
D<-resid(fit.2g)
plot(D)

AIC(fit.2g)
AICc(fit.2g)

lm.2g<-glm(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)

##################################################################################################################################
##################### Add Summer PDO AND swap lag 1 FHS for lag 2  #############################################################################################

fit.2gg<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2gg)
D<-resid(fit.2gg)
plot(D)

AIC(fit.2gg)
AICc(fit.2gg)

lm.2gg<-glm(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+pdos.ed)
step(lm.2gg)
###############################################################################################################################
##################### Add interaction term between Summer PDO and minus 75 degree component ########################################################

fit.2h<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2h)
D<-resid(fit.2h)
plot(D)

AIC(fit.2h)
AICc(fit.2h)

###############################################################################################################################
##################### AS EXPERIMENT Add interaction term between Summer PDO and minus 75 degree component and drop interaction between FHS lag 2 and neg 75 ########################################################

fit.2hh<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2hh)
D<-resid(fit.2hh)
plot(D)

AIC(fit.2hh)
AICc(fit.2hh)

###############################################################################################################################
##################### AS EXPERIMENT Add interaction term between Summer PDO and minus 75 degree component and drop interaction between FHS lag 2 and neg 75 and swap in lag 1 FHS########################################################

fit.2hhh<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2hhh)
D<-resid(fit.2hhh)
plot(D)

AIC(fit.2hhh)
AICc(fit.2hhh)

###############################################################################################################################
#################### Reinstate interaction term between minus 75 and AO #######################################################
fit.2i<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.75.May.June.ed+
              AO.ed*D.neg.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.2i)
D<-resid(fit.2i)
plot(D)
AIC(fit.2i)
AICc(fit.2i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.2j<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2j)
D<-resid(fit.2j)
plot(D)
AIC(fit.2j)
AICc(fit.2j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.2jj<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2jj)
D<-resid(fit.2jj)
plot(D)
AIC(fit.2jj)
AICc(fit.2jj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.2k<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2k)
D<-resid(fit.2k)
plot(D)

AIC(fit.2k)
AICc(fit.2k)

###############################################################################################################################
##################### Add May-June SST and swap in lag 1 FHS  #########################################################################################

fit.2kk<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2kk)
D<-resid(fit.2kk)
plot(D)

AIC(fit.2kk)
AICc(fit.2kk)
################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.2m<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.neg.75.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.2m)
D<-resid(fit.2m)
plot(D)

AIC(fit.2m)
AICc(fit.2m)

################################################################################################################################
#####################  Drop all interaction terms and AO and swap in lag 1 FHS ##############################################################################

fit.2n<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.1.ed^2)+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.2n)
D<-resid(fit.2n)
plot(D)
AIC(fit.2n)
AICc(fit.2n)

################################################################################################################################
#####################  Drop all interaction terms add lag 1 males ##############################################################################

fit.2o<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+lag.1.males.ed,correlation=corAR1(),method="ML")
summary(fit.2o)
D<-resid(fit.2o)
plot(D)
AIC(fit.2o)
AICc(fit.2o)

################################################################################################################################
#####################  Drop all interaction terms add lag 1 and 2 males ################################################################

fit.2p<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+lag.1.males.ed+lag.2.males.ed,correlation=corAR1(),method="ML")
summary(fit.2p)
D<-resid(fit.2p)
plot(D)
AIC(fit.2p)
AICc(fit.2p)
################################################################################################################################
#####################  Add interaction between lag 1 and 2 males ################################################################

fit.2q<-gls(r6~D.neg.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+lag.1.males.ed+lag.2.males.ed+lag.1.males.ed*lag.2.males.ed,correlation=corAR1(),method="ML")
summary(fit.2q)
D<-resid(fit.2q)
plot(D)
AIC(fit.2q)
AICc(fit.2q)


###########################################################################################################################################
##################### Exchange neg 60 for neg 75 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### Minus 60 degree component and lag 1 FHS with interaction term  #######################################

fit.3aa<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)

AIC(fit.3aa)
AICc(fit.3aa)
aic.1<-AIC(fit.3aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.3a<-gls(r6~D.neg.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3a)
D<-resid(fit.3a)
plot(D)
AIC(fit.3a)
AICc(fit.3a)
aic.1a<-AIC(fit.3a)

############################################################################################################################
################### Drop neg 60 and add AO #########################################################################################

fit.3aa<-gls(r6~AO.ed,correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)
AIC(fit.3aa)
AICc(fit.3aa)
aic.1a<-AIC(fit.3aa)
############################################################################################################################
################### Drop minus 60 degree and add lag 1 FHS #################################################################

fit.3b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3b)
D<-resid(fit.3b)
plot(D)
AIC(fit.3b)
AICc(fit.3b)
aic.1b<-AIC(fit.3b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.3bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3bb)
D<-resid(fit.3bb)
plot(D)
AIC(fit.3bb)
AICc(fit.3bb)
aic.1b<-AIC(fit.3bb)
############################################################################################################################
################### Add minus 60 degree component #########################################################################

fit.3c<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3c)
D<-resid(fit.3c)
plot(D)
AIC(fit.3c)
AICc(fit.3c)

aic.1c<-AIC(fit.3c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between minus 60 degree component and lag 2 FHS ####################################

fit.3d<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3d)
D<-resid(fit.3d)
plot(D)

AIC(fit.3d)
AICc(fit.3d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.3e<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3e)
D<-resid(fit.3e)
plot(D)

AIC(fit.3e)
AICc(fit.3e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.3ee<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3ee)
D<-resid(fit.3ee)
plot(D)

AIC(fit.3ee)
AICc(fit.3ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.3e3<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3e3)
D<-resid(fit.3e3)
plot(D)

AIC(fit.3e3)
AICc(fit.3e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 60 degree component ###########################

fit.3f<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.neg.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3f)
D<-resid(fit.3f)
plot(D)

AIC(fit.3f)
AICc(fit.3f)
#####################################################################################################################################################
##################### Add interaction term between arctic oscillation and minus 60 degree component AND swap in FHS lag 1 ###########################

fit.3ff<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.neg.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3ff)
D<-resid(fit.3ff)
plot(D)

AIC(fit.3ff)
AICc(fit.3ff)
########################################################################################################################################################
##################### Add Summer PDO ####################################################################################################################

fit.3g<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3g)
D<-resid(fit.3g)
plot(D)

AIC(fit.3g)
AICc(fit.3g)

lm.2g<-glm(r~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)
################################################################################################################################################
##################### Add uinteraction term between Summer PDO and minus 60 degree component ########################################################

fit.3h<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3h)
D<-resid(fit.3h)
plot(D)

AIC(fit.3h)
AICc(fit.3h)

###############################################################################################################################
#################### Reinstate interaction term between minus 60 and AO #######################################################
fit.3i<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.neg.60.May.June.ed
            +AO.ed*D.neg.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3i)
D<-resid(fit.3i)
plot(D)
AIC(fit.3i)
AICc(fit.3i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.3j<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3j)
D<-resid(fit.3j)
plot(D)
AIC(fit.3j)
AICc(fit.3j)

################################################################################################################################
##################### AS EXPERIMENT: Drop pdos ##############################################################################

fit.3jj<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3jj)
D<-resid(fit.3jj)
plot(D)
AIC(fit.3jj)
AICc(fit.3jj)
################################################################################################################################
#####################  Drop all interaction terms and swap FHS lag 1 in ##############################################################################

fit.3jjj<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3jjj)
D<-resid(fit.3jjj)
plot(D)
AIC(fit.3jjj)
AICc(fit.3jjj)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.3k<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3k)
D<-resid(fit.3k)
plot(D)

AIC(fit.3k)
AICc(fit.3k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.3m<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3m)
D<-resid(fit.3m)
plot(D)

AIC(fit.3m)
AICc(fit.3m)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 and swap in lag 1 FHS #########################################################################################

fit.3mm<-gls(r6~D.neg.60.May.June.ed+I(FHS.lag.1.ed^2)+May.June.ed+AO.ed+pdos.ed+D.neg.60.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3mm)
D<-resid(fit.3mm)
plot(D)

AIC(fit.3mm)
AICc(fit.3mm)

###########################################################################################################################################
##################### Exchange POSITIVE 90 for neg 60 ###########################################################################################
#############################################################################################################################################

#############################################################################################################################
###################### POSITIVE 90 degree component and lag 1 FHS with interaction term  #######################################

fit.4aa<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+D.90.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.4aa)
D<-resid(fit.4aa)
plot(D)

AIC(fit.4aa)
AICc(fit.4aa)
aic.1<-AIC(fit.4aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.4a<-gls(r6~D.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.4a)
D<-resid(fit.4a)
plot(D)
AIC(fit.4a)
AICc(fit.4a)
aic.1a<-AIC(fit.4a)

############################################################################################################################
################### Drop 90 degree and add lag 1 FHS #################################################################

fit.4b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.4b)
D<-resid(fit.4b)
plot(D)
AIC(fit.4b)
AICc(fit.4b)
aic.1b<-AIC(fit.3b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.4bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.4bb)
D<-resid(fit.4bb)
plot(D)
AIC(fit.4bb)
AICc(fit.4bb)
aic.1b<-AIC(fit.4bb)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.4c<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.4c)
D<-resid(fit.4c)
plot(D)
AIC(fit.4c)
AICc(fit.4c)

aic.1c<-AIC(fit.4c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic
############################################################################################################################
################### Add minus 90 degree component and swap in lag 1 FHS #########################################################################

fit.4cc<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.4cc)
D<-resid(fit.4cc)
plot(D)
AIC(fit.4cc)
AICc(fit.4cc)

aic.1c<-AIC(fit.4c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic
###################################################################################################################################
##################### Add interaction term between 90 degree component and lag 2 FHS ##############################################

fit.4d<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.4d)
D<-resid(fit.4d)
plot(D)

AIC(fit.4d)
AICc(fit.4d)

###################################################################################################################################
##################### Add arctic oscillation  #####################################################################################

fit.4e<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.4e)
D<-resid(fit.4e)
plot(D)

AIC(fit.4e)
AICc(fit.4e)

#####################################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.4ee<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+D.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.4ee)
D<-resid(fit.4ee)
plot(D)

AIC(fit.4ee)
AICc(fit.4ee)

###############################################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  ##############################################################################
#####################       NOTE: BAD model         ###########################################################################################

fit.4e3<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+D.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.4e3)
D<-resid(fit.4e3)
plot(D)

AIC(fit.4e3)
AICc(fit.4e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component ###########################

fit.4f<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.4f)
D<-resid(fit.4f)
plot(D)

AIC(fit.4f)
AICc(fit.4f)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and 90 degree component AND swap in FHS lag 1 ###########################

fit.4ff<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+D.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.4ff)
D<-resid(fit.4ff)
plot(D)

AIC(fit.4ff)
AICc(fit.4ff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.4g<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.4g)
D<-resid(fit.4g)
plot(D)

AIC(fit.4g)
AICc(fit.4g)

lm.2g<-glm(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)
##################################################################################################################################
##################### Add Summer PDO and swap in lag 1 FHS #############################################################################################

fit.4gg<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+D.90.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.4gg)
D<-resid(fit.4gg)
plot(D)

AIC(fit.4gg)
AICc(fit.4gg)

lm.2g<-glm(r6~D.neg.60.May.June.ed+I(FHS.lag.2.ed^2)+D.neg.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)
###############################################################################################################################
##################### Add uinteraction term between Summer PDO and 90 degree component ########################################################

fit.4h<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.4h)
D<-resid(fit.4h)
plot(D)

AIC(fit.4h)
AICc(fit.4h)

###############################################################################################################################
#################### Reinstate interaction term between 90 and AO #######################################################
fit.4i<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+D.90.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.90.May.June.ed+AO.ed*D.90.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.4i)
D<-resid(fit.4i)
plot(D)
AIC(fit.4i)
AICc(fit.4i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.4j<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.4j)
D<-resid(fit.4j)
plot(D)
AIC(fit.4j)
AICc(fit.4j)

################################################################################################################################
#####################  Drop all interaction terms and swap in FHS lag 1##############################################################################

fit.4jj<-gls(r6~D.90.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.4jj)
D<-resid(fit.4jj)
plot(D)
AIC(fit.4jj)
AICc(fit.4jj)
###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.4k<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.4k)
D<-resid(fit.4k)
plot(D)

AIC(fit.4k)
AICc(fit.4k)

################################################################################################################################
##################### Add interaction term for FHS and May-June 90 #########################################################################################

fit.4m<-gls(r6~D.90.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.90.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.4m)
D<-resid(fit.4m)
plot(D)

AIC(fit.4m)
AICc(fit.4m)




############################################################################################################################
###################### Drop 90 degrees and use POSITIVE 75 degrees ###################################################
##############################################################################################################################

fit.5aa<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+D.75.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.5aa)
D<-resid(fit.5aa)
plot(D)

AIC(fit.5aa)
AICc(fit.5aa)
aic.5<-AIC(fit.5aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.5a<-gls(r6~D.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.5a)
D<-resid(fit.5a)
plot(D)
AIC(fit.5a)
AICc(fit.5a)
aic.5a<-AIC(fit.5a)

############################################################################################################################
################### Drop minus 90 degree and add lag 1 FHS #################################################################

fit.5b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.5b)
D<-resid(fit.5b)
plot(D)
AIC(fit.5b)
AICc(fit.5b)
aic.5b<-AIC(fit.5b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.5bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.5bb)
D<-resid(fit.5bb)
plot(D)
AIC(fit.5bb)
AICc(fit.5bb)
aic.5bb<-AIC(fit.5bb)
############################################################################################################################
################### Add minus 90 degree component #########################################################################

fit.5c<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.5c)
D<-resid(fit.5c)
plot(D)
AIC(fit.5c)
AICc(fit.5c)
aic.1c<-AIC(fit.5c)


model.5.aic<-c(aic.5,aic.5a,aic.5b,aic.5c)
model.5.aic

##################################################################################################################################
##################### Add interaction term between minus 90 degree component and lag 2 FHS ####################################

fit.5d<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.5d)
D<-resid(fit.5d)
plot(D)

AIC(fit.5d)
AICc(fit.5d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.5e<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.5e)
D<-resid(fit.5e)
plot(D)

AIC(fit.5e)
AICc(fit.5e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.5ee<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+D.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.5ee)
D<-resid(fit.5ee)
plot(D)

AIC(fit.5ee)
AICc(fit.5ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.5e3<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+D.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.5e3)
D<-resid(fit.5e3)
plot(D)

AIC(fit.5e3)
AICc(fit.5e3)


###############################################################################################################################
##################### Add interaction term between arctic oscillation and 75 degree component ###########################

fit.5f<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.5f)
D<-resid(fit.5f)
plot(D)

AIC(fit.5f)
AICc(fit.5f)

###############################################################################################################################
##################### Add interaction term between arctic oscillation and minus 90 degree component AND swap in FHS lag 1 ###########################

fit.5ff<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+D.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.5ff)
D<-resid(fit.5ff)
plot(D)

AIC(fit.5ff)
AICc(fit.5ff)

###############################################################################################################################
##################### AS EXPERIMENT: remove all interaction terms and AND swap in FHS lag 1 ###########################

fit.5fff<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.5fff)
D<-resid(fit.5fff)
plot(D)

AIC(fit.5fff)
AICc(fit.5fff)
##################################################################################################################################
##################### Add Summer PDO #############################################################################################

fit.5g<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.5g)
D<-resid(fit.5g)
plot(D)

AIC(fit.5g)
AICc(fit.5g)

lm.1g<-glm(r6~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)

##################################################################################################################################
##################### Add Summer PDO and swap in lag 1 FHS #############################################################################################

fit.5gg<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+D.75.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.5gg)
D<-resid(fit.5gg)
plot(D)

AIC(fit.5gg)
AICc(fit.5gg)

lm.1g<-glm(r6~D.75.May.June+I(FHS.lag.2^2)+D.75.May.June*I(FHS.lag.2^2)+AO+pdos)
step(lm.1g)

###############################################################################################################################
##################### Add interaction term between Summer PDO and minus 90 degree component ########################################################

fit.5h<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.5h)
D<-resid(fit.5h)
plot(D)

AIC(fit.5h)
AICc(fit.5h)

###############################################################################################################################
#################### Reinstate interaction term between minus 90 and AO #######################################################
fit.5i<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+D.75.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.75.May.June.ed+AO.ed*D.75.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.5i)
D<-resid(fit.5i)
plot(D)
AIC(fit.5i)
AICc(fit.5i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.5j<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.5j)
D<-resid(fit.5j)
plot(D)
AIC(fit.5j)
AICc(fit.5j)

################################################################################################################################
#####################  Drop all interaction terms and swap in lag 1 FHS ##############################################################################

fit.5j<-gls(r6~D.75.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.5j)
D<-resid(fit.5j)
plot(D)
AIC(fit.5j)
AICc(fit.5j)

###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.5k<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.5k)
D<-resid(fit.5k)
plot(D)

AIC(fit.5k)
AICc(fit.5k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.5m<-gls(r6~D.75.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.75.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.5m)
D<-resid(fit.5m)
plot(D)

AIC(fit.5m)
AICc(fit.5m)







###########################################################################################################################################
##################### Exchange 60 degree component for 75 ###########################################################################################
#############################################################################################################################################






#############################################################################################################################
###################### 60 degree component and lag 1 FHS with interaction term  #######################################

fit.3aa<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+D.60.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)

AIC(fit.3aa)
AICc(fit.3aa)
aic.1<-AIC(fit.3aa)



############################################################################################################################
################### Drop lag 1 FHS #########################################################################################

fit.3a<-gls(r6~D.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3a)
D<-resid(fit.3a)
plot(D)
AIC(fit.3a)
AICc(fit.3a)
aic.1a<-AIC(fit.3a)

############################################################################################################################
################### Drop neg 60 and add AO #########################################################################################

fit.3aa<-gls(r6~AO.ed,correlation=corAR1(),method="ML")
summary(fit.3aa)
D<-resid(fit.3aa)
plot(D)
AIC(fit.3aa)
AICc(fit.3aa)
aic.1a<-AIC(fit.3aa)
############################################################################################################################
################### Drop minus 60 degree and add lag 1 FHS #################################################################

fit.3b<-gls(r6~I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3b)
D<-resid(fit.3b)
plot(D)
AIC(fit.3b)
AICc(fit.3b)
aic.1b<-AIC(fit.3b)
############################################################################################################################
################### Drop lag 1 FHS and add lag 2 FHS #################################################################

fit.3bb<-gls(r6~I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3bb)
D<-resid(fit.3bb)
plot(D)
AIC(fit.3bb)
AICc(fit.3bb)
aic.1b<-AIC(fit.3bb)
############################################################################################################################
################### Add 60 degree component #########################################################################

fit.3c<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3c)
D<-resid(fit.3c)
plot(D)
AIC(fit.3c)
AICc(fit.3c)

aic.1c<-AIC(fit.3c)


model.1.aic<-c(aic.1,aic.1a,aic.1b,aic.1c)
model.1.aic

##################################################################################################################################
##################### Add interaction term between 60 degree component and lag 2 FHS ####################################

fit.3d<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3d)
D<-resid(fit.3d)
plot(D)

AIC(fit.3d)
AICc(fit.3d)

###############################################################################################################################
##################### Add arctic oscillation  #################################################################################

fit.3e<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3e)
D<-resid(fit.3e)
plot(D)

AIC(fit.3e)
AICc(fit.3e)

###############################################################################################################################
##################### swap lag 1 FHS for lag 2 FHS  #################################################################################

fit.3ee<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+D.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3ee)
D<-resid(fit.3ee)
plot(D)

AIC(fit.3ee)
AICc(fit.3ee)

###############################################################################################################################
##################### Add lag 2 FHS with interaction with lag 1  #################################################################################
#####################       NOTE: BAD model         ###########################################################################

fit.3e3<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+D.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+I(FHS.lag.2.ed^2)*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3e3)
D<-resid(fit.3e3)
plot(D)

AIC(fit.3e3)
AICc(fit.3e3)
###############################################################################################################################
##################### Add interaction term between arctic oscillation and 60 degree component ###########################

fit.3f<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+AO.ed*D.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3f)
D<-resid(fit.3f)
plot(D)

AIC(fit.3f)
AICc(fit.3f)
#####################################################################################################################################################
##################### Add interaction term between arctic oscillation and 60 degree component AND swap in FHS lag 1 ###########################

fit.3ff<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+D.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+AO.ed*D.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3ff)
D<-resid(fit.3ff)
plot(D)

AIC(fit.3ff)
AICc(fit.3ff)
########################################################################################################################################################
##################### Add Summer PDO ####################################################################################################################

fit.3g<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3g)
D<-resid(fit.3g)
plot(D)

AIC(fit.3g)
AICc(fit.3g)

lm.2g<-glm(r~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)

########################################################################################################################################################
##################### Add Summer PDO and swap in lag 2 FHS ####################################################################################################################

fit.3gg<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+D.60.May.June.ed*I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3gg)
D<-resid(fit.3gg)
plot(D)

AIC(fit.3gg)
AICc(fit.3gg)

lm.2g<-glm(r~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed)
step(lm.2g)
################################################################################################################################################
##################### Add uinteraction term between Summer PDO and 60 degree component ########################################################

fit.3h<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3h)
D<-resid(fit.3h)
plot(D)

AIC(fit.3h)
AICc(fit.3h)

###############################################################################################################################
#################### Reinstate interaction term between 60 and AO #######################################################
fit.3i<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+D.60.May.June.ed*I(FHS.lag.2.ed^2)+AO.ed+pdos.ed+pdos.ed*D.60.May.June.ed
            +AO.ed*D.60.May.June.ed,correlation=corAR1(),method="ML")
summary(fit.3i)
D<-resid(fit.3i)
plot(D)
AIC(fit.3i)
AICc(fit.3i)

################################################################################################################################
#####################  Drop all interaction terms ##############################################################################

fit.3j<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3j)
D<-resid(fit.3j)
plot(D)
AIC(fit.3j)
AICc(fit.3j)

################################################################################################################################
##################### AS EXPERIMENT: Drop pdos ##############################################################################

fit.3jj<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3jj)
D<-resid(fit.3jj)
plot(D)
AIC(fit.3jj)
AICc(fit.3jj)
################################################################################################################################
#####################  Drop all interaction terms and swap FHS lag 1 in ##############################################################################

fit.3jjj<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3jjj)
D<-resid(fit.3jjj)
plot(D)
AIC(fit.3jjj)
AICc(fit.3jjj)

################################################################################################################################
#####################  Drop all interaction terms and 60 degree and swap FHS lag 1 in ##############################################################################

fit.3jjjj<-gls(r6~I(FHS.lag.1.ed^2)+AO.ed,correlation=corAR1(),method="ML")
summary(fit.3jjjj)
D<-resid(fit.3jjjj)
plot(D)
AIC(fit.3jjjj)
AICc(fit.3jjjj)
###############################################################################################################################
##################### Add May-June SST #########################################################################################

fit.3k<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed,correlation=corAR1(),method="ML")
summary(fit.3k)
D<-resid(fit.3k)
plot(D)

AIC(fit.3k)
AICc(fit.3k)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 #########################################################################################

fit.3m<-gls(r6~D.60.May.June.ed+I(FHS.lag.2.ed^2)+May.June.ed+AO.ed+pdos.ed+D.60.May.June.ed*I(FHS.lag.2.ed^2),correlation=corAR1(),method="ML")
summary(fit.3m)
D<-resid(fit.3m)
plot(D)

AIC(fit.3m)
AICc(fit.3m)

################################################################################################################################
##################### Add interaction term for FHS and May-June -90 and swap in lag 1 FHS #########################################################################################

fit.3mm<-gls(r6~D.60.May.June.ed+I(FHS.lag.1.ed^2)+May.June.ed+AO.ed+pdos.ed+D.60.May.June.ed*I(FHS.lag.1.ed^2),correlation=corAR1(),method="ML")
summary(fit.3mm)
D<-resid(fit.3mm)
plot(D)

AIC(fit.3mm)
AICc(fit.3mm)



























#################################################################################################################################
fit.4<-gls(r~D.neg.90.May.June+I(FHS.lag.2^2)+May.June.anom+D.neg.75.May.June+D.neg.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.4)
D<-resid(fit.4)
plot(D)

############################################################################################################################
##################### target D.neg.75 in interaction term + drop SST terms + add lag 1 FHS effect ##########################

fit.5<-gls(r~I(FHS.lag.2^2)+I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.5)
D<-resid(fit.5)
plot(D)

################################################################################################################################
##################### add lag 1 and 2 FHS interaction term #####################################################################

fit.6<-gls(r~I(FHS.lag.2^2)*I(FHS.lag.1^2)+I(FHS.lag.2^2)+I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.6)
D<-resid(fit.6)
plot(D)

############################################################################################################################
###################### remove lag 2 FHS ###################################################################################

fit.7<-gls(r~I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.75.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.7)
D<-resid(fit.7)
plot(D)
aic.7<-AIC(fit.7)
aic.7


###############################################################################################################################
##################### Drop minus 75 degree component-FHS lag 1 interaction, add minus 60 degree component #####################


fit.7<-gls(r~I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.60.May.June,correlation=corAR1(),method="ML")
summary(fit.7)
D<-resid(fit.7)
plot(D)


###############################################################################################################################
##################### drop minus 75 and 60 degree components and add May-June SST anomaly and male 3 year rolling average #######
fit.8<-gls(r~I(FHS.lag.2^2)+May.June.anom+RA.3.males,correlation=corAR1())
summary(fit.8)




###############################################################################################################################
################### Drop minus 90 degree component and add minus 75 degree component #########################################
fit.10<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.10)
D<-resid(fit.10)
plot(D)

###############################################################################################################################
################### Add interaction term between minus 75 degree component and lag 2 FHS ######################################
fit.11<-gls(r~D.neg.75.May.June+I(FHS.lag.2^2)+D.neg.75.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.11)
D<-resid(fit.11)
plot(D)

###############################################################################################################################
################### drop neg 75 and interaction term and add minus 60 degree component #######################################

fit.12<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.12)
D<-resid(fit.12)
plot(D)

###############################################################################################################################
################### Add interaction term between minus 60 degree component and lag 2 FHS ######################################
fit.13<-gls(r~D.neg.60.May.June+I(FHS.lag.2^2)+D.neg.60.May.June*I(FHS.lag.2^2),correlation=corAR1(),method="ML")
summary(fit.13)
D<-resid(fit.13)
plot(D)


###############################################################################################################################
################## Substitute lag 1 FHS for lag 2 FHS #########################################################################
fit.16<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^1),correlation=corAR1(),method="ML")
summary(fit.16)
D<-resid(fit.16)
plot(D)



###############################################################################################################################
fit.19<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.19)
D<-resid(fit.19)
plot(D)


###############################################################################################################################
fit.20<-gls(r~D.neg.60.May.June+I(FHS.lag.1^2)+D.neg.60.May.June*I(FHS.lag.1^1),correlation=corAR1(),method="ML")
summary(fit.20)
D<-resid(fit.20)
plot(D)

#################################################################################################################################
####################### 90 degree component and lag 1 FHS  ######################################################################

fit.21<-gls(r~D.90.May.June+I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.21)
D<-resid(fit.21)
plot(D)

###############################################################################################################################
################## Add interaction term between 90 degree and lag 1 FHS #################################################
fit.22<-gls(r~D.90.May.June+I(FHS.lag.1^2)+D.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.22)
D<-resid(fit.22)
plot(D)


###############################################################################################################################
################## Add lag 2 FHS with interaction terms with 90 degree component and lag 1 FHS #####################################
fit.23<-gls(r~I(FHS.lag.2^2)*I(FHS.lag.1^2)+I(FHS.lag.2^2)+I(FHS.lag.1^2)+D.90.May.June+D.90.May.June*I(FHS.lag.2^2)+D.90.May.June*I(FHS.lag.1^2),correlation=corAR1(),method="ML")
summary(fit.23)
D<-resid(fit.23)
plot(D)













##############################################################################################################################
#################################### GLM models ##############################################################################
fit.1<-glm(r~D.neg.90.May.June+I(FHS.lag.2^2)+May.June.anom+D.neg.75.May.June)
summary(fit.1)
step(fit.1)
D<-resid(fit.1)
plot(D)

fit.2<-glm(r~I(FHS.lag.2^2)+May.June.anom+D.neg.75.May.June+)
summary(fit.2)
step(fit.2)
D<-resid(fit.2)
plot(D)

fit.3<-glm(r~I(FHS.lag.2^2)+D.neg.75.May.June+D.neg.60.May.June)
summary(fit.3)
step(fit.3)
D<-resid(fit.3)
plot(D)

fit.4<-glm(r~I(FHS.lag.2^2)+I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.60.May.June)
summary(fit.4)
step(fit.4)
D<-resid(fit.4)
plot(D)

fit.5<-glm(r~I(FHS.lag.2^2)+I(FHS.lag.1^2)+D.neg.75.May.June+D.neg.60.May.June+RA.3.males)
summary(fit.5)
step(fit.5)
D<-resid(fit.5)
plot(D)


fit.6<-glm(r~D.neg.90.May.June+I(FHS.lag.1^2)+D.neg.90.May.June*I(FHS.lag.1^2))
summary(fit.6)
step(fit.6)
D<-resid(fit.6)
plot(D)
AIC(step(fit.6))

fit.8<-lm(r~I(FHS.lag.1^2)+I(FHS.lag.1^2)+AO+AO*I(FHS.lag.1^2))
summary(fit.8)
D<-resid(fit.8)
plot(D)

AIC(fit.8)
AICc(fit.8)
