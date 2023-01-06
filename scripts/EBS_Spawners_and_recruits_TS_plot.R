
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
sp_max_cs<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_THREE_QRT+Spawners$NUM_FEMALE_FULL))))
sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))

colnames(sp_ovig)<-c("Year", "EBS_ovigfem_abun")
colnames(sp_max_cs)<-c("Year", "EBS_max_cs_fem_abun")
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")

e166_sp_ovig<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_OVIGEROUS))))
e166_sp_max_cs<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_THREE_QRT+E166Spawners$NUM_FEMALE_FULL))))
e166_sp_sc3<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3))))

colnames(e166_sp_ovig)<-c("Year", "E166_ovigfem_abun")
colnames(e166_sp_max_cs)<-c("Year", "E166_max_cs_fem_abun")
colnames(e166_sp_sc3)<-c("Year", "E166_sc3fem_abun")

w166_sp_ovig<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_OVIGEROUS))))
w166_sp_max_cs<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_THREE_QRT+W166Spawners$NUM_FEMALE_FULL))))
w166_sp_sc3<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3))))

colnames(w166_sp_ovig)<-c("Year", "W166_ovigfem_abun")
colnames(w166_sp_max_cs)<-c("Year", "W166_max_cs_fem_abun")
colnames(w166_sp_sc3)<-c("Year", "W166_sc3fem_abun")

Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
Spawners_Max_CS<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_THREE_QRT+Spawners$NUM_FEMALE_FULL),(E166Spawners$NUM_FEMALE_THREE_QRT+E166Spawners$NUM_FEMALE_FULL),(W166Spawners$NUM_FEMALE_THREE_QRT+W166Spawners$NUM_FEMALE_FULL))))
Spawners_ovig<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_OVIGEROUS),(E166Spawners$NUM_FEMALE_OVIGEROUS),(W166Spawners$NUM_FEMALE_OVIGEROUS))))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3","E166_SC3","W166_SC3")
colnames(Spawners_Max_CS)<-c("Year", "EBS_Max_CS","E166_Max_CS","W166_Max_CS")
colnames(Spawners_ovig)<-c("Year", "EBS_ovig","E166_ovig","W166_ovig")

Spawners_sc4<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC4),(E166Spawners$NUM_FEMALE_SC4),(W166Spawners$NUM_FEMALE_SC4))))

Spawners_os<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4),(E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4),(W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4))))

colnames(Spawners_sc4)<-c("Year", "EBS_sc4fem_abun","E166_sc4fem_abun","W166_sc4fem_abun")
colnames(Spawners_os)<-c("Year", "EBS_sc3and4fem_abun", "E166_sc3and4fem_abun", "W166_sc3and4fem_abun")
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
###############################################################################################################################

###############################################################################################################################
########################################## Shell condition 3 #################################################################
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

############################################################################################################################################################
########################################## Three quarters full and full clutch size females #################################################################
colnames(Spawners_ovig)<-c("Year", "EBS_ovig","E166_ovig","W166_ovig")

Spawners_Max_CS$Year
S<-Spawners_Max_CS$EBS_Max_CS[4:45]
########################################## Full Timeseries #############################################################

acf(S,main="Max CS females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title
#?acf()

########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008

acf(S3,main="Max CS females 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S1,main="Max CS females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="Max CS females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)


############################################################################################################################################################
########################################## All ovigerous females #################################################################

Spawners_ovig$Year
S<-Spawners_ovig$EBS_ovig[4:45]
########################################## Full Timeseries #############################################################

acf(S,main="Ovigerous females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title
#?acf()

########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008

acf(S3,main="Ovigerous females 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S1,main="Ovigerous females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="Ovigerous females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
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
###################################### Plot for paper ###########################################
y.range<-range(0,4e+08)
dev.new()
par(mfrow=c(2,1),cex.main=1.2, cex.lab=1.2,cex.axis=1.2,cex=1.2) #configure axis labels

#########################################################################################################################
################################## Plot all juveniles together ##########################################################
EBS.juv.a<-rec_30to50$EBS_Abun_30to50
EBS.juv<-EBS.juv.a[4:length(EBS.juv.a)]
rec_30to50$Year
plot(EBS.juv,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)
#abline(h=0)

axis(1, at=1:42,lab=rec_30to50$Year[4:length(rec_30to50$Year)],tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
#legend(14,9e+08,c("Total EBS", "Western area", "Eastern area"),cex=1.25,col=c(1,1,1),pch=c(2,18,22))
#legend(15,9e+08,c("EBS juvenile estimate"),cex=1,col=c(4,1,2),pch=c(16,18,22))
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main = "EBS juvenile abundance by year")
#title(main=" Juvenile abundance by year")


############################## Plot Spawners ###########################################

EBS.fem<-Spawners_SC3$EBS_SC3[4:length(Spawners_SC3$EBS_SC3)]
EBS.fem
plot(EBS.fem,type="b",ylim=y.range, pch=2,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1:42,lab=Spawners_SC3$Year[4:length(Spawners_SC3$Year)],tck=0.03)
axis(2, las=1, at=1e+08*0:4,labels=c("0","100","200","300","400"),tck=0.03)
#abline(h=0)
box()

title(xlab="Year")
title(ylab="Abundance (Millions)")
title(main=" Oldshell female abundance by year")



