
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
y.range<-range(0,3.5e+08)

plot(EBS.juv)
plot(EBS.juv,type="b",ylim=y.range, pch=16,col=4,axes=FALSE,ann=FALSE)
lines(W166.juv,type="b",pch=18,col=1)
lines(E166.juv,type="b",pch=22,col=2)
#abline(h=0)

axis(1, at=1:n,lab=rec_30to50$Year,tck=0.02)
axis(2, las=1, at=1e+08*0:9,labels=c("0","100","200","300","400","500","600","700","800","900"),tck=0.02)
legend(14,3.5e+08,c("Total EBS", "Western area", "Eastern area"),cex=1,col=c(4,1,2),pch=c(16,18,22))
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


####################################Fit gls models with 1st and 2nd order autocorrelation######################
#par(mfrow=c(1,1))
Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main= "Lag 2 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)

r1<-resid(Fit3)
#plot(1r,pch=16)
r.fit3<-loess(r1~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 2 S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)


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

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
#par(mfrow=c(1,1))

par(mfrow=c(2,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)



r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()
dev.new()
par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
r6<-resid(Fit3)
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
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="SC3 female abundance (Millions)")
title(ylab="EBS log-survival")
title(main="Lag-3 yrs")


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

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

Fit5<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit5)

plot(log(R/S)~S, xlab="Reproductive female abundance",ylab="log(R/S)",main="Lag 4 S-R model",pch=16)
abline(Fit5)
#plot(Fit3)

coef(Fit5)

par(mfrow=c(1,1))
r5<-resid(Fit3)
#plot(r5,pch=16)
r.fit5<-loess(r5~xyear.k)
plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 4 S/R residuals")
lines(xyear.k,fitted(r.fit5),col=4)



##############################################################################################################################################
########################################### SHELL CONDITION 3 AND 4 COMBINED FEMALES  ########################################################
###############################################################################################################################################

########################################### Lag=2 ###########################################################################################
names(Spawners_os)
########################################### Recreate data series for actual analysis ########################################################
par(mfrow=c(3,1))
n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_os$EBS_sc3and4fem_abun[4:n]

recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_os$Year[4:n]
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


####################################Fit gls models with 1st and 2nd order autocorrelation######################
#par(mfrow=c(1,1))
Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)

plot(log(R/S)~S, xlab="SC3 and 4 reproductive female abundance",ylab="log(R/S)",main= "Lag 2 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)

r1<-resid(Fit3)
#plot(1r,pch=16)
r.fit3<-loess(r1~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 2 S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)


########################################################################################################
######################################### Lag=3 ########################################################
#par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_os$EBS_sc3and4fem_abun[4:n]

R
S
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_os$Year[4:n]
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


#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="SC3 and 4 reproductive female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()

par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels

r6<-resid(Fit3)
xyear<-xyear.k
y.range<-range(0,9e+08)

###raw estimates

#plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)
#axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
#axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
#box()

#title(xlab="EBS SC3 and SC4 female abundance (Millions)")
#title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")

### log-survival
plot(log(R/S)~S,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="EBS SC3 and SC4 female abundance (Millions)")
title(ylab="EBS log-survival")
title(main="SC3 and SC4 females")


### residual trends
##r.fit3<-loess(r6~xyear,span=0.30)
##plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
##lines(xyear,fitted(r.fit3),col=1)

r.fit3<-loess(r6~xyear,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()
title(xlab="Hatch year")
title(ylab="EBS S-R residuals")

##############################################################################################################################
################################################# Lag=4 ######################################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_os$EBS_sc3and4fem_abun[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_os$Year[4:n]

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

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

Fit5<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit5)

plot(log(R/S)~S, xlab="EBS SC3 and 4 reproductive female abundance",ylab="log(R/S)",main="Lag 4 S-R model",pch=16)
abline(Fit5)
#plot(Fit3)

coef(Fit5)

par(mfrow=c(1,1))
r5<-resid(Fit3)
#plot(r5,pch=16)
r.fit5<-loess(r5~xyear.k)
plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 4 S/R residuals")
lines(xyear.k,fitted(r.fit5),col=4)

##############################################################################################################################################
########################################### Three quarters full and full clutch size females  ################################################
##############################################################################################################################################

########################################### Recreate data series for actual analysis ########################################################
par(mfrow=c(3,1))

########################################################################################################
######################################### Lag=3 - equivalent to lag 2 for SC model########################################################
#par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_Max_CS$EBS_Max_CS[4:n]

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


#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Large clutch size female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()

par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
box()

title(xlab="EBS Large clutch size female abundance (Millions)")
title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")


plot(log(R/S)~S,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Large clutch size Female abundance (Millions)")
title(ylab="EBS log-survival")

r6<-resid(Fit3)
xyear<-xyear.k
#plot(r6,pch=16)
r.fit3<-loess(r6~xyear,span=0.30)
plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
lines(xyear,fitted(r.fit3),col=1)

r.fit3<-loess(r6~xyear,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")

##############################################################################################################################
################################################# Lag=4 - equivalent to lag 3 for SC models######################################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_Max_CS$EBS_Max_CS[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_Max_CS$Year[4:n]

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

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

#################################Fit gls model with 1st order autocorrelation#######################################################

Fit7<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit7)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Large clutch size female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit7)

#plot(Fit3)
coef(Fit7)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r7<-resid(Fit7)
#plot(r3,pch=16)
r.fit7<-loess(r7~xyear.k)
plot(r.fit7,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit7),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()
dev.new()
par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
r7<-resid(Fit7)
xyear<-xyear.k


#### Raw abundance
#plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)
#axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
#axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
#box()
#title(xlab="EBS Large clutch size female abundance (Millions)")
#title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")

#### log-survival
plot(log(R/S)~S,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit7)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Large clutch size Female abundance (Millions)")
title(ylab="EBS log-survival")
title(main="Large clutch size Females")


#### residual trends
#plot(r6,pch=16)
#r.fit3<-loess(r6~xyear,span=0.30)
#plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
#lines(xyear,fitted(r.fit3),col=1)

r.fit7<-loess(r7~xyear,span=0.5)
plot(r.fit7,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit7),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")


##############################################################################################################################
################################################# Lag=5 - equivalent to lag 4 for SC models######################################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_Max_CS$EBS_Max_CS[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_Max_CS$Year[4:n]

xi<-S
xyear<-spYear[-5]
yi<-R
yyear<-recYear[-9]
lag <- 5
n <- length(yi)			##########note change from xi as previously employed
xi.k <- xi[1:(n-lag)]       # Select reproductive female ests
S<-xi.k
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k 
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

#################################Fit gls model with 1st order autocorrelation#######################################################

Fit7<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit7)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Large clutch size female abundance",ylab="log(R/S)",main="Lag 6 S-R model",pch=16)
abline(Fit7)

#plot(Fit3)
coef(Fit7)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r7<-resid(Fit7)
#plot(r3,pch=16)
r.fit7<-loess(r7~xyear.k)
plot(r.fit7,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit7),col=4)
##############################################################################################################################################
########################################### Ovigerous females  ################################################
##############################################################################################################################################
colnames(Spawners_ovig)

########################################### Recreate data series for actual analysis ########################################################
par(mfrow=c(3,1))


########################################################################################################
######################################### Lag = 3, equivalent to lag = 2 for SC models ########################################################
#par(mfrow=c(1,1),cex.main=1.25, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_ovig$EBS_ovig[4:n]

R
S
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_ovig$Year[4:n]
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


#################################Fit gls model with 1st order autocorrelation#######################################################

Fit3<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit3)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Ovigerous female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit3)

#plot(Fit3)
coef(Fit3)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r3<-resid(Fit3)
#plot(r3,pch=16)
r.fit3<-loess(r3~xyear.k)
plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit3),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()

par(mfrow=c(1,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
plot(R2~S2,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)


axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
box()

title(xlab="EBS ovigerous female abundance (Millions)")
title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")


plot(log(R2/S2)~S2,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit3)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Ovigerous female abundance (Millions)")
title(ylab="EBS log-survival")

r6<-resid(Fit3)
xyear<-xyear.k
#plot(r6,pch=16)
r.fit3<-loess(r6~xyear,span=0.30)
plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
lines(xyear,fitted(r.fit3),col=1)

r.fit3<-loess(r6~xyear,span=0.5)
plot(r.fit3,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit3),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")

##############################################################################################################################
################################################# Lag=4--equivalent to lag 3 in SC models ######################################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_ovig$EBS_ovig[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_ovig$Year[4:n]

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

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

Fit8<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit8)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Large clutch size female abundance",ylab="log(R/S)",main="Lag 3 S-R model",pch=16)
abline(Fit8)

#plot(Fit3)
coef(Fit8)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r8<-resid(Fit8)
#plot(r3,pch=16)
r.fit8<-loess(r8~xyear.k)
plot(r.fit8,pch=16,xlab="Hatch year", ylab="Lag 3 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit8),col=4)



##################################################################################################################################
################################################## Plot above for MEPS #########################################################
#,mai=c(1,1.25,1,1)#add to below to format margins if necessary--in inches, for in lines, use mar()
dev.new()
par(mfrow=c(2,1),cex.lab=1.55,cex.axis=1.25,cex=1.25) #configure axis labels
y.range<-range(0,9e+08)
r8<-resid(Fit8)
xyear<-xyear.k


#### Raw abundance
#plot(R~S,ylim=y.range, pch=16,col=1,axes=FALSE,ann=FALSE)
#axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
#axis(2, las=1, at=2e+08*c(0,1,2,3,4,5),labels=c("0","200","400","600","800","1000"),tck=0.02)
#box()
#title(xlab="EBS Large clutch size female abundance (Millions)")
#title(ylab="EBS juvenile abundance (Millions)")
#title(main=" Juvenile abundance by year")

#### log-survival
plot(log(R/S)~S,pch=16,col=1,axes=FALSE,ann=FALSE)
abline(Fit8)
axis(1, at=1e+07*c(0,5,10,15,20),labels=c("0","50","100","150","200"),tck=0.02)
axis(2, las=1, at=c(-1,0,1,2,3,4,5),labels=c("-1","0","1","2","3","4","5"),tck=0.02)
box()
title(xlab="Large clutch size Female abundance (Millions)")
title(ylab="EBS log-survival")
title(main="Ovigerous Females")


#### residual trends
#plot(r6,pch=16)
#r.fit3<-loess(r6~xyear,span=0.30)
#plot(r.fit3,pch=16,xlab="Hatch year", ylab=" S-R residuals")
#lines(xyear,fitted(r.fit3),col=1)

r.fit8<-loess(r8~xyear,span=0.5)
plot(r.fit8,pch=16,col=1,axes=FALSE,ann=FALSE)
lines(xyear,fitted(r.fit8),col=1)
axis(1, at=c(1980,1985,1990,1995,2000,2005,2010,2015),labels=c("1980","1985","1990","1995","2000","2005","2010","2015"),tck=0.02)
axis(2, las=1,at=c(-3,-2,-1,0,1,2,3),labels=c("-3","-2","-1","0","1","2","3"),tck=0.02)
box()

title(xlab="Hatch year")
title(ylab="EBS S-R residuals")



##############################################################################################################################
################################################# Lag=5--equivalent to lag 4 in SC models ######################################################################

n<-nrow(rec_30to50)
R<-rec_30to50$EBS_Abun_30to50[4:n]
S<-Spawners_ovig$EBS_ovig[4:n]

R<-R[-9]
S<-S[-5]
recYear<-rec_30to50$Year[4:n]
spYear<-Spawners_ovig$Year[4:n]

xi<-S
xyear<-spYear[-5]
yi<-R
yyear<-recYear[-9]
lag <- 5
n <- length(yi)			##########note change from xi as previously employed
xi.k <- xi[1:(n-lag)]       # Select reproductive female ests
S<-xi.k
xyear.k<-xyear[1:(n-lag)]   # Select corresponding years
xyear.k
yi.k <- yi[(lag+1):n]       # Select Juvenile recruitment estimates'lag' years later
R<-yi.k 
yyear.k<-yyear[(lag+1):n]   # Select years corresponding to juvenile recruitment estimates
yyear.k

####################################################################################################################################
############################################## Fit gls model with 1st order autocorrelation ########################################

Fit8<-gls(log(R/S)~S,correlation=corAR1(), na.action=na.omit)
summary(Fit8)
#par(mfrow=c(1,1))

plot(log(R/S)~S, xlab="Large clutch size female abundance",ylab="log(R/S)",main="Lag 5 S-R model",pch=16)
abline(Fit8)

#plot(Fit3)
coef(Fit8)


par(mfrow=c(1,1),cex.main=1.5, cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels
r8<-resid(Fit8)
#plot(r3,pch=16)
r.fit8<-loess(r8~xyear.k)
plot(r.fit8,pch=16,xlab="Hatch year", ylab="Lag 5 S/R residuals",main="Interdecadal trend in EBS group S/R residuals")
lines(xyear.k,fitted(r.fit8),col=4)