
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

Spawners_SC0<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC0,E166Spawners$NUM_FEMALE_SC0,W166Spawners$NUM_FEMALE_SC0)))
Spawners_SC1<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC1,E166Spawners$NUM_FEMALE_SC1,W166Spawners$NUM_FEMALE_SC1)))
Spawners_SC2<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2)))
Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
Spawners_SC4<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC4,E166Spawners$NUM_FEMALE_SC4,W166Spawners$NUM_FEMALE_SC4)))
Spawners_SC5<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC5,E166Spawners$NUM_FEMALE_SC5,W166Spawners$NUM_FEMALE_SC5)))
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

Spawners_analysis<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,
                                                 Spawners$NUM_FEMALE_SC0,E166Spawners$NUM_FEMALE_SC0,W166Spawners$NUM_FEMALE_SC0,
                                                 Spawners$NUM_FEMALE_SC1,E166Spawners$NUM_FEMALE_SC1,W166Spawners$NUM_FEMALE_SC1,
                                                 Spawners$NUM_FEMALE_SC2,E166Spawners$NUM_FEMALE_SC2,W166Spawners$NUM_FEMALE_SC2,
                                                 Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3,
                                                 Spawners$NUM_FEMALE_SC4,E166Spawners$NUM_FEMALE_SC4,W166Spawners$NUM_FEMALE_SC4,
                                                 Spawners$NUM_FEMALE_SC5,E166Spawners$NUM_FEMALE_SC5,W166Spawners$NUM_FEMALE_SC5,
                                                 (Spawners$NUM_FEMALE_SC0+Spawners$NUM_FEMALE_SC1+Spawners$NUM_FEMALE_SC2+Spawners$NUM_FEMALE_SC3+Spawners$NUM_FEMALE_SC4+Spawners$NUM_FEMALE_SC5),
                                                 (E166Spawners$NUM_FEMALE_SC0+E166Spawners$NUM_FEMALE_SC1+E166Spawners$NUM_FEMALE_SC2+E166Spawners$NUM_FEMALE_SC3+E166Spawners$NUM_FEMALE_SC4+E166Spawners$NUM_FEMALE_SC5),
                                                 (W166Spawners$NUM_FEMALE_SC0+W166Spawners$NUM_FEMALE_SC1+W166Spawners$NUM_FEMALE_SC2+W166Spawners$NUM_FEMALE_SC3+W166Spawners$NUM_FEMALE_SC4+W166Spawners$NUM_FEMALE_SC5),
                                                 Spawners$NUM_FEMALE_OVIGEROUS,E166Spawners$NUM_FEMALE_OVIGEROUS,W166Spawners$NUM_FEMALE_OVIGEROUS,
                                                 Spawners_OS$EBS_os,Spawners_OS$E166_os, Spawners_OS$W166_os)))

colnames(Spawners_analysis)<-c("Year","EBS_SC0","E166_SC0","W166_SC0",
                               "EBS_SC1","E166_SC1","W166_SC1",
                               "EBS_SC2","E166_SC2","W166_SC2",
                               "EBS_SC3","E166_SC3","W166_SC3",
                               "EBS_SC4","E166_SC4","W166_SC4",
                               "EBS_SC5","E166_SC5","W166_SC5",
                               "EBS_SC_TOTAL","E166_SC_TOTAL","W166_SC_TOTAL",
                               "EBS_ovig","E166_ovig","W166_ovig",
                               "EBS_os","E166_os","W166_os")

#Spawners_analysis
#write.csv(Spawners_analysis,"./data/Female_Tanner_Crab_Series_for_analysis.csv")