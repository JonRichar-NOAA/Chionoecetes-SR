###########################################################################################################

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

opies<-read.csv("data/co_male_female_group_abun_ebs.csv")
opies

opies_edit<-read.csv("data/co_male_female_group_abun_eb_noCI_CVs.csv")
opies_edit

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

colnames(Spawners_SC3)<-c("Year", "EBS_CBfem_SC3")

########################################################################################################
juvs<-cbind(rec_30to40[,2],rec_40to50[,2],rec_50to60[,2],rec_30to50[,2],rec_30to60[,2],rec_40to60[,2])
colnames(juvs)<-c("juvs30to40","juvs40to50","juvs50to60","juvs30to50","juvs30to60","juvs40to60")
sub.juvs<-juvs[6:nrow(juvs),]
sub<-Spawners_SC3$EBS_CBfem_SC3[6:nrow(Spawners_SC3)]

sub<-as.matrix(sub)
sub
colnames(sub)<-"EBS_CBfem_SC3"
X<-cbind(opies,sub,sub.juvs)
X

cor(X)
cor2 <- cor(X, use="complete.obs")


corrplot(cor2)

t<-cor2[,81:ncol(cor2)]
t
corrplot(t)
############################EDITED OPILIO SERIES#############################################
X<-cbind(opies_edit,sub,sub.juvs)
X

cor(X)
cor2 <- cor(X, use="complete.obs")


corrplot(cor2)

t<-cor2[,1:ncol(cor2)]
t
corrplot(t)

################################### 3 year lag###########################################
sub.juvs<-juvs[9:nrow(juvs),]
sub<-Spawners_SC3$EBS_CBfem_SC3[6:(nrow(Spawners_SC3)-3)]
opies_edit2<-opies_edit[1:(nrow(opies_edit)-3),]

X<-cbind(opies_edit2,sub,sub.juvs)
X

cor(X)
cor2 <- cor(X, use="complete.obs")


corrplot(cor2)

t<-cor2[,1:ncol(cor2)]
t
corrplot(t)

################################### 2 year lag###########################################
sub.juvs<-juvs[9:nrow(juvs),]
sub<-Spawners_SC3$EBS_CBfem_SC3[7:(nrow(Spawners_SC3)-2)]
opies_edit2<-opies_edit[2:(nrow(opies_edit)-2),]

X<-cbind(opies_edit2,sub,sub.juvs)
X

cor(X)
cor2 <- cor(X, use="complete.obs")


corrplot(cor2)

t<-cor2[,1:ncol(cor2)]
t
corrplot(t)

################################### 1 year lag###########################################
sub.juvs<-juvs[9:nrow(juvs),]
sub<-Spawners_SC3$EBS_CBfem_SC3[8:(nrow(Spawners_SC3)-1)]
opies_edit2<-opies_edit[3:(nrow(opies_edit)-1),]

X<-cbind(opies_edit2,sub,sub.juvs)
X

cor(X)
cor2 <- cor(X, use="complete.obs")


corrplot(cor2)

t<-cor2[,1:ncol(cor2)]
t
corrplot(t)