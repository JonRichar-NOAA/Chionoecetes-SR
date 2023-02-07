
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

sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))
colnames(sp_sc3)<-c("Year", "EBS_sc3fem_abun")

e166_sp_sc3<-as.data.frame(as.matrix(cbind(E166Spawners$SURVEY_YEAR,(E166Spawners$NUM_FEMALE_SC3))))
colnames(e166_sp_sc3)<-c("Year", "E166_sc3fem_abun")

w166_sp_sc3<-as.data.frame(as.matrix(cbind(W166Spawners$SURVEY_YEAR,(W166Spawners$NUM_FEMALE_SC3))))
colnames(w166_sp_sc3)<-c("Year", "W166_sc3fem_abun")

Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3,E166Spawners$NUM_FEMALE_SC3,W166Spawners$NUM_FEMALE_SC3)))
colnames(Spawners_SC3)<-c("Year", "EBS_SC3","E166_SC3","W166_SC3")


######################################### juvenile indices ####################################################
rec_30to50$Year[4:45] #check years
par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

R<-rec_30to50$EBS_Abun_30to50[4:45]


########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
R1<-R[1:21]  #1978-1998
R2<-R[22:42] #1999-2019
R3<-R[1:31]  #1978-2008



###############################################################################################################################
########################################## female indices ######################################################
###############################################################################################################################

###############################################################################################################################
########################################## Shell condition 3 #################################################################
Spawners_SC3$Year
S<-Spawners_SC3$EBS_SC3[4:45]
########################################## Full Timeseries #############################################################


########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008


#########################################################################################################################
################################## Plot all juveniles together ##########################################################
EBS.juv.a<-rec_30to50$EBS_Abun_30to50
EBS.juv<-EBS.juv.a[4:length(EBS.juv.a)]
year<-rec_30to50$Year[4:length(rec_30to50$Year)]
EBS.fem<-Spawners_SC3$EBS_SC3[4:length(Spawners_SC3$EBS_SC3)]
juv.dat<-data.frame(Year=year, abundance = EBS.juv, name ="Juvenile")

fem.dat<-data.frame(Year=year, abundance = EBS.fem, name ="SC3 female")

plot_dat <-data.frame(rbind(juv.dat,fem.dat))

theme_set(theme_bw())

ggplot(plot_dat, aes(Year, abundance)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  theme(axis.title.x = element_blank()) +
  ylab("Abundance (millions)")

ggsave("./figs/juvenile_female_abundance.png",
        width = 6,
        height = 6,
        units = 'in')
getwd()
