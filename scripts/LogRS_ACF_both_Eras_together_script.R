

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
library(MuMIn)
#?nlstools

########## Import data and define variables####


Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")
names(Recruits)
names(Spawners)


#########################################################################################################
########################### Create recruitment time series for analysis #################################
#########################################################################################################
names(Recruits)

########################### EBS ####################################################################################

rec_30to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO50 + Recruits$NUM_FEMALE_30TO50))))


colnames(rec_30to50)<-c("Year", "EBS_Abun_30to50")


#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################
names(Spawners)


sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))
Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3)))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3")
Spawners_SC3



########################################################################################################################
########################################## Create juvenile series ####################################################
rec_30to50$Year[4:45] #check years
R<-rec_30to50$EBS_Abun_30to50[4:45]

########################################## Full Timeseries #############################################################
log_R <-log(rec_30to50$EBS_Abun_30to50[4:45])
rec_30to50$EBS_Abun_30to50[4:45]
## loop through different window lengths to evaluate support for 2008 breakpoint ----------------------

theme_set(theme_bw())

# create data frame
dat <- data.frame(year = rec_30to50$Year[4:45],
                  R = rec_30to50$EBS_Abun_30to50[4:45],
                  log_R = log(rec_30to50$EBS_Abun_30to50[4:45]))





dev.new()
plot_acf <- data.frame(Lag = as.factor(1:10),
                       ACF_1981_1998 = acf(dat$log_R[dat$year %in% 1981:1998])$acf[1:10],
                       ACF_1999_2019 = acf(dat$log_R[dat$year %in% 1999:2019])$acf[1:10]) %>%
  pivot_longer(cols = -Lag)
  
ggplot(plot_acf, aes(Lag, value, fill = name)) +
  geom_col(position = "dodge",width =0.5) +
  geom_hline(yintercept = 0) +
  ylab("Autocorrelation")

ggsave("./figs/log_R_ACF_both_Stanzas_together.png", width = 6, height = 4, units = 'in')

