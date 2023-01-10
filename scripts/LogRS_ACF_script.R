

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


#########################################################################################################
########################### Create recruitment time series for analysis #################################
#########################################################################################################
names(Recruits)

########################### EBS ####################################################################################

rec_30to50<-as.data.frame(as.matrix(cbind(Recruits$SURVEY_YEAR,(Recruits$NUM_MALE_30TO50 + Recruits$NUM_FEMALE_30TO50))))


colnames(rec_30to50)<-c("Year", "EBS_Abun_30to50","E166_Abun_30to50","W166_Abun_30to50")


#########################################################################################################
########################### Create spawner series #######################################################
#########################################################################################################
names(Spawners)


sp_sc3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,(Spawners$NUM_FEMALE_SC3))))




Spawners_SC3<-as.data.frame(as.matrix(cbind(Spawners$SURVEY_YEAR,Spawners$NUM_FEMALE_SC3)))

colnames(Spawners_SC3)<-c("Year", "EBS_SC3")
Spawners_SC3


########################################################################################################################
########################################## ACF analysis of series by stanza ############################################
########################################################################################################################

par(mfrow=c(2,2),cex.lab=1.25,cex.axis=1.25,cex=1.25) #configure axis labels)

########################################################################################################################
########################################## ACF for juvenile indices ####################################################
rec_30to50$Year[4:45] #check years


R<-rec_30to50$EBS_Abun_30to50[4:45]

########################################## Full Timeseries #############################################################
dev.new()
output <- acf(R,main="Juvenile index, 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title

str(output)
output$acf
#?acf()

## loop through different window lengths to evaluate support for 2008 breakpoint ----------------------
library(tidyverse)
library(MuMIn)
theme_set(theme_bw())

# create data frame
dat <- data.frame(year = rec_30to50$Year[4:45],
                  R = rec_30to50$EBS_Abun_30to50[4:45],
                  log_R = log(rec_30to50$EBS_Abun_30to50[4:45]),
                  R_lag7 = c(rep(NA, 7), rec_30to50$EBS_Abun_30to50[4:38]),
                  log_R_lag7 = c(rep(NA, 7), log(rec_30to50$EBS_Abun_30to50[4:38])))

acf(dat$R[dat$year %in% 1978:2008])
acf(dat$log_R[dat$year %in% 1978:2008])

# loop through and fit AR model by different window lengths

length(1985:2019) # 35 years with lag 7 data available

# 0.2 of 35 = 7 years (minimum window length)

window_end <- 1991:2012

# create data frame to capture results
output <- data.frame()

# loop through each candidate window, fit linear model, and record AICc value

for(i in 1:length(window_end)){
  
  # i <- 1
  dat$era <- if_else(dat$year <= window_end[i], "one", "two")
  
  # fit models
  mod1 <- lm(R ~ R_lag7:era, data = dat)
  
  mod2 <- lm(log_R ~ log_R_lag7:era, data = dat)
  
  # record AICc scores
  output <- rbind(output,
                  data.frame(window_end = window_end[i],
                             R_AICc = AICc(mod1),
                             log_R_AICc = AICc(mod2)))

}

# arrange and plot
plot <- output %>%
  mutate(delta_AICc_R = R_AICc - min(R_AICc),
         delta_AICc_log_R = log_R_AICc - min(log_R_AICc)) %>%
  select(-R_AICc, -log_R_AICc) %>%
  pivot_longer(cols = -window_end)

ggplot(plot, aes(window_end, value, color = name)) +
  geom_line() +
  geom_point() +
  labs(y = "AICc",x="Analysis window end")

ggsave("./figs/R_autocorrelation_candidate_windows.png", width = 6, height = 4, units = 'in')

# 1998 is the best-supported window end

dev.new()
plot_acf <- data.frame(Lag = as.factor(1:10),
                       ACF_85_98 = acf(dat$log_R[dat$year %in% 1985:1998])$acf[1:10],
                       ACF_99_19 = acf(dat$log_R[dat$year %in% 1999:2019])$acf[1:10]) %>%
  pivot_longer(cols = -Lag)
  
ggplot(plot_acf, aes(Lag, value, fill = name)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0) +
  ylab("Autocorrelation")

ggsave("./figs/log_R_autocorrelation_time_blocks.png", width = 6, height = 4, units = 'in')

########################################## Divide into two equal stanzas and original from grad school work (1978 to 2008)#
R1<-R[1:21]  #1978-1998
R2<-R[22:42] #1999-2019
R3<-R[1:31]  #1978-2008
R4<-R[32:42]  #2009-2019

acf(R3,main="Juvenile index, 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R4,main="Juvenile index, 2009 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R1,main="Juvenile index, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R2,main="Juvenile index, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

###############################################################################################################################
########################################## Now do ACF for female indices ######################################################
Spawners_SC3$Year
S<-Spawners_SC3$EBS_SC3[4:45]

########################################## Full Timeseries #############################################################

acf(S,main="SC3 females 1978 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)##### using "" removes main title


########################################## Divide into two equal stanzas and Original (1978 to 2008)##############################################
S1<-S[1:21]  #1978-1998
S2<-S[22:42] #1999-2019
S3<-S[1:31]  #1978-2008
S4<-S[32:42]  #2009-2019

acf(S3,main="SC3 females 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S4,main="SC3 females 2009 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S1,main="SC3 females 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(S2,main="SC3 females 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

## evaluate evidence for changing ar(7) characteristics across the time series -----------------------

# first, examine ln(R)
R1<-log(R[1:21])  #1978-1998
R2<-log(R[22:42]) #1999-2019
R3<-log(R[1:31])  #1978-2008
R4<-log(R[32:42])  #2009-2019

acf(R3,main="Juvenile index, 1978 to 2008",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R4,main="Juvenile index, 2009 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R1,main="Juvenile index, 1978 to 1998",cex.lab=1.25,cex.axis=1.25,cex=1.25)
acf(R2,main="Juvenile index, 1999 to 2019",cex.lab=1.25,cex.axis=1.25,cex=1.25)

# not seeing a change from ACF on raw data - seems odd?

recr <- data.frame(year=1978:2019,
                   R=R,
                   R.7=c(rep(NA,7), R[1:35]))

# run a loop calculating AIC over different windows, while keeping 
# at least 15 years in each period

mod.out <- data.frame()

for(i in 1992:2004){
  
  recr$era <- ifelse(recr$year <= i, "early", "late")
  
  mod <- lm(R ~ R.7:era, data=recr)
  
  temp.out <- data.frame(year=i,
                         AIC = MuMIn::AICc(mod))
  
  mod.out <- rbind(mod.out, temp.out)
}

mod.out$delta_AICc <- mod.out$AIC - min(mod.out$AIC)

# add null (stationary) model as comparison

mod.null <- lm(R ~ R.7, data=recr)

# rerun using GLS procedures to allow accounting for autocorrelation
mod.null2<-gls(R ~ R.7, data=recr,correlation=corAR1(),na.action=na.omit)

null <- data.frame(year="NULL",
                   delta_AICc = MuMIn::AICc(mod.null)-min(mod.out$AIC))

library(tidyverse)
theme_set(theme_bw())

ggplot(mod.out, aes(year, delta_AICc)) +
  geom_line() +
  geom_hline(yintercept = null$delta_AICc, lty = 2) +
  labs(title = "AICc for different change points: AR(7) model for recruitment", subtitle = "Dashed line = stationary model") +
  theme(axis.title.x = element_blank()) 

ggsave("./figs/AR7_changepoint_AIC_plot_R.png", width = 6, height = 4, units = 'in')

# best model is 2000
# get model summaries

i <- 2000
recr$era <- ifelse(recr$year <= i, "early", "late")
mod <- lm(R ~ R.7:era, data=recr)
summary(mod)

# rerun using GLS procedures to allow accounting for autocorrelation
mod2<-gls(R ~ R.7:era, data=recr,correlation=corAR1(),na.action=na.omit)
summary(mod2)


# and plot eras
# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

i <- 2000
recr$era <- ifelse(recr$year <= i, "early", "late")

dev.new()
ggplot(recr, aes(R.7, R, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = cb[c(2,6)]) +
  labs(title = "Early = 2000 and earlier", 
       x = "Recruitment year-7",
       y = "Recruitment year 0")

ggsave("./figs/AR7_changepoint_2000_scatter_plot_R.png", width = 5, height = 3, units = 'in')

anova(mod.null, mod) # note this isn't correct n/c df don't account for autocorrelation!

