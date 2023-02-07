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
library(magick)
library(here)
library(magrittr)
#?nlstools

########## Import data and define variables####
Recruits<-read.csv("data/cb_ebs_pop_juvenile.csv")
Spawners<-read.csv("data/cb_ebs_pop_sc_cs.csv")

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
output2 <- data.frame()
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
  output2 <- rbind(output2,
                data.frame(window_end = window_end[i],
                           log_R_AICc = AICc(mod2)))
}

######################################################################################################
########### Log R only for paper #############################################################################
dev.new()
plot <- output2 %>%
  mutate(
         delta_AICc_log_R = log_R_AICc - min(log_R_AICc)) %>%
  select(-log_R_AICc) %>%
  pivot_longer(cols = -window_end)

ggplot(plot, aes(window_end, value, color = name)) +
  geom_line() +
  geom_point() +
  labs(y = "AICc",x="Analysis window end")

ggsave("./figs/Log_R_autocorrelation_candidate_windows.png", width = 6, height = 4, units = 'in')

AC_window_plot <- image_read("./figs/Log_R_autocorrelation_candidate_windows.png")
map_plot <- image_read("./figs/EBS_BaseMap_forpaper_nowhitespace_resize25pct.jpg") 

final_plot <- image_append(image_scale(c(map_plot, AC_window_plot)), stack = TRUE)
dev.new()
plot(final_plot)

ggsave("./figs/Log_R_autocorrelation_candidate_windows.png", width = 6, height = 4, units = 'in')
############################################################################################
