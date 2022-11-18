library(Hmisc)
library(ggplot2)
library(tidyverse)
library(corrplot)
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

########################################################################################################

dat<-read.csv("./data/envars.csv")
names(dat)
dat<-as.data.frame(dat)
n<-ncol(dat)
m<-nrow(dat)
dat2<-dat[,2:n]
dat3<-dat2[8:m,]
eda.norm(dat3)

cor(dat3)
cor1 <- cor(dat3, use="complete.obs")
corrplot(cor1)

###########################################################################################################################################################
################################ Examine data for correlations between environmental variables using reduced dataset ##############################################
reduced_envars<-read.csv("data/envars_ebs_only.csv")
ncol(reduced_envars)
#cor_vars<-reduced_envars[8:44,2:15]
cor_vars<-reduced_envars[,2:15]
cor2 <- cor(cor_vars, use="complete.obs")
corrplot(cor2)
dimnames(cor_vars)
#########################################################################################################################################
######################## Age 3 to 7 Pacific cod annual estimates ##########################################################################
ggplot(dat, aes(Year, Age3to7Pcodabun)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

#########################Age 3 to 7 Pacific cod annual annual anomalies ####################################################################
ggplot(dat, aes(Year, Age3to7Pcodabun_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

############################################################################################################################################
################################ Pcod 3 year rolling averages on mid year###################################################################
ggplot(dat, aes(Year, Age3to7Pcodabun_RA3_mid)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun_RA3_mid)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ Pcod 3 year rolling averages on mid end year - anomalies#####################################################
ggplot(dat, aes(Year, Age3to7Pcodabun_RA3_mid_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun_RA3_mid_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ Pcod 3 year rolling averages on end year #####################################################
ggplot(dat, aes(Year, Age3to7Pcodabun_RA3_end)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun_RA3_end)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ Pcod 3 year rolling averages on end year - anomalies#####################################################
ggplot(dat, aes(Year, Age3to7Pcodabun_RA3_end_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(Age3to7Pcodabun_RA3_end_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye





##########################################################################################################################################
################################ FHS TBM #####################################################
ggplot(dat, aes(Year, FHS_TBM)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(FHS_TBM)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ FHS TBM - anomalies#####################################################
ggplot(dat, aes(Year, FHS_TBM_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(FHS_TBM_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye



##########################################################################################################################################
################################ EBS Mean NBT #####################################################
ggplot(dat, aes(Year, EBS_mean_NBT)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_mean_NBT)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ EBS Mean NBT - anomalies#####################################################
ggplot(dat, aes(Year, EBS_NBT_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_NBT_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye



##########################################################################################################################################
################################ EBS_NBT_RA3_final_year #####################################################
ggplot(dat, aes(Year, EBS_NBT_RA3_final_year)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_NBT_RA3_final_year)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ EBS_NBT_RA3_final_year - anomalies#####################################################
ggplot(dat, aes(Year, EBS_NBT_RA3_final_year_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_NBT_RA3_final_year_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

##########################################################################################################################################
################################ EBS_NBT_RA3_mid_year #####################################################
ggplot(dat, aes(Year, EBS_NBT_RA3_midyear)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_NBT_RA3_midyear)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye

################################ EBS_NBT_RA3_final_year - anomalies#####################################################
ggplot(dat, aes(Year, EBS_NBT_RA3_midyear_anom)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(EBS_NBT_RA3_midyear_anom)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ AO winter  #####################################################
ggplot(dat, aes(Year, AO_jfm)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(AO_jfm)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ Southeast wind.May-Sep  #####################################################
ggplot(dat, aes(Year, SE.wind.May.Sep)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(SE.wind.May.Sep)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ Northwest wind.May-Sep  #####################################################
ggplot(dat, aes(Year, NW.wind.May.Sep)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(NW.wind.May.Sep)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ PDO - Summer  #####################################################
ggplot(dat, aes(Year, PDO_jja)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(PDO_jja)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ PDO - Winter  #####################################################
ggplot(dat, aes(Year, PDO_djf)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(PDO_djf)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye


##########################################################################################################################################
################################ Winter ice area  #####################################################
ggplot(dat, aes(Year, ice.area.jfma)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

ggplot(dat, aes(ice.area.jfma)) + geom_histogram() #+ facet_wrap(~Year, scales="free") #little better to my eye



##############################################################################################################################################
################################# Envars index for selecting data ###########################################################################
envar$Year          #45 rows...1975 is row 1, 2019 is row 45. 2016 (final hatch year is 42)
# First hatch year will be 1983(9), final will be 2016(42)
#Stock recruit residuals are lage three years from hatch: first will be 1986(12), last will be 2019(45)
# For embryonic year effects envar$[8:41]
# for hatch year effects: envar$[9:42]
# for effects one year after hatch(lag 1): envar$[10:43]
# for effects two years after hatch(lag 2): envar$[11:44]
# for two-year rolling averages with embryonic and hatch year effects : envar$[9:42]
# for two-year rolling averages with hatch year and following year effects : envar[10:43]
# for three-year rolling averages on END YEAR with embryonic year and following 2 years effects : envar$[10:43]
# for three-year rolling averages on END YEAR with hatch year and following 2 years effects : envar$[11:44]
# for three-year rolling averages on MID YEAR with embryonic year and following 2 years effects : envar$[9:42]
# for three-year rolling averages on MID YEAR with hatch year and following 2 years effects : envar$[10:43]



############################################################################################################################################################
####################################### Plot resids versus environmental data ##############################################################################
##############################################################################################################################################################
names(envar)

############################################################################################################################################################
####################################### Arctic oscillation (winter) hatch year ##############################################################################################
dat<-as.data.frame(cbind(envar$AO_jfm[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("AO_jfm","SR_residuals")
dat
ggplot(dat,aes(AO_jfm,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
############################################################################################################################################################
####################################### Arctic oscillation (winter) 1st year ##############################################################################################
dat<-as.data.frame(cbind(envar$AO_jfm[10:43],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("AO_jfm","SR_residuals")
dat
ggplot(dat,aes(AO_jfm,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Arctic oscillation (winter) 2nd year ##############################################################################################
dat<-as.data.frame(cbind(envar$AO_jfm[11:44],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("AO_jfm","SR_residuals")
dat
ggplot(dat,aes(AO_jfm,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
############################################################################################################################################################
####################################### Southeast Wind (May-September) - lagged for hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$SE.wind.May.Sep[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("SE.wind","SR_residuals")
dat
ggplot(dat,aes(SE.wind,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Northwest Wind (May-September) - lagged for hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$NW.wind.May.Sep[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("NW.wind","SR_residuals")
dat
ggplot(dat,aes(NW.wind,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()


############################################################################################################################################################
####################################### Summer PDO - lagged for hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$PDO_jja[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("PDO_jja","SR_residuals")
dat
ggplot(dat,aes(PDO_jja,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()


############################################################################################################################################################
####################################### Winter PDO - lagged for hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$PDO_djf[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("PDO_djf","SR_residuals")
dat
ggplot(dat,aes(PDO_djf,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Winter PDO - lagged for 1st year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$PDO_djf[10:43],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("PDO_djf","SR_residuals")
dat
ggplot(dat,aes(PDO_djf,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Winter PDO - lagged for 2st year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$PDO_djf[11:44],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("PDO_djf","SR_residuals")
dat
ggplot(dat,aes(PDO_djf,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
############################################################################################################################################################
####################################### Winter ice area - lagged for hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$ice.area.jfma[9:42],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("ice.area.jfma","SR_residuals")
dat
ggplot(dat,aes(ice.area.jfma,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()


############################################################################################################################################################
####################################### Winter ice area - lagged for year following hatch year effect ##############################################################################################
dat<-as.data.frame(cbind(envar$ice.area.jfma[10:43],envar$Lag3_SR_residuals[12:45]))
colnames(dat)<-c("ice.area.jfma","SR_residuals")
dat
ggplot(dat,aes(ice.area.jfma,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
#################################################################################################################################################################
####################################### PACIFIC COD ##########################################################################################################
##################################################################################################################################################################

############################################################################################################################################################
####################################### Lag 3 year (release year) predation by Pacific cod ##############################################################################################
dat<-as.data.frame(cbind(envar$Age3to7Pcodabun[9:42],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("Pcod_lag3","SR_residuals")
dat
ggplot(dat,aes(Pcod_lag3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 2 year (one year post release year) predation by Pacific cod ##############################################################################################
dat<-as.data.frame(cbind(envar$Age3to7Pcodabun[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("Pcod_lag2","SR_residuals")
dat
ggplot(dat,aes(Pcod_lag2,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 1 year (two year post release year) predation by Pacific cod ##############################################################################################
dat<-as.data.frame(cbind(envar$Age3to7Pcodabun[11:44],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("Pcod_lag1","SR_residuals")
dat
ggplot(dat,aes(Pcod_lag1,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Pacific cod 3 year RA-end year for predation beginning release year ##############################################################################################
dat<-as.data.frame(cbind(envar$Age3to7Pcodabun_RA3_end[11:44],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("Pcod_RA3","SR_residuals")
dat
ggplot(dat,aes(Pcod_RA3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Pacific cod 3 year RA-mid year for predation beginning release year ##############################################################################################
dat<-as.data.frame(cbind(envar$Age3to7Pcodabun_RA3_mid[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("Pcod_RA3","SR_residuals")
dat
ggplot(dat,aes(Pcod_RA3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()



#################################################################################################################################################################
####################################### FLATHEAD SOLE ##########################################################################################################
##################################################################################################################################################################


############################################################################################################################################################
####################################### Lag 3 year (release year) predation by FHS ##############################################################################################
dat<-as.data.frame(cbind(envar$FHS_TBM[9:42],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("FHS_TBM_lag3","SR_residuals")
dat
ggplot(dat,aes(FHS_TBM_lag3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 2 year (one year post release year) predation by FHS ##############################################################################################
dat<-as.data.frame(cbind(envar$FHS_TBM[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("FHS_TBM_lag2","SR_residuals")
dat
ggplot(dat,aes(FHS_TBM_lag2,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 1 year (two year post release year) predation by FHS ##############################################################################################
dat<-as.data.frame(cbind(envar$FHS_TBM[11:44],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("FHS_TBM_lag1","SR_residuals")
dat
ggplot(dat,aes(FHS_TBM_lag1,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()



#################################################################################################################################################################
####################################### NBT ##########################################################################################################
##################################################################################################################################################################

############################################################################################################################################################
####################################### Lag 4 year (embryonic year)  ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_anom[8:41],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_anom_lag3","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_anom_lag3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()


############################################################################################################################################################
####################################### Lag 3 year (release year)  ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_anom[9:42],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_anom_lag3","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_anom_lag3,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 2 year (one year post release year)  ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_anom[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_anom_lag2","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_anom_lag2,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()

############################################################################################################################################################
####################################### Lag 1 year (two year post release year)  ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_anom[11:44],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_anom_lag1","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_anom_lag1,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()



############################################################################################################################################################
#######################################3 year RA-end year for effect beginning embryonic year ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_RA3_final_year_anom[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_RA3_final_year_anom","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_RA3_final_year_anom,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
#######################################3 year RA-end year for effect beginning release year ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_RA3_final_year_anom[11:44],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_RA3_final_year_anom","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_RA3_final_year_anom,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()
############################################################################################################################################################
####################################### 3 year RA-mid year for effect beginning release year ##############################################################################################
dat<-as.data.frame(cbind(envar$EBS_NBT_RA3_midyear_anom[10:43],envar$Lag3_SR_residuals[12:45]))

colnames(dat)<-c("EBS_NBT_RA3_midyear_anom","SR_residuals")
dat
ggplot(dat,aes(EBS_NBT_RA3_midyear_anom,SR_residuals)) + geom_point() + theme(legend.position = "none") +
  geom_smooth()



###############################################################################################################################################
####################################### All Envars #####################################################

dat4<-read.csv("data/all_envars.csv")
names(dat4)
dat<-as.data.frame(dat4)
n<-ncol(dat4)
m<-nrow(dat4)
n
m
dat5<-dat[,2:n]
dat6<-dat4[8:44,]
eda.norm(dat6)
names(dat6)
cor(dat6)
cor2 <- cor(dat6, use="complete.obs")
dev.new()
corrplot(cor2)

dat6
cod<-dat6$Age3to7Pcodabun
fhs<-dat6$FHS_TBM
#pdo_s<-dat6$PDO_jja
pdo_w<-dat6$PDO_djf
sst_MJ<-dat6$May.July.anom
AO<-dat6$AO_jfm
NBT<-dat6$EBS_NBT_RA3_final_year
dat7<-cbind(cod,fhs,pdo_w,sst_MJ,ao,nbt)
cor(dat7)
cor3 <- cor(dat7, use="complete.obs")
dev.new()
corrplot(cor3)

write.csv(cor3, "")