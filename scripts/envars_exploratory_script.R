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

dat<-read.csv("data/envars.csv")
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


######################################################################################################
names(envar)
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

cor(dat6)
cor2 <- cor(dat6, use="complete.obs")
corrplot(cor2)

dat6
