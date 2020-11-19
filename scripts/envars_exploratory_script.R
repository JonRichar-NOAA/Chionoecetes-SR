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
