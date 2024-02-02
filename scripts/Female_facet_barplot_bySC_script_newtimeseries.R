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
library(voxel)
#?nlstools

#======================== Import data ================================================================
getwd()

dat <- read.csv("data/Female_Tanner_Crab_Series_for_analysis.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

#========================= females dataseriews ==========================================================================
SC0<-as.data.frame(cbind(dat$Year,dat$EBS_SC0))
SC1<-as.data.frame(cbind(dat$Year,dat$EBS_SC1))
SC2<-as.data.frame(cbind(dat$Year,dat$EBS_SC2))
SC3<-as.data.frame(cbind(dat$Year,dat$EBS_SC3))
SC4<-as.data.frame(cbind(dat$Year,dat$EBS_SC4))
SC5<-as.data.frame(cbind(dat$Year,dat$EBS_SC5))
SC_tot<-as.data.frame(cbind(dat$Year,dat$EBS_SC_TOTAL))

SC0$SC<-"SC0"
SC1$SC<-"SC1"
SC2$SC<-"SC2"
SC3$SC<-"SC3"
SC4$SC<-"SC4"
SC5$SC<-"SC5"

#view(SC2)
#view(SC3)
#view(SC4)

colnames(SC0)<-c("Year","NumFem","SC")
colnames(SC1)<-c("Year","NumFem","SC")
colnames(SC2)<-c("Year","NumFem","SC")
colnames(SC3)<-c("Year","NumFem","SC")
colnames(SC4)<-c("Year","NumFem","SC")
colnames(SC5)<-c("Year","NumFem","SC")

plot_dat_SC2_sameyr<-rbind(SC0,SC1,SC2,SC3,SC4,SC5)


#plot_dat<-rbind(SC2b,SC3,SC4)

#view(plot_dat)
view(plot_dat_SC2_sameyr)
#######################################################################################################################
################################ Facet plot ###########################################################################

#facet_abun<- ggplot(data=plot_dat, aes(x=SC, y=NumFem))+
 # geom_bar(stat="identity") +
  #labs(y="Abundance(millions)",x="Shell condition")+
  #facet_wrap(~Year, scales = "free_y") +
  #labs(y = "Abundance(millions)", main = "Female Tanners") 

#ggsave(plot = facet_abun, "./Female_bairdi_by SC_barplot_facet.jpeg", height=10, width=7, units="in")


#facet_abun2<- ggplot(data=plot_dat, aes(x=SC, y=NumFem))+
#  geom_bar(stat="identity") +
#  labs(y="Abundance(millions)",x="Shell condition")+
#  facet_wrap(~Year, scales = "free_x") +
#  labs(y = "Abundance(millions)", main = "Female Tanners") 
#ggsave(plot = facet_abun2, "./Female_bairdi_by SC_barplot_facet_fixed_yrange.jpeg", height=10, width=7, units="in")

#######################################################################################################################
################################ Facet plots using same year SC2 ###########################################################################

facet_abun<- ggplot(data=plot_dat_SC2_sameyr, aes(x=SC, y=NumFem))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)",x="Shell condition")+
  facet_wrap(~Year, scales = "free_y") +
  labs(y = "Abundance(millions)", main = "Female Tanners") 

ggsave(plot = facet_abun, "./Female_bairdi_by SC_barplot_facet_sameyear_SC2.jpeg", height=10, width=7, units="in")


facet_abun2<- ggplot(data=plot_dat_SC2_sameyr, aes(x=SC, y=NumFem))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)",x="Shell condition")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Abundance(millions)", main = "Female Tanners") 
ggsave(plot = facet_abun2, "./Female_bairdi_by SC_barplot_facet_fixed_yrange_sameyr_SC2.jpeg", height=10, width=7, units="in")
#================================================================================================================
names(dat)

SC2<-as.data.frame(cbind(dat$releaseyear,dat$SC2ReproductiveFemales))
SC2b<-as.data.frame(cbind(dat$releaseyear,dat$SC2RepFem_SameYear))
SC3<-as.data.frame(cbind(dat$releaseyear,dat$SC3ReproductiveFemales))
SC4<-as.data.frame(cbind(dat$releaseyear,dat$SC3_SC4_ReproductiveFemales-dat$SC3ReproductiveFemales))

SC0<-as.data.frame(cbind(dat$Year,dat$EBS_SC0))
SC1<-as.data.frame(cbind(dat$Year,dat$EBS_SC1))
SC2<-as.data.frame(cbind(dat$Year,dat$EBS_SC2))
SC3<-as.data.frame(cbind(dat$Year,dat$EBS_SC3))
SC4<-as.data.frame(cbind(dat$Year,dat$EBS_SC4))
SC5<-as.data.frame(cbind(dat$Year,dat$EBS_SC5))
SC_tot<-as.data.frame(cbind(dat$Year,dat$EBS_SC_TOTAL))

SC_total<-as.data.frame(cbind(dat$Year,
                              dat$EBS_SC0,
                              dat$EBS_SC1,
                              dat$EBS_SC2,
                              dat$EBS_SC3,
                              dat$EBS_SC4,
                              dat$EBS_SC5,
                              (dat$EBS_SC3+dat$EBS_SC4),
                              dat$EBS_SC_TOTAL,
                              dat$EBS_SC2/(dat$EBS_SC_TOTAL),
                              dat$EBS_SC3/(dat$EBS_SC_TOTAL),
                              dat$EBS_SC4/(dat$EBS_SC_TOTAL),
                              dat$EBS_SC5/(dat$EBS_SC_TOTAL),
                              (dat$EBS_SC3+dat$EBS_SC4)/(dat$EBS_SC_TOTAL)))
                          

colnames(SC_total)<-c("Year","SC0Fem","SC1Fem","SC2Fem","SC3Fem","SC4Fem","SC5Fem","SC3_SC4Fem","Total_Fem","propSC2","propSC3","propSC4","propSC5","propSC3_SC4")
#view(SC_total)
SC_total%>% select(Year,propSC2,propSC3,propSC4,propSC5)->prop_dat
prop_dat

prop_dat%>%pivot_longer(cols=propSC2:propSC5)->piv_dat
piv_dat
colnames(piv_dat)<-c("Year","Class","Proportion")
piv_dat

facet_abun2<- ggplot(data=piv_dat, aes(y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Proportion SC class", main = "Female Tanners") 

ggsave(plot = facet_abun2, "./Female_bairdi_by SC_barplot_fixed_yrange_sameyr_SC2.jpeg", height=10, width=7, units="in")

ggplot(data=piv_dat, aes(x=Year,y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Proportion SC class", main = "Female Tanners") 

ggplot(data=piv_dat, aes(x=Year,y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  labs(y = "Proportion SC class", main = "Female Tanners") 
