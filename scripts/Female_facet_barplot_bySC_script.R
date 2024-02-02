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

dat <- read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis_reducedTime_updated.csv", row.names = 1) #Has 1982 female bairdi and corresponding juvenile and environmental data removed

head(dat)

#========================= females dataseriews ==========================================================================
SC2<-as.data.frame(cbind(dat$releaseyear,dat$SC2ReproductiveFemales))
SC2b<-as.data.frame(cbind(dat$releaseyear,dat$SC2RepFem_SameYear))
SC3<-as.data.frame(cbind(dat$releaseyear,dat$SC3ReproductiveFemales))
SC3and4<-as.data.frame(cbind(dat$releaseyear,dat$SC3_SC4_ReproductiveFemales))
SC4<-as.data.frame(cbind(dat$releaseyear,(dat$SC3_SC4_ReproductiveFemales-dat$SC3ReproductiveFemales)))

SC2$SC<-"SC2"
SC2b$SC<-"SC2"
SC3$SC<-"SC3"
SC4$SC<-"SC4"
view(SC2)
view(SC3)
view(SC4)

colnames(SC2)<-c("Year","NumFem","SC")
colnames(SC2b)<-c("Year","NumFem","SC")
colnames(SC3)<-c("Year","NumFem","SC")
colnames(SC4)<-c("Year","NumFem","SC")

plot_dat<-rbind(SC2,SC3,SC4)


plot_dat_SC2_sameyr<-rbind(SC2b,SC3,SC4)

view(plot_dat)
view(plot_dat_SC2_sameyr)
#######################################################################################################################
################################ Facet plot ###########################################################################

facet_abun<- ggplot(data=plot_dat, aes(x=SC, y=NumFem))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)",x="Shell condition")+
  facet_wrap(~Year, scales = "free_y") +
  labs(y = "Abundance(millions)", main = "Female Tanners") 

ggsave(plot = facet_abun, "./Female_bairdi_by SC_barplot_facet.jpeg", height=10, width=7, units="in")


facet_abun2<- ggplot(data=plot_dat, aes(x=SC, y=NumFem))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)",x="Shell condition")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Abundance(millions)", main = "Female Tanners") 
ggsave(plot = facet_abun2, "./Female_bairdi_by SC_barplot_facet_fixed_yrange.jpeg", height=10, width=7, units="in")

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

SC_total<-as.data.frame(cbind(dat$releaseyear,
                              dat$SC2RepFem_SameYear,
                              dat$SC3ReproductiveFemales,
                              dat$SC3_SC4_ReproductiveFemales-dat$SC3ReproductiveFemales,
                              dat$SC3_SC4_ReproductiveFemales,
                              dat$SC2RepFem_SameYear/(dat$SC2RepFem_SameYear+dat$SC3_SC4_ReproductiveFemales),
                              dat$SC3ReproductiveFemales/(dat$SC2RepFem_SameYear+dat$SC3_SC4_ReproductiveFemales),
                              (dat$SC3_SC4_ReproductiveFemales-dat$SC3ReproductiveFemales)/(dat$SC2RepFem_SameYear+dat$SC3_SC4_ReproductiveFemales),
                              dat$SC3_SC4_ReproductiveFemales/(dat$SC2RepFem_SameYear+dat$SC3_SC4_ReproductiveFemales)))




colnames(SC_total)<-c("Year","SC2Fem","SC3Fem","SC4Fem","SC3_SC4Fem","propSC2","propSC3","propSC4","propSC3_SC4")
view(SC_total)
SC_total%>% select(Year,propSC2,propSC3,propSC4)->prop_dat
prop_dat

prop_dat%>%pivot_longer(cols=propSC2:propSC4)->piv_dat
colnames(piv_dat)<-c("Year","Class","Proportion")

facet_abun2<- ggplot(data=piv_dat, aes(y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Proportion SC class", main = "Female Tanners") 
ggsave(plot = facet_abun2, "./Female_bairdi_by SC_barplot_facet_fixed_yrange_sameyr_SC2.jpeg", height=10, width=7, units="in")

ggplot(data=piv_dat, aes(x=Year,y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  facet_wrap(~Year, scales = "free_x") +
  labs(y = "Proportion SC class", main = "Female Tanners") 

ggplot(data=piv_dat, aes(x=Year,y=Proportion,fill = Class))+
  geom_bar(stat="identity") +
  labs(y="Abundance(millions)")+
  labs(y = "Proportion SC class", main = "Female Tanners") 
