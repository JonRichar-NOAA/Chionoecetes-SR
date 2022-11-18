library(Hmisc)
library(ggplot2)
library(tidyverse)
library(corrplot)



###################Final envar dataset ###############
final_dat<-read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv")
names(final_dat)
cor <- cor(final_dat, use="complete.obs")
dev.new()
corrplot(cor)

Bairdi_females<-final_dat$ReproductiveFemales
cod<-final_dat$PCod_RA3
fhs_RA2<-final_dat$FHS_RA2
FHS_lag2<-final_dat$FHS_lag2
opilio_female<-final_dat$Ovig_female_CO
pdo<-final_dat$PDO_RA3
sst_MJ<-final_dat$SST_May_July
AO<-final_dat$AO_RA3
NBT<-final_dat$NBT_3RA
SE_wind<-final_dat$SE.wind
NE_wind<-final_dat$NE.wind
dat<-cbind(Bairdi_females,cod,fhs_RA2,FHS_lag2,opilio_female,pdo,sst_MJ,AO,NBT,SE_wind,NE_wind)

cor(dat, use="complete.obs")
cor2 <- cor(dat, use="complete.obs")
dev.new()
corrplot(cor2)

write.csv(cor2, "data/CorrelationMatrixforFinalEnvironmentalDataset.csv")
