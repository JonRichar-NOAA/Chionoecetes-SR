library(Hmisc)
library(ggplot2)
library(tidyverse)
library(corrplot)



###################Final envar dataset ###############
final_dat<-read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis_reducedTime_updated.csv")
names(final_dat)
final_dat
cor <- cor(final_dat, use="complete.obs")
dev.new()
corrplot(cor)

Bairdi_females<-final_dat$ReproductiveFemales
cod<-final_dat$Pcod_lag1
cod2<-final_dat$PCod_RA2
cod3<-final_dat$PCod_RA3
fhs2<-final_dat$FHS_RA2
FHS<-final_dat$FHS_lag2
opilio_female<-final_dat$Ovig_female_CO
pdo3<-final_dat$PDO_RA3
pdo2<-final_dat$PDO_RA2
sst_MJ<-final_dat$SST_May_July
AO3<-final_dat$AO_RA3
AO2<-final_dat$AO_RA2
NBT3<-final_dat$NBT_3RA
SE_wind<-final_dat$SE.wind
NW_wind<-final_dat$NW.wind

dat<-cbind(Bairdi_females,cod,cod2,cod3,fhs2,FHS,opilio_female,pdo3,pdo2,sst_MJ,AO3,AO2,NBT3,SE_wind,NW_wind)

cor(dat, use="complete.obs")
cor2 <- cor(dat, use="complete.obs")
dev.new()
corrplot(cor2)

write.csv(cor2, "data/CorrelationMatrixforFinalEnvironmentalDataset.csv")


plot(opilio_female~AO3)
plot(opilio_female~AO2)
