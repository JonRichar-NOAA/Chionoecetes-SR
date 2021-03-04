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
#?nlstools

########## Import data and define variables####

Lim_Dat<-read.csv("data/EBS_Crab_and_envar_data_for_analysis.csv")
All_Dat<-read.csv("data/EBS_Crab_and_envar_data_full_extent_for_analysis.csv") #Has 1982 female bairdi and corresponding juvenile and environmental data removed
SR_Dat<-read.csv("data/SR_data.csv")
All_Dat
names(All_Dat)
names(SR_Dat)
###### Create data variables
r<-All_Dat$Lag3_recruits
s<-All_Dat$ReproductiveFemales
lnRS<-All_Dat$logRS
Fem_CO<-All_Dat$Ovig_female_CO
cod<-All_Dat$Pcod_lag1
cod2<-All_Dat$PCod_RA2
cod3<-All_Dat$PCod_RA3
fhs<-All_Dat$FHS_lag2
fhs2<-All_Dat$FHS_RA2
nbt<-All_Dat$NBT_3RA
wind<-All_Dat$SE.wind
ao<-All_Dat$AO_jfm
ao2<-All_Dat$AO_RA2
ao3<-All_Dat$AO_RA3
pdo<-All_Dat$PDO_djf
pdo2<-All_Dat$PDO_RA2
pdo3<-All_Dat$PDO_RA3
sst<-All_Dat$SST_May_July
year<-All_Dat$releaseyear


year.o<-All_Dat$releaseyear[3:length(All_Dat$releaseyear)]   # matches with female opilio annual estimates
year.fhs<-All_Dat$releaseyear[4:length(All_Dat$releaseyear)] # matches flathead sole TBM annual estimates
year.fhs2<-All_Dat$releaseyear[5:length(All_Dat$releaseyear)]  #matches with 2 year rolling average of flathead sole TBM
year
year.o
year.fhs
year.fhs2
length(fhs2)
length(year.fhs2)
fhs2
##############################################################################################################################################
####################################### Now use GAM models ##################################################################################

mod1 <- gam(r ~ s(s, k = 4),
            data = dat, family="Gamma",method = "REML")

mod2 <- gam(r ~ s(Fem_CO, k = 4),
            data = dat, family="Gamma",method = "REML")

mod3 <- gam(r ~ s(cod, k = 4),
            data = dat, family="Gamma",method = "REML")

mod4 <- gam(r ~ s(cod2, k = 4),
            data = dat, family="Gamma",method = "REML")

mod5 <- gam(r ~ s(cod3, k = 4),
            data = dat, family="Gamma",method = "REML")

mod6 <- gam(r ~ s(fhs, k = 4),
            data = dat, family="Gamma",method = "REML")

mod7 <- gam(r ~ s(fhs2, k = 4),
            data = dat, family="Gamma",method = "REML")

mod8 <- gam(r ~ s(nbt, k = 4),
            data = dat, family="Gamma",method = "REML")

mod9 <- gam(r ~ s(wind, k = 4),
            data = dat, family="Gamma",method = "REML")

mod10 <- gam(r ~ s(ao, k = 4),
            data = dat, family="Gamma",method = "REML")

mod11 <- gam(r ~ s(ao2, k = 4),
             data = dat, family="Gamma",method = "REML")

mod12 <- gam(r ~ s(ao3, k = 4),
             data = dat, family="Gamma",method = "REML")

mod13 <- gam(r ~ s(pdo, k = 4),
            data = dat, family="Gamma",method = "REML")

mod14 <- gam(r ~ s(pdo2, k = 4),
             data = dat, family="Gamma",method = "REML")

mod15 <- gam(r ~ s(pdo3, k = 4),
             data = dat, family="Gamma",method = "REML")

mod16 <- gam(r ~ s(sst, k = 4),
            data = dat, family="Gamma",method = "REML")

######################################################################################################################################################################
################################################ Plot GAM models ############################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod1 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod2, resid=T, pch=16, shade=T, rug=F, main = "Opilio females")

plot(mod3, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod")
plot(mod4, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod - 2 year rolling average")
plot(mod5, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod - 3 year rolling average")

plot(mod6, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole")
plot(mod7, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole - 2 year rolling average")

plot(mod8, resid=T, pch=16, shade=T, rug=F, main = "NBT")
plot(mod9, resid=T, pch=16, shade=T, rug=F, main = "SE wind")

plot(mod10, resid=T, pch=16, shade=T, rug=F, main = "AO")
plot(mod11, resid=T, pch=16, shade=T, rug=F, main = "AO - 2 year rolling average")
plot(mod12, resid=T, pch=16, shade=T, rug=F, main = "AO - 3 year rolling average")

plot(mod13, resid=T, pch=16, shade=T, rug=F, main = "PDO")
plot(mod14, resid=T, pch=16, shade=T, rug=F, main = "PDO - 2 year rolling average")
plot(mod15, resid=T, pch=16, shade=T, rug=F, main = "PDO - 3 year rolling average")

plot(mod16, resid=T, pch=16, shade=T, rug=F, main = "SST")

r1<-resid(mod1)
r2<-resid(mod2)
r3<-resid(mod3)
r4<-resid(mod4)
r5<-resid(mod5)
r6<-resid(mod6)
r7<-resid(mod7)
r8<-resid(mod8)
r9<-resid(mod9)
r10<-resid(mod10)
r11<-resid(mod11)
r12<-resid(mod12)
r13<-resid(mod13)
r14<-resid(mod14)
r15<-resid(mod15)
r16<-resid(mod16)

dev.new()
par(mfrow=c(3,3))

plot(r1~year,pch=16 , main = "Reproductive female Bairdi")
plot(r2~year.o,pch=16 , main = "Reproductive female Opilio")
plot(r3~year,pch=16 , main = "Pacific cod")
plot(r4~year,pch=16 , main = "Pacific cod - 2 year rolling average")
plot(r5~year,pch=16 , main = "Pacific cod - 3 year rolling average")

plot(r6~year.fhs,pch=16 , main = "Flathead sole")
plot(r7~year.fhs2,pch=16 , main = "Flathead sole - 2 year rolling average")
plot(r8~year,pch=16 , main = "NBT")

dev.new()
par(mfrow=c(3,3))

plot(r9~year,pch=16 , main = "SE wind")
plot(r10~year,pch=16 , main = "AO")
plot(r11~year,pch=16 , main = "AO - 2 year rolling average")
plot(r12~year,pch=16 , main = "AO - 3 year rolling average")

plot(r13~year,pch=16 , main = "PDO")
plot(r14~year,pch=16 , main = "PDO - 2 year rolling average")
plot(r15~year,pch=16 , main = "PDO - 3 year rolling average")
plot(r16~year,pch=16 , main = "SST")


length(r2)


################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit1<-loess(r1~year,span=0.4)
r.fit2<-loess(r2~year.o,span=0.4)
r.fit3<-loess(r3~year,span=0.4)
r.fit4<-loess(r4~year,span=0.4)
r.fit5<-loess(r5~year,span=0.4)
r.fit6<-loess(r6~year.fhs,span=0.4)
r.fit7<-loess(r7~year.fhs2,span=0.4)
r.fit8<-loess(r8~year,span=0.4)
r.fit9<-loess(r9~year,span=0.4)

r.fit10<-loess(r10~year,span=0.4)
r.fit11<-loess(r11~year,span=0.4)
r.fit12<-loess(r12~year,span=0.4)
r.fit13<-loess(r13~year,span=0.4)
r.fit14<-loess(r14~year,span=0.4)
r.fit15<-loess(r15~year,span=0.4)
r.fit16<-loess(r16~year,span=0.4)
dev.new()
par(mfrow=c(3,3))





length(year)

plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female Bairdi")
lines(year,fitted(r.fit1),col=4)

plot(r.fit2,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female opilio")
lines(year.o,fitted(r.fit2),col=4)

plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod")
lines(year,fitted(r.fit3),col=4)

plot(r.fit4,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod - 2 year rolling average")
lines(year,fitted(r.fit4),col=4)

plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod - 3 year rolling average")
lines(year,fitted(r.fit5),col=4)

plot(r.fit6,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole TBM")
lines(year.fhs,fitted(r.fit6),col=4)

plot(r.fit7,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole TBM - 2 year rolling average")
lines(year.fhs2,fitted(r.fit7),col=4)

plot(r.fit8,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="NBT - 3 year rolling average")
lines(year,fitted(r.fit8),col=4)

plot(r.fit9,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="SE wind")
lines(year,fitted(r.fit9),col=4)

plot(r.fit10,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="AO")
lines(year,fitted(r.fit10),col=4)


plot(r.fit11,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="AO - 2 year rolling average")
lines(year,fitted(r.fit11),col=4)

plot(r.fit12,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="AO - 3 year rolling average")
lines(year,fitted(r.fit12),col=4)

plot(r.fit13,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="PDO")
lines(year,fitted(r.fit13),col=4)

plot(r.fit14,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="PDO - 2 year rolling average")
lines(year,fitted(r.fit14),col=4)

plot(r.fit15,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="PDO - 3 year rolling average")
lines(year,fitted(r.fit15),col=4)

plot(r.fit16,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="May-July SST")
lines(year,fitted(r.fit16),col=4)

##############################################################################################################################################
####################################### Now use GAM models with log(R/S) as response##################################################################################

mod20 <- gam(lnRS ~ s(s, k = 4),
             data = dat,method = "REML")

mod21 <- gam(lnRS ~ s(Fem_CO, k = 4),
             data = dat, method = "REML")

mod22 <- gam(lnRS ~ s(cod, k = 4),
             data = dat, method = "REML")


mod23 <- gam(lnRS ~ s(cod2, k = 4),
             data = dat, method = "REML")


mod24 <- gam(lnRS ~ s(cod3, k = 4),
             data = dat, method = "REML")

mod25 <- gam(lnRS ~ s(fhs, k = 4),
             
             data = dat, method = "REML")
mod26 <- gam(lnRS ~ s(fhs2, k = 4),
             
             data = dat, method = "REML")

mod27 <- gam(lnRS ~ s(nbt, k = 4),
             data = dat, method = "REML")

mod28 <- gam(lnRS ~ s(wind, k = 4),
             data = dat, method = "REML")

mod29 <- gam(lnRS ~ s(ao, k = 4),
             data = dat, method = "REML")

mod30 <- gam(lnRS ~ s(ao2, k = 4),
             data = dat, method = "REML")

mod31 <- gam(lnRS ~ s(ao3, k = 4),
             data = dat, method = "REML")

mod32 <- gam(lnRS ~ s(pdo, k = 4),
             data = dat,method = "REML")

mod33 <- gam(lnRS ~ s(pdo2, k = 4),
             data = dat,method = "REML")

mod34 <- gam(lnRS ~ s(pdo3, k = 4),
             data = dat, method = "REML")

mod35 <- gam(lnRS ~ s(sst, k = 4),
             data = dat, method = "REML")

###########################################################################################################################################################################
################################################ Model summaries ##########################################################################################################

summary(mod20)
summary(mod21)
summary(mod22)

summary(mod23)
summary(mod24)
summary(mod25)

summary(mod26)
summary(mod27)
summary(mod28)

summary(mod29)

summary(mod30)
summary(mod31)
summary(mod32)
summary(mod33)
summary(mod34)
summary(mod35)
######################################################################################################################################################################
################################################ Plot GAM models ############################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod20 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod21, resid=T, pch=16, shade=T, rug=F, main = "Opilio females")
plot(mod22, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod")
plot(mod23, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod - 2 year rolling average")
plot(mod24, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod - 3 year rolling average")
plot(mod25, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole TBM")
plot(mod26, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole TBM - 2 year rolling average")
plot(mod27, resid=T, pch=16, shade=T, rug=F, main = "NBT - 3 year rolling average")
plot(mod28, resid=T, pch=16, shade=T, rug=F, main = "SE Wind")

dev.new()
par(mfrow=c(3,3))

plot(mod29, resid=T, pch=16, shade=T, rug=F, main = "AO")
plot(mod30, resid=T, pch=16, shade=T, rug=F, main = "AO - 2 year rolling average")
plot(mod31, resid=T, pch=16, shade=T, rug=F, main = "AO - 3 year rolling average")
plot(mod32, resid=T, pch=16, shade=T, rug=F, main = "PDO")
plot(mod33, resid=T, pch=16, shade=T, rug=F, main = "PDO - 2 year rolling average")
plot(mod34, resid=T, pch=16, shade=T, rug=F, main = "PDO -3 year rolling average")
plot(mod35, resid=T, pch=16, shade=T, rug=F, main = "SST")

r20<-resid(mod20)
r21<-resid(mod21)
r22<-resid(mod22)
r23<-resid(mod23)
r24<-resid(mod24)
r25<-resid(mod25)
r26<-resid(mod26)
r27<-resid(mod27)
r28<-resid(mod28)
r29<-resid(mod29)

r30<-resid(mod30)
r31<-resid(mod31)
r32<-resid(mod32)
r33<-resid(mod33)
r34<-resid(mod34)
r35<-resid(mod35)

dev.new()
par(mfrow=c(3,3))

plot(r20~year,pch=16 , main = "Reproductive female Bairdi")
plot(r21~year.o,pch=16 , main = "Reproductive female Opilio")
plot(r22~year,pch=16 , main = "Pacific cod")
plot(r23~year,pch=16 , main = "Pacific cod - 2 year rolling average")
plot(r24~year,pch=16 , main = "Pacific cod - 3 year rolling average")
plot(r25~year.fhs,pch=16 , main = "Flathead sole TBM")
plot(r26~year.fhs2,pch=16 , main = "Flathead sole TBM - 2 year rolling average")
plot(r27~year,pch=16 , main = "NBT - 3 year rolling average")
plot(r28~year,pch=16 , main = "SE wind")
plot(r29~year,pch=16 , main = "AO")
plot(r30~year,pch=16 , main = "AO - 2 year rolling average")
plot(r31~year,pch=16 , main = "AO - 3 year rolling average")
plot(r32~year,pch=16 , main = "PDO")
plot(r33~year,pch=16 , main = "PDO - 2 year rolling average")
plot(r34~year,pch=16 , main = "PDO - 3 year rolling average")
plot(r35~year,pch=16 , main = "SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit20<-loess(r20~year,span=0.4)
r.fit21<-loess(r21~year.o,span=0.4)
r.fit22<-loess(r22~year,span=0.4)
r.fit23<-loess(r23~year,span=0.4)
r.fit24<-loess(r24~year,span=0.4)
r.fit25<-loess(r25~year.fhs,span=0.4)
r.fit26<-loess(r26~year.fhs2,span=0.4)
r.fit27<-loess(r27~year,span=0.4)
r.fit28<-loess(r28~year,span=0.4)
r.fit29<-loess(r29~year,span=0.4)
r.fit30<-loess(r30~year,span=0.4)
r.fit31<-loess(r31~year,span=0.4)
r.fit32<-loess(r32~year,span=0.4)
r.fit33<-loess(r33~year,span=0.4)
r.fit34<-loess(r34~year,span=0.4)
r.fit35<-loess(r35~year,span=0.4)
dev.new()
par(mfrow=c(3,3))


plot(r.fit20,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi")
lines(year,fitted(r.fit20),col=4)

length(year.o)
plot(r.fit21,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female opilio")
lines(year.o,fitted(r.fit21),col=4)


plot(r.fit22,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific cod")
lines(year,fitted(r.fit22),col=4)

plot(r.fit23,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific cod - 2 year rolling average")
lines(year,fitted(r.fit23),col=4)

plot(r.fit24,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific cod - 3 year rolling average")
lines(year,fitted(r.fit24),col=4)

plot(r.fit25,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Flathead sole")
lines(year.fhs,fitted(r.fit25),col=4)

plot(r.fit26,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Flathead sole - 2 year rolling average")
lines(year.fhs2,fitted(r.fit26),col=4)

plot(r.fit27,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="NBT - 3 year rolling average")
lines(year,fitted(r.fit27),col=4)

plot(r.fit28,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="SE wind")
lines(year,fitted(r.fit28),col=4)

plot(r.fit29,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Arctic Oscillation")
lines(year,fitted(r.fit29),col=4)

plot(r.fit30,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Arctic Oscillation - 2 year rolling average")
lines(year,fitted(r.fit30),col=4)

plot(r.fit31,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Arctic Oscillation - 3 year rolling average")
lines(year,fitted(r.fit31),col=4)

plot(r.fit32,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit32),col=4)

plot(r.fit33,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific Decadal Oscillation - 2 year rolling average")
lines(year,fitted(r.fit33),col=4)

plot(r.fit34,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific Decadal Oscillation - 3 year rolling average")
lines(year,fitted(r.fit34),col=4)

plot(r.fit35,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="May-July SST")
lines(year,fitted(r.fit35),col=4)

#########################################################################################################################################################################################################
####################################### Now use GAM models with log(R/S) as response and with additional environmental variable ##################################################################################

mod40 <- gam(lnRS ~ s(s, k = 4),
             data = dat, method = "REML")

mod41 <- gam(lnRS ~s(s, k = 4)+ s(Fem_CO, k = 4),
             data = dat, method = "REML")

mod42 <- gam(lnRS ~ s(s, k = 4)+ s(cod, k = 4),
             data = dat, method = "REML")

mod43 <- gam(lnRS ~ s(s, k = 4)+ s(cod2, k = 4),
             data = dat, method = "REML")

mod44 <- gam(lnRS ~ s(s, k = 4)+ s(cod3, k = 4),
             data = dat, method = "REML")

mod45 <- gam(lnRS ~ s(s, k = 4)+ s(fhs, k = 4),
             data = dat, method = "REML")

mod46 <- gam(lnRS ~ s(s, k = 4)+ s(fhs2, k = 4),
             data = dat, method = "REML")

mod47 <- gam(lnRS ~ s(s, k = 4)+ s(nbt, k = 4),
             data = dat, method = "REML")

mod48 <- gam(lnRS ~ s(s, k = 4)+ s(wind, k = 4),
             data = dat, method = "REML")

mod49 <- gam(lnRS ~ s(s, k = 4)+ s(ao, k = 4),
             data = dat, method = "REML")

mod50 <- gam(lnRS ~ s(s, k = 4)+ s(ao2, k = 4),
             data = dat, method = "REML")

mod51 <- gam(lnRS ~ s(s, k = 4)+ s(ao3, k = 4),
             data = dat, method = "REML")

mod52 <- gam(lnRS ~ s(s, k = 4)+ s(pdo, k = 4),
             data = dat, method = "REML")

mod53 <- gam(lnRS ~ s(s, k = 4)+ s(pdo2, k = 4),
             data = dat, method = "REML")

mod54 <- gam(lnRS ~ s(s, k = 4)+ s(pdo3, k = 4),
             data = dat, method = "REML")

mod55 <- gam(lnRS ~ s(s, k = 4)+ s(sst, k = 4),
             data = dat, method = "REML")

dwtest(mod40)
dwtest(mod41)
dwtest(mod42)
dwtest(mod43)
dwtest(mod44)
dwtest(mod45)
dwtest(mod46)
dwtest(mod47)
dwtest(mod48)
dwtest(mod49)
dwtest(mod50)
dwtest(mod51)
dwtest(mod52)
dwtest(mod53)
dwtest(mod54)
dwtest(mod55)


summary(mod40)
summary(mod41)
summary(mod42)
summary(mod43)
summary(mod44)
summary(mod45)
summary(mod46)
summary(mod47)
summary(mod48)
summary(mod49)
summary(mod50)
summary(mod51)
summary(mod52)
summary(mod53)
summary(mod54)
summary(mod55)
####################################################################################################################################################################################################
################################################ Plot GAM models ###################################################################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod40 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod41, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Opilio females")
plot(mod42, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Pacific cod")
plot(mod43, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA2 Pacific cod")

dev.new()
par(mfrow=c(3,3))
plot(mod44, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA3 Pacific cod")
plot(mod45, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Flathead sole")
plot(mod46, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA2 Flathead sole")
plot(mod47, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and NBT")

# each of these models except first will have two plots (one for each variable), so need to set up second graphics window
dev.new()
par(mfrow=c(3,3))

plot(mod48, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and SE wind")
plot(mod49, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and AO")
plot(mod50, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA2 AO")
plot(mod51, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA3 AO")

# each of these models except first will have two plots (one for each variable), so need to set up second graphics window
dev.new()
par(mfrow=c(3,3))

plot(mod52, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and PDO")
plot(mod53, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA2 PDO")
plot(mod54, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and RA3 PDO")
plot(mod55, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and SST")

r40<-resid(mod40)
r41<-resid(mod41)
r42<-resid(mod42)
r43<-resid(mod43)
r44<-resid(mod44)
r45<-resid(mod45)
r46<-resid(mod46)
r47<-resid(mod47)
r48<-resid(mod48)

r49<-resid(mod49)
r50<-resid(mod50)
r51<-resid(mod51)
r52<-resid(mod52)
r53<-resid(mod53)
r54<-resid(mod54)
r55<-resid(mod55)


dev.new()
par(mfrow=c(3,3))

plot(r40~year,pch=16 , main = "Reproductive female Bairdi")
acf(r40)
plot(r41~year.o,pch=16 , main = "Reproductive female Bairdi+Reproductive female Opilio")
plot(r42~year,pch=16 , main = "Reproductive female Bairdi+Pacific cod")
plot(r43~year,pch=16 , main = "Reproductive female Bairdi+RA2 Pacific cod")
plot(r44~year,pch=16 , main = "Reproductive female Bairdi+RA3 Pacific cod")
plot(r45~year.fhs,pch=16 , main = "Reproductive female Bairdi+Flathead sole")
plot(r46~year.fhs2,pch=16 , main = "Reproductive female Bairdi+RA2 Flathead sole")
plot(r47~year,pch=16 , main = "Reproductive female Bairdi+NBT")
plot(r48~year,pch=16 , main = "Reproductive female Bairdi+SE wind")
plot(r49~year,pch=16 , main = "Reproductive female Bairdi+AO")
plot(r50~year,pch=16 , main = "Reproductive female Bairdi+RA2 AO")
plot(r51~year,pch=16 , main = "Reproductive female Bairdi+RA3 AO")
plot(r52~year,pch=16 , main = "Reproductive female Bairdi+PDO")
plot(r53~year,pch=16 , main = "Reproductive female Bairdi+RA2 PDO")
plot(r54~year,pch=16 , main = "Reproductive female Bairdi+RA3 PDO")
plot(r55~year,pch=16 , main = "Reproductive female Bairdi+SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit40<-loess(r40~year,span=0.4)
r.fit41<-loess(r41~year.o,span=0.4)
r.fit42<-loess(r42~year,span=0.4)
r.fit43<-loess(r43~year,span=0.4)
r.fit44<-loess(r44~year,span=0.4)
r.fit45<-loess(r45~year.fhs,span=0.4)
r.fit46<-loess(r46~year.fhs2,span=0.4)
r.fit47<-loess(r47~year,span=0.4)
r.fit48<-loess(r48~year,span=0.4)
r.fit49<-loess(r49~year,span=0.4)
r.fit50<-loess(r50~year,span=0.4)
r.fit51<-loess(r51~year,span=0.4)
r.fit52<-loess(r52~year,span=0.4)
r.fit53<-loess(r53~year,span=0.4)
r.fit54<-loess(r54~year,span=0.4)
r.fit55<-loess(r55~year,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit40,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi")
lines(year,fitted(r.fit40),col=4)

plot(r.fit41,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Female opilio")
lines(year.o,fitted(r.fit41),col=4)

plot(r.fit42,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Pacific cod")
lines(year,fitted(r.fit42),col=4)

plot(r.fit43,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA2Pacific cod")
lines(year,fitted(r.fit43),col=4)

plot(r.fit44,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA3Pacific cod")
lines(year,fitted(r.fit44),col=4)

plot(r.fit45,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Flathead sole")
lines(year.fhs,fitted(r.fit45),col=4)

plot(r.fit46,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA2 Flathead sole")
lines(year.fhs2,fitted(r.fit46),col=4)

plot(r.fit47,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA3 NBT")
lines(year,fitted(r.fit47),col=4)

plot(r.fit48,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+SE wind")
lines(year,fitted(r.fit48),col=4)

plot(r.fit49,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Arctic Oscillation")
lines(year,fitted(r.fit49),col=4)

plot(r.fit50,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA2 Arctic Oscillation")
lines(year,fitted(r.fit50),col=4)

plot(r.fit51,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA3 Arctic Oscillation")
lines(year,fitted(r.fit51),col=4)

plot(r.fit52,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Pacific Decadal Oscillation")
lines(year,fitted(r.fit52),col=4)

plot(r.fit53,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA2 Pacific Decadal Oscillation")
lines(year,fitted(r.fit53),col=4)

plot(r.fit54,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+RA3 Pacific Decadal Oscillation")
lines(year,fitted(r.fit54),col=4)

plot(r.fit55,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+May-July SST")
lines(year,fitted(r.fit55),col=4)

##############################################################################################################################################
########################## Add juveniles 6 years prior to recruitment as predictor ###########################################################
###### Create data variables
n<-7
lnRS_6<-All_Dat$logRS[n:length(r)]
Fem_CO_6<-All_Dat$Ovig_female_CO[n:length(r)]
cod_6<-All_Dat$Pcod_lag1[n:length(r)]
cod2_6<-All_Dat$PCod_RA2[n:length(r)]
cod3_6<-All_Dat$PCod_RA3[n:length(r)]
fhs_6<-All_Dat$FHS_lag2[n:length(r)]
fhs2_6<-All_Dat$FHS_RA2[n:length(r)]
nbt_6<-All_Dat$NBT_3RA[n:length(r)]
wind_6<-All_Dat$SE.wind[n:length(r)]
ao_6<-All_Dat$AO_jfm[n:length(r)]
ao2_6<-All_Dat$AO_RA2[n:length(r)]
ao3_6<-All_Dat$AO_RA3[n:length(r)]
pdo_6<-All_Dat$PDO_djf[n:length(r)]
pdo2_6<-All_Dat$PDO_RA2[n:length(r)]
pdo3_6<-All_Dat$PDO_RA3[n:length(r)]
sst_6<-All_Dat$SST_May_July[n:length(r)]
year_6<-All_Dat$releaseyear[n:length(r)]


R6<-r[n:length(r)]
Rminus6<-r[1:(length(r)-6)]
length(Rminus6)
length(R6)
S6<-s[n:length(r)]
#logR7<-log(R7/S7)

acf(r)



mod_6_test <- gam(lnRS_6 ~ s(S6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")
summary(mod_6_test)

r_6<-resid(mod_6)
r.fit_6<-loess(r_6~year_6,span=0.4)

dev.new()
plot(r.fit_6,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Females and juveniles6 years prior")
lines(year_6,fitted(r.fit_6),col=4)
cor(R6,S6)
cor(R6,Rminus6)
cor(S6,Rminus6)

dwtest(mod_6_test)


########################################################################################################################################################################################
##################################### Run models adding juvenile abundance 6 years previously to multifactor models ##################################################################
#########################################################################################################################################################################################################
####################################### Now use GAM models with log(R/S) as response and with additional environmental variable ##################################################################################

mod60 <- gam(lnRS_6 ~ s(S6, k = 4)+s(Rminus6, k = 4),
           data = dat, method = "REML")

mod61 <- gam(lnRS_6 ~s(S6, k = 4)+ s(Fem_CO_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod62 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(cod_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod63 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(cod2_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod64 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(cod3_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod65 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(fhs_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod66 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(fhs2_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod67 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(nbt_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod68 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(wind_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod69 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(ao_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod70 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(ao2_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod71 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(ao3_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod72 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(pdo_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod73 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(pdo2_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod74 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(pdo3_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

mod75 <- gam(lnRS_6 ~ s(S6, k = 4)+ s(sst_6, k = 4)+s(Rminus6, k = 4),
             data = dat, method = "REML")

dwtest(mod60)
dwtest(mod61)
dwtest(mod62)
dwtest(mod63)
dwtest(mod64)
dwtest(mod65)
dwtest(mod66)
dwtest(mod67)
dwtest(mod68)
dwtest(mod69)
dwtest(mod70)
dwtest(mod71)
dwtest(mod72)
dwtest(mod73)
dwtest(mod74)
dwtest(mod75)
###############################################################################################################################################################################################
############################################### Model Summary ############################################################################

summary(mod60)
summary(mod61)
summary(mod62)
summary(mod63)
summary(mod64)
summary(mod65)
summary(mod66)
summary(mod67)
summary(mod68)
summary(mod69)
summary(mod70)
summary(mod71)
summary(mod72)
summary(mod73)
summary(mod74)
summary(mod75)
####################################################################################################################################################################################################
################################################ Plot GAM models ###################################################################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod60 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females with Juveniles 6 yrs prior")
plot(mod61, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs priorand Opilio females")
plot(mod62, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and Pacific cod")

dev.new()
par(mfrow=c(3,3))

plot(mod63, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA2 Pacific cod")
plot(mod64, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA3 Pacific cod")
plot(mod65, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and Flathead sole")

dev.new()
par(mfrow=c(3,3))

plot(mod66, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA2 Flathead sole")
plot(mod67, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and NBT")
plot(mod68, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and SE wind")

dev.new()
par(mfrow=c(3,3))

plot(mod69, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and AO")
plot(mod70, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA2 AO")
plot(mod71, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA3 AO")

# each of these models except first will have two plots (one for each variable), so need to set up second graphics window
dev.new()
par(mfrow=c(3,3))

plot(mod72, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and PDO")
plot(mod73, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA2 PDO")
plot(mod74, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and RA3 PDO")

dev.new()
par(mfrow=c(3,3))

plot(mod75, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 6 yrs prior and SST")

r60<-resid(mod60)
r61<-resid(mod61)
r62<-resid(mod62)
r63<-resid(mod63)
r64<-resid(mod64)
r65<-resid(mod65)
r66<-resid(mod66)
r67<-resid(mod67)
r68<-resid(mod68)
r69<-resid(mod69)
r70<-resid(mod70)
r71<-resid(mod71)
r72<-resid(mod72)
r73<-resid(mod73)
r74<-resid(mod74)
r75<-resid(mod75)

acf(r60)

dev.new()
par(mfrow=c(3,3))



plot(r60~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior")
plot(r61~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Reproductive female Opilio")
plot(r62~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Pacific cod")
plot(r63~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 Pacific cod")
plot(r64~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 Pacific cod")
plot(r65~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Flathead sole")
plot(r66~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 Flathead sole")
plot(r67~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and NBT")
plot(r68~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and SE wind")

par(mfrow=c(3,3))

plot(r69~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and AO")
plot(r70~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 AO")
plot(r71~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 AO")
plot(r72~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and PDO")
plot(r73~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 PDO")
plot(r74~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 PDO")
plot(r75~year_6,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit60<-loess(r60~year_6,span=0.4)
r.fit61<-loess(r61~year_6,span=0.4)
r.fit62<-loess(r62~year_6,span=0.4)
r.fit63<-loess(r63~year_6,span=0.4)
r.fit64<-loess(r64~year_6,span=0.4)
r.fit65<-loess(r65~year_6,span=0.4)
r.fit66<-loess(r66~year_6,span=0.4)
r.fit67<-loess(r67~year_6,span=0.4)
r.fit68<-loess(r68~year_6,span=0.4)
r.fit69<-loess(r69~year_6,span=0.4)
r.fit70<-loess(r70~year_6,span=0.4)
r.fit71<-loess(r71~year_6,span=0.4)
r.fit72<-loess(r72~year_6,span=0.4)
r.fit73<-loess(r73~year_6,span=0.4)
r.fit74<-loess(r74~year_6,span=0.4)
r.fit75<-loess(r75~year_6,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit60,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi +Juveniles 6 yrs prior")
lines(year_6,fitted(r.fit60),col=4)

plot(r.fit61,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and Female opilio")
lines(year_6,fitted(r.fit61),col=4)

plot(r.fit62,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and Pacific cod")
lines(year_6,fitted(r.fit62),col=4)

plot(r.fit63,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA2Pacific cod")
lines(year_6,fitted(r.fit63),col=4)

plot(r.fit64,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA3Pacific cod")
lines(year_6,fitted(r.fit64),col=4)

plot(r.fit65,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and Flathead sole")
lines(year_6,fitted(r.fit65),col=4)

plot(r.fit66,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA2 Flathead sole")
lines(year_6,fitted(r.fit66),col=4)

plot(r.fit67,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA3 NBT")
lines(year_6,fitted(r.fit67),col=4)

plot(r.fit68,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and SE wind")
lines(year_6,fitted(r.fit68),col=4)

dev.new()
par(mfrow=c(3,3))

plot(r.fit69,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and Arctic Oscillation")
lines(year_6,fitted(r.fit69),col=4)

plot(r.fit70,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA2 Arctic Oscillation")
lines(year_6,fitted(r.fit70),col=4)

plot(r.fit71,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA3 Arctic Oscillation")
lines(year_6,fitted(r.fit71),col=4)

plot(r.fit72,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and Pacific Decadal Oscillation")
lines(year_6,fitted(r.fit72),col=4)

plot(r.fit73,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA2 Pacific Decadal Oscillation")
lines(year_6,fitted(r.fit73),col=4)

plot(r.fit74,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and RA3 Pacific Decadal Oscillation")
lines(year_6,fitted(r.fit74),col=4)

plot(r.fit75,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 6 yrs prior and May-July SST")
lines(year_6,fitted(r.fit75),col=4)


##############################################################################################################################################
########################## Add juveniles 7 years prior to recruitment as predictor ###########################################################
###### Create data variables
n<-8
lnRS_7<-All_Dat$logRS[n:length(r)]
Fem_CO_7<-All_Dat$Ovig_female_CO[n:length(r)]
cod_7<-All_Dat$Pcod_lag1[n:length(r)]
cod2_7<-All_Dat$PCod_RA2[n:length(r)]
cod3_7<-All_Dat$PCod_RA3[n:length(r)]
fhs_7<-All_Dat$FHS_lag2[n:length(r)]
fhs2_7<-All_Dat$FHS_RA2[n:length(r)]
nbt_7<-All_Dat$NBT_3RA[n:length(r)]
wind_7<-All_Dat$SE.wind[n:length(r)]
ao_7<-All_Dat$AO_jfm[n:length(r)]
ao2_7<-All_Dat$AO_RA2[n:length(r)]
ao3_7<-All_Dat$AO_RA3[n:length(r)]
pdo_7<-All_Dat$PDO_djf[n:length(r)]
pdo2_7<-All_Dat$PDO_RA2[n:length(r)]
pdo3_7<-All_Dat$PDO_RA3[n:length(r)]
sst_7<-All_Dat$SST_May_July[n:length(r)]
year_7<-All_Dat$releaseyear[n:length(r)]


R7<-r[n:length(r)]
Rminus7<-r[1:(length(r)-7)]
length(Rminus7)
length(R7)
S7<-s[n:length(r)]
#logR7<-log(R7/S7)

acf(r)



mod_7_test <- gam(lnRS_7 ~ s(S7, k = 4)+s(Rminus7, k = 4),
                  data = dat, method = "REML")
summary(mod_7_test)

r_7<-resid(mod_7)
r.fit_7<-loess(r_7~year_7,span=0.4)

dev.new()
plot(r.fit_7,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Females and juveniles 7 years prior")
lines(year_7,fitted(r.fit_7),col=4)
cor(R7,S7)
cor(R6,S6)
cor(R7,Rminus7)
cor(S7,Rminus7)
cor(S6,Rminus6)

dwtest(mod_7_test)


##############################################################################################################################
##################################### Run models adding juvenile abundance 7 years previously to multifactor models ##########
#########################################################################################################################################################################################################
####################################### Now use GAM models with log(R/S) as response and with additional environmental variable ##################################################################################

mod80 <- gam(lnRS_7 ~ s(S7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod81 <- gam(lnRS_7 ~s(S7, k = 4)+ s(Fem_CO_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod82 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(cod_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod83 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(cod2_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod84 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(cod3_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod85 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(fhs_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod86 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(fhs2_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod87 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(nbt_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod88 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(wind_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod89 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(ao_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod90 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(ao2_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod91 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(ao3_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod92 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(pdo_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod93 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(pdo2_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod94 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(pdo3_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

mod95 <- gam(lnRS_7 ~ s(S7, k = 4)+ s(sst_7, k = 4)+s(Rminus7, k = 4),
             data = dat, method = "REML")

dwtest(mod80)
dwtest(mod81)
dwtest(mod82)
dwtest(mod83)
dwtest(mod84)
dwtest(mod85)
dwtest(mod86)
dwtest(mod87)
dwtest(mod88)
dwtest(mod89)
dwtest(mod90)
dwtest(mod91)
dwtest(mod92)
dwtest(mod93)
dwtest(mod94)
dwtest(mod95)


summary(mod80)
summary(mod81)
summary(mod82)
summary(mod83)
summary(mod84)
summary(mod85)
summary(mod86)
summary(mod87)
summary(mod88)
summary(mod89)
summary(mod90)
summary(mod91)
summary(mod92)
summary(mod93)
summary(mod94)
summary(mod95)
  `####################################################################################################################################################################################################
################################################ Plot GAM models ###################################################################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod80 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females with Juveniles 7 yrs prior")
plot(mod81, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs priorand Opilio females")
plot(mod82, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and Pacific cod")

dev.new()
par(mfrow=c(3,3))

plot(mod83, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA2 Pacific cod")
plot(mod84, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA3 Pacific cod")
plot(mod85, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and Flathead sole")

dev.new()
par(mfrow=c(3,3))

plot(mod86, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA2 Flathead sole")
plot(mod87, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and NBT")
plot(mod88, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and SE wind")

dev.new()
par(mfrow=c(3,3))

plot(mod89, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and AO")
plot(mod90, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA2 AO")
plot(mod91, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA3 AO")

# each of these models except first will have two plots (one for each variable), so need to set up second graphics window
dev.new()
par(mfrow=c(3,3))

plot(mod82, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and PDO")
plot(mod83, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA2 PDO")
plot(mod84, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and RA3 PDO")

dev.new()
par(mfrow=c(3,3))

plot(mod85, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi with Juveniles 7 yrs prior and SST")

r80<-resid(mod80)
r81<-resid(mod81)
r82<-resid(mod82)
r83<-resid(mod83)
r84<-resid(mod84)
r85<-resid(mod85)
r86<-resid(mod86)
r87<-resid(mod87)
r88<-resid(mod88)
r89<-resid(mod89)
r90<-resid(mod90)
r91<-resid(mod91)
r92<-resid(mod92)
r93<-resid(mod93)
r94<-resid(mod94)
r95<-resid(mod95)

acf(r80)

dev.new()
par(mfrow=c(3,3))



plot(r80~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior")
plot(r81~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Reproductive female Opilio")
plot(r82~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Pacific cod")
plot(r83~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 Pacific cod")
plot(r84~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 Pacific cod")
plot(r85~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and Flathead sole")
plot(r86~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 Flathead sole")
plot(r87~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and NBT")
plot(r88~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and SE wind")

dev.new()
par(mfrow=c(3,3))

plot(r89~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and AO")
plot(r90~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 AO")
plot(r91~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 AO")
plot(r92~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and PDO")
plot(r93~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA2 PDO")
plot(r94~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and RA3 PDO")
plot(r95~year_7,pch=16 , main = "Reproductive female Bairdi+Juveniles 6 yrs prior and SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit80<-loess(r80~year_7,span=0.4)
r.fit81<-loess(r81~year_7,span=0.4)
r.fit82<-loess(r82~year_7,span=0.4)
r.fit83<-loess(r83~year_7,span=0.4)
r.fit84<-loess(r84~year_7,span=0.4)
r.fit85<-loess(r85~year_7,span=0.4)
r.fit86<-loess(r86~year_7,span=0.4)
r.fit87<-loess(r87~year_7,span=0.4)
r.fit88<-loess(r88~year_7,span=0.4)
r.fit89<-loess(r89~year_7,span=0.4)
r.fit90<-loess(r90~year_7,span=0.4)
r.fit91<-loess(r91~year_7,span=0.4)
r.fit92<-loess(r92~year_7,span=0.4)
r.fit93<-loess(r93~year_7,span=0.4)
r.fit94<-loess(r94~year_7,span=0.4)
r.fit95<-loess(r95~year_7,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit80,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior")
lines(year_7,fitted(r.fit80),col=4)

plot(r.fit81,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and Female opilio")
lines(year_7,fitted(r.fit81),col=4)

plot(r.fit82,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and Pacific cod")
lines(year_7,fitted(r.fit82),col=4)

plot(r.fit83,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA2Pacific cod")
lines(year_7,fitted(r.fit83),col=4)

plot(r.fit84,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA3Pacific cod")
lines(year_7,fitted(r.fit84),col=4)

plot(r.fit85,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and Flathead sole")
lines(year_7,fitted(r.fit85),col=4)

plot(r.fit86,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA2 Flathead sole")
lines(year_7,fitted(r.fit86),col=4)

plot(r.fit87,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA3 NBT")
lines(year_7,fitted(r.fit87),col=4)

plot(r.fit88,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and SE wind")
lines(year_7,fitted(r.fit88),col=4)

dev.new()
par(mfrow=c(3,3))

plot(r.fit89,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and Arctic Oscillation")
lines(year_7,fitted(r.fit89),col=4)

plot(r.fit90,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA2 Arctic Oscillation")
lines(year_7,fitted(r.fit90),col=4)

plot(r.fit91,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA3 Arctic Oscillation")
lines(year_7,fitted(r.fit91),col=4)

plot(r.fit92,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and Pacific Decadal Oscillation")
lines(year_7,fitted(r.fit92),col=4)

plot(r.fit93,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA2 Pacific Decadal Oscillation")
lines(year_7,fitted(r.fit93),col=4)

plot(r.fit94,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and RA3 Pacific Decadal Oscillation")
lines(year_7,fitted(r.fit94),col=4)

plot(r.fit95,pch=16,xlab="Hatch year", ylab="Ln(R/S) residual",main="Female Bairdi+Juveniles 7 yrs prior and May-July SST")
lines(year_7,fitted(r.fit95),col=4)