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

All_Dat<-read.csv("data/EBS_Crab_and_envar_data_for_analysis.csv")
SR_Dat<-read.csv("data/SR_data.csv")

names(All_Dat)
names(SR_Dat)
###### Create data variables
r<-All_Dat$Lag3_recruits
s<-All_Dat$ReproductiveFemales
lnRS<-All_Dat$logRS
Fem_CO<-All_Dat$Ovig_female_CO
cod<-All_Dat$Pcod_lag1
fhs<-All_Dat$FHS_lag2
nbt<-All_Dat$NBT_3RA
wind<-All_Dat$SE.wind
ao<-All_Dat$AO_jfm
pdo<-All_Dat$PDO_djf
sst<-All_Dat$SST_May_July
year<-All_Dat$releaseyear


############################################################################################################################################################################
################################################ LINEAR MODELS #############################################################################################################

########################### Spawners ######################################################################
#dev.new()
#par(mfrow=c(2,2))
lm1<-lm(r~s)
summary(lm1)

#plot(lm1)

#dev.new()
#par(mfrow=c(1,1))
res<-resid(lm1)

#plot(res~year,pch=16 )

r.fit1<-loess(res~year,span=0.4)
#plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using female bairdi as predictor")
#lines(year,fitted(r.fit1),col=4)
########################### female opilio ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm2<-lm(r~Fem_CO)
summary(lm2)
res<-resid(lm2)
#plot(lm2)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit2<-loess(res~year,span=0.4)
#plot(r.fit2,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using female opilio as predictor")
#lines(year,fitted(r.fit2),col=4)
########################### Pacific cod######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm3<-lm(r~cod)
summary(lm3)
res<-resid(lm3)
#plot(lm3)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit3<-loess(res~year,span=0.4)
#plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using Pacific cod as predictor")
#lines(year,fitted(r.fit3),col=4)
########################### Flathead sole ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm4<-lm(r~fhs)
summary(lm4)
res<-resid(lm4)
#plot(lm4)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit4<-loess(res~year,span=0.4)
#plot(r.fit4,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using flathead sole as predictor")
#lines(year,fitted(r.fit4),col=4)
########################### NBT######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm5<-lm(r~nbt)
summary(lm5)
res<-resid(lm5)
#plot(lm5)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit5<-loess(res~year,span=0.4)
#plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using NBT as predictor")
#lines(year,fitted(r.fit5),col=4)
########################### wind ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm6<-lm(r~wind)
summary(lm6)
res<-resid(lm6)
#plot(lm6)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit6<-loess(res~year,span=0.4)
#plot(r.fit6,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using SE wind as predictor")
#lines(year,fitted(r.fit6),col=4)
########################### AO ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm7<-lm(r~ao)
summary(lm7)
res<-resid(lm7)
#plot(lm7)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit7<-loess(res~year,span=0.4)
#plot(r.fit7,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using AO as predictor")
#lines(year,fitted(r.fit7),col=4)
############################ PDO ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm8<-lm(r~pdo)
summary(lm8)
res<-resid(lm8)
#plot(lm8)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit8<-loess(res~year,span=0.4)
#plot(r.fit8,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using PDO as predictor")
#lines(year,fitted(r.fit8),col=4)
########################### SST ######################################################################
#dev.new()
#par(mfrow=c(2,2))

lm9<-lm(r~sst)
summary(lm9)
res<-resid(lm9)
#plot(lm9)

#dev.new()
#par(mfrow=c(1,1))
#plot(res~year,pch=16 )

r.fit9<-loess(res~year,span=0.4)
#plot(r.fit9,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Interdecadal trend in EBS group residuals using SST as predictor")
#lines(year,fitted(r.fit9),col=4)


################################################################################################################################################
#################################### plot all together #############################################################################################

dev.new()
par(mfrow=c(3,3))


plot(r.fit1,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female Bairdi")
lines(year,fitted(r.fit1),col=4)

plot(r.fit2,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female opilio")
lines(year,fitted(r.fit2),col=4)

plot(r.fit3,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod")
lines(year,fitted(r.fit3),col=4)

plot(r.fit4,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole")
lines(year,fitted(r.fit4),col=4)

plot(r.fit5,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="NBT")
lines(year,fitted(r.fit5),col=4)

plot(r.fit6,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="SE wind")
lines(year,fitted(r.fit6),col=4)

plot(r.fit7,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Arctic Oscillation")
lines(year,fitted(r.fit7),col=4)

plot(r.fit8,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit8),col=4)

plot(r.fit9,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="May-July SST")
lines(year,fitted(r.fit9),col=4)

##############################################################################################################################################
####################################### Now use GAM models ##################################################################################

mod1 <- gam(r ~ s(s, k = 4),
            data = dat, family="Gamma",method = "REML")

mod2 <- gam(r ~ s(Fem_CO, k = 4),
            data = dat, family="Gamma",method = "REML")

mod3 <- gam(r ~ s(cod, k = 4),
            data = dat, family="Gamma",method = "REML")

mod4 <- gam(r ~ s(fhs, k = 4),
            data = dat, family="Gamma",method = "REML")

mod5 <- gam(r ~ s(nbt, k = 4),
            data = dat, family="Gamma",method = "REML")

mod6 <- gam(r ~ s(wind, k = 4),
            data = dat, family="Gamma",method = "REML")

mod7 <- gam(r ~ s(ao, k = 4),
            data = dat, family="Gamma",method = "REML")

mod8 <- gam(r ~ s(pdo, k = 4),
            data = dat, family="Gamma",method = "REML")

mod9 <- gam(r ~ s(sst, k = 4),
            data = dat, family="Gamma",method = "REML")

######################################################################################################################################################################
################################################ Plot GAM models ############################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod1 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod2, resid=T, pch=16, shade=T, rug=F, main = "Opilio females")
plot(mod3, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod")
plot(mod4, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole")
plot(mod5, resid=T, pch=16, shade=T, rug=F, main = "NBT")
plot(mod6, resid=T, pch=16, shade=T, rug=F, main = "SE wind")
plot(mod7, resid=T, pch=16, shade=T, rug=F, main = "AO")
plot(mod8, resid=T, pch=16, shade=T, rug=F, main = "PDO")
plot(mod9, resid=T, pch=16, shade=T, rug=F, main = "SST")

r1<-resid(mod1)
r2<-resid(mod2)
r3<-resid(mod3)
r4<-resid(mod4)
r5<-resid(mod5)
r6<-resid(mod6)
r7<-resid(mod7)
r8<-resid(mod8)
r9<-resid(mod9)

dev.new()
par(mfrow=c(3,3))

plot(r1~year,pch=16 , main = "Reproductive female Bairdi")
plot(r2~year,pch=16 , main = "Reproductive female Opilio")
plot(r3~year,pch=16 , main = "Pacific cod")
plot(r4~year,pch=16 , main = "Flathead sole")
plot(r5~year,pch=16 , main = "NBT")
plot(r6~year,pch=16 , main = "SE wind")
plot(r7~year,pch=16 , main = "AO")
plot(r8~year,pch=16 , main = "PDO")
plot(r9~year,pch=16 , main = "SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit10<-loess(r1~year,span=0.4)
r.fit11<-loess(r2~year,span=0.4)
r.fit12<-loess(r3~year,span=0.4)
r.fit13<-loess(r4~year,span=0.4)
r.fit14<-loess(r5~year,span=0.4)
r.fit15<-loess(r6~year,span=0.4)
r.fit16<-loess(r7~year,span=0.4)
r.fit17<-loess(r8~year,span=0.4)
r.fit18<-loess(r9~year,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit10,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female Bairdi")
lines(year,fitted(r.fit10),col=4)

plot(r.fit11,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female opilio")
lines(year,fitted(r.fit11),col=4)

plot(r.fit12,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod")
lines(year,fitted(r.fit12),col=4)

plot(r.fit13,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole")
lines(year,fitted(r.fit13),col=4)

plot(r.fit14,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="NBT")
lines(year,fitted(r.fit14),col=4)

plot(r.fit15,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="SE wind")
lines(year,fitted(r.fit15),col=4)

plot(r.fit16,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Arctic Oscillation")
lines(year,fitted(r.fit16),col=4)

plot(r.fit17,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit17),col=4)

plot(r.fit18,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="May-July SST")
lines(year,fitted(r.fit18),col=4)

###########################################################################################################################################
########################################## GLS Models #####################################################################################

fit21<-gls(r ~ s, correlation = corAR1(), method="ML")
fit22<-gls(r ~ Fem_CO, correlation = corAR1(), method="ML")
fit23<-gls(r ~ cod, correlation = corAR1(), method="ML")
fit24<-gls(r ~ fhs, correlation = corAR1(), method="ML")
fit25<-gls(r ~ nbt, correlation = corAR1(), method="ML")
fit26<-gls(r ~ wind, correlation = corAR1(), method="ML")
fit27<-gls(r ~ ao, correlation = corAR1(), method="ML")
fit28<-gls(r ~ pdo, correlation = corAR1(), method="ML")
fit29<-gls(r ~ sst, correlation = corAR1(), method="ML")

r21<-resid(fit21)
r22<-resid(fit22)
r23<-resid(fit23)
r24<-resid(fit24)
r25<-resid(fit25)
r26<-resid(fit26)
r27<-resid(fit27)
r28<-resid(fit28)
r29<-resid(fit29)

r.fit21<-loess(r21~year,span=0.4)
r.fit22<-loess(r22~year,span=0.4)
r.fit23<-loess(r23~year,span=0.4)
r.fit24<-loess(r24~year,span=0.4)
r.fit25<-loess(r25~year,span=0.4)
r.fit26<-loess(r26~year,span=0.4)
r.fit27<-loess(r27~year,span=0.4)
r.fit28<-loess(r28~year,span=0.4)
r.fit29<-loess(r29~year,span=0.4)

dev.new()
par(mfrow=c(3,3))

plot(r.fit21,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female Bairdi")
lines(year,fitted(r.fit21),col=4)

plot(r.fit22,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female opilio")
lines(year,fitted(r.fit22),col=4)

plot(r.fit23,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod")
lines(year,fitted(r.fit23),col=4)

plot(r.fit24,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole")
lines(year,fitted(r.fit24),col=4)

plot(r.fit25,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="NBT")
lines(year,fitted(r.fit25),col=4)

plot(r.fit26,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="SE wind")
lines(year,fitted(r.fit26),col=4)

plot(r.fit27,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Arctic Oscillation")
lines(year,fitted(r.fit27),col=4)

plot(r.fit28,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit28),col=4)

plot(r.fit29,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="May-July SST")
lines(year,fitted(r.fit29),col=4)
