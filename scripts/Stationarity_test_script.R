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

##############################################################################################################################################
####################################### Now use GAM models with log(R/S) as response##################################################################################

mod41 <- gam(lnRS ~ s(s, k = 4),
            data = dat,method = "REML")

mod42 <- gam(lnRS ~ s(Fem_CO, k = 4),
            data = dat, method = "REML")

mod43 <- gam(lnRS ~ s(cod, k = 4),
            data = dat, method = "REML")

mod44 <- gam(lnRS ~ s(fhs, k = 4),
            data = dat, method = "REML")

mod45 <- gam(lnRS ~ s(nbt, k = 4),
            data = dat, method = "REML")

mod46 <- gam(lnRS ~ s(wind, k = 4),
            data = dat, method = "REML")

mod47 <- gam(lnRS ~ s(ao, k = 4),
            data = dat, method = "REML")

mod48 <- gam(lnRS ~ s(pdo, k = 4),
            data = dat, method = "REML")

mod49 <- gam(lnRS ~ s(sst, k = 4),
            data = dat, method = "REML")

######################################################################################################################################################################
################################################ Plot GAM models ############################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod41 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod42, resid=T, pch=16, shade=T, rug=F, main = "Opilio females")
plot(mod43, resid=T, pch=16, shade=T, rug=F, main = "Pacific cod")
plot(mod44, resid=T, pch=16, shade=T, rug=F, main = "Flathead sole")
plot(mod45, resid=T, pch=16, shade=T, rug=F, main = "NBT")
plot(mod46, resid=T, pch=16, shade=T, rug=F, main = "SE wind")
plot(mod47, resid=T, pch=16, shade=T, rug=F, main = "AO")
plot(mod48, resid=T, pch=16, shade=T, rug=F, main = "PDO")
plot(mod49, resid=T, pch=16, shade=T, rug=F, main = "SST")

r41<-resid(mod41)
r42<-resid(mod42)
r43<-resid(mod43)
r44<-resid(mod44)
r45<-resid(mod45)
r46<-resid(mod46)
r47<-resid(mod47)
r48<-resid(mod48)
r49<-resid(mod49)

dev.new()
par(mfrow=c(3,3))

plot(r41~year,pch=16 , main = "Reproductive female Bairdi")
plot(r42~year,pch=16 , main = "Reproductive female Opilio")
plot(r43~year,pch=16 , main = "Pacific cod")
plot(r44~year,pch=16 , main = "Flathead sole")
plot(r45~year,pch=16 , main = "NBT")
plot(r46~year,pch=16 , main = "SE wind")
plot(r47~year,pch=16 , main = "AO")
plot(r48~year,pch=16 , main = "PDO")
plot(r49~year,pch=16 , main = "SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit41<-loess(r41~year,span=0.4)
r.fit42<-loess(r42~year,span=0.4)
r.fit43<-loess(r43~year,span=0.4)
r.fit44<-loess(r44~year,span=0.4)
r.fit45<-loess(r45~year,span=0.4)
r.fit46<-loess(r46~year,span=0.4)
r.fit47<-loess(r47~year,span=0.4)
r.fit48<-loess(r48~year,span=0.4)
r.fit49<-loess(r49~year,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit41,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi")
lines(year,fitted(r.fit41),col=4)

plot(r.fit42,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female opilio")
lines(year,fitted(r.fit42),col=4)

plot(r.fit43,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific cod")
lines(year,fitted(r.fit43),col=4)

plot(r.fit44,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Flathead sole")
lines(year,fitted(r.fit44),col=4)

plot(r.fit45,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="NBT")
lines(year,fitted(r.fit45),col=4)

plot(r.fit46,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="SE wind")
lines(year,fitted(r.fit46),col=4)

plot(r.fit47,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Arctic Oscillation")
lines(year,fitted(r.fit47),col=4)

plot(r.fit48,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit48),col=4)

plot(r.fit49,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="May-July SST")
lines(year,fitted(r.fit49),col=4)

###############################################################################################################################################################
####################################### Now use GAM models with log(R/S) as response and with additional environmental variable ##################################################################################

mod51 <- gam(lnRS ~ s(s, k = 4),
             data = dat, method = "REML")

mod52 <- gam(lnRS ~s(s, k = 4)+ s(Fem_CO, k = 4),
             data = dat, method = "REML")

mod53 <- gam(lnRS ~ s(s, k = 4)+ s(cod, k = 4),
             data = dat, method = "REML")

mod54 <- gam(lnRS ~ s(s, k = 4)+ s(fhs, k = 4),
             data = dat, method = "REML")

mod55 <- gam(lnRS ~ s(s, k = 4)+ s(nbt, k = 4),
             data = dat, method = "REML")

mod56 <- gam(lnRS ~ s(s, k = 4)+ s(wind, k = 4),
             data = dat, method = "REML")

mod57 <- gam(lnRS ~ s(s, k = 4)+ s(ao, k = 4),
             data = dat, method = "REML")

mod58 <- gam(lnRS ~ s(s, k = 4)+ s(pdo, k = 4),
             data = dat, method = "REML")

mod59 <- gam(lnRS ~ s(s, k = 4)+ s(sst, k = 4),
             data = dat, method = "REML")

######################################################################################################################################################################
################################################ Plot GAM models ############################################################################################
dev.new()
par(mfrow=c(3,3))

plot(mod51 ,resid=T, pch=16, shade=T, rug=F, main = "Reproductive females")
plot(mod52, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Opilio females")
plot(mod53, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Pacific cod")
plot(mod54, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and Flathead sole")
plot(mod55, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and NBT")

# each of these models except first will have two plots (one for each variable), so need to set up second graphics window
dev.new()
par(mfrow=c(3,3))

plot(mod56, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and SE wind")
plot(mod57, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and AO")
plot(mod58, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and PDO")
plot(mod59, resid=T, pch=16, shade=T, rug=F, main = "Reproductive female bairdi and SST")

r51<-resid(mod51)
r52<-resid(mod52)
r53<-resid(mod53)
r54<-resid(mod54)
r55<-resid(mod55)
r56<-resid(mod56)
r57<-resid(mod57)
r58<-resid(mod58)
r59<-resid(mod59)

dev.new()
par(mfrow=c(3,3))

plot(r51~year,pch=16 , main = "Reproductive female Bairdi")
plot(r52~year,pch=16 , main = "Reproductive female Bairdi+Reproductive female Opilio")
plot(r53~year,pch=16 , main = "Reproductive female Bairdi+Pacific cod")
plot(r54~year,pch=16 , main = "Reproductive female Bairdi+Flathead sole")
plot(r55~year,pch=16 , main = "Reproductive female Bairdi+NBT")
plot(r56~year,pch=16 , main = "Reproductive female Bairdi+SE wind")
plot(r57~year,pch=16 , main = "Reproductive female Bairdi+AO")
plot(r58~year,pch=16 , main = "Reproductive female Bairdi+PDO")
plot(r59~year,pch=16 , main = "Reproductive female Bairdi+SST")





################################################################################################################################################
#################################### plot all together #############################################################################################
r.fit51<-loess(r51~year,span=0.4)
r.fit52<-loess(r52~year,span=0.4)
r.fit53<-loess(r53~year,span=0.4)
r.fit54<-loess(r54~year,span=0.4)
r.fit55<-loess(r55~year,span=0.4)
r.fit56<-loess(r56~year,span=0.4)
r.fit57<-loess(r57~year,span=0.4)
r.fit58<-loess(r58~year,span=0.4)
r.fit59<-loess(r59~year,span=0.4)

dev.new()
par(mfrow=c(3,3))


plot(r.fit51,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi")
lines(year,fitted(r.fit51),col=4)

plot(r.fit52,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Female opilio")
lines(year,fitted(r.fit52),col=4)

plot(r.fit53,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Pacific cod")
lines(year,fitted(r.fit53),col=4)

plot(r.fit54,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Flathead sole")
lines(year,fitted(r.fit54),col=4)

plot(r.fit55,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+NBT")
lines(year,fitted(r.fit55),col=4)

plot(r.fit56,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+SE wind")
lines(year,fitted(r.fit56),col=4)

plot(r.fit57,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Arctic Oscillation")
lines(year,fitted(r.fit57),col=4)

plot(r.fit58,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+Pacific Decadal Oscillation")
lines(year,fitted(r.fit58),col=4)

plot(r.fit59,pch=16,xlab="Hatch year", ylab="Ln(R/S)",main="Female Bairdi+May-July SST")
lines(year,fitted(r.fit59),col=4)

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


###########################################################################################################################################
########################################## GLS Models with spawners + single environmental variable  #####################################################################################

fit31<-gls(r ~ s, correlation = corAR1(), method="ML")
fit32<-gls(r ~ s+ Fem_CO, correlation = corAR1(), method="ML")
fit33<-gls(r ~ s+ cod, correlation = corAR1(), method="ML")
fit34<-gls(r ~ s+ fhs, correlation = corAR1(), method="ML")
fit35<-gls(r ~ s+ nbt, correlation = corAR1(), method="ML")
fit36<-gls(r ~ s+ wind, correlation = corAR1(), method="ML")
fit37<-gls(r ~ s+ ao, correlation = corAR1(), method="ML")
fit38<-gls(r ~ s+ pdo, correlation = corAR1(), method="ML")
fit39<-gls(r ~ s+ sst, correlation = corAR1(), method="ML")

r31<-resid(fit31)
r32<-resid(fit32)
r33<-resid(fit33)
r34<-resid(fit34)
r35<-resid(fit35)
r36<-resid(fit36)
r37<-resid(fit37)
r38<-resid(fit38)
r39<-resid(fit39)

r.fit31<-loess(r31~year,span=0.4)
r.fit32<-loess(r32~year,span=0.4)
r.fit33<-loess(r33~year,span=0.4)
r.fit34<-loess(r34~year,span=0.4)
r.fit35<-loess(r35~year,span=0.4)
r.fit36<-loess(r36~year,span=0.4)
r.fit37<-loess(r37~year,span=0.4)
r.fit38<-loess(r38~year,span=0.4)
r.fit39<-loess(r39~year,span=0.4)

dev.new()
par(mfrow=c(3,3))

plot(r.fit31,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female Bairdi")
lines(year,fitted(r.fit31),col=4)

plot(r.fit32,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Female opilio")
lines(year,fitted(r.fit32),col=4)

plot(r.fit33,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific cod")
lines(year,fitted(r.fit33),col=4)

plot(r.fit34,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Flathead sole")
lines(year,fitted(r.fit34),col=4)

plot(r.fit35,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="NBT")
lines(year,fitted(r.fit35),col=4)

plot(r.fit36,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="SE wind")
lines(year,fitted(r.fit36),col=4)

plot(r.fit37,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Arctic Oscillation")
lines(year,fitted(r.fit37),col=4)

plot(r.fit38,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="Pacific Decadal Oscillation")
lines(year,fitted(r.fit38),col=4)

plot(r.fit39,pch=16,xlab="Hatch year", ylab="Lag 3 residuals",main="May-July SST")
lines(year,fitted(r.fit39),col=4)
