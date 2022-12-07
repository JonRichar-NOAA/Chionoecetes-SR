dat0 <- read.csv("data/bairdi_strata_areas.csv") 
names(dat0)
plot(dat0$SURVEY_AREA~dat0$SURVEY_YEAR,pch=16)
