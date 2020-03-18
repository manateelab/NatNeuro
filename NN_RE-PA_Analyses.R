
A<-Sys.getenv("LOGNAME")
setwd(paste("/Users/", A, "/Dropbox/Cornell/Research/WCMCseed/Analysis/GPX_AnalysisReneau",sep=""))
movesEMA = read.table("moves_daily_masterfile_v18_NDVI_Dimension_entropy.csv", header=T, sep=",",as.is=TRUE,strip.white=TRUE,fill=TRUE)
movesEMA$RE <- movesEMA$RE/log(648000000)
movesEMA$RE_NoInterpolation <- movesEMA$RE_nointerpolation/log(648000000)

moves2a <- subset(movesEMA, select = c(subject.id,
                                        date,
                                        cohort,
                                        mean.temp,
                                        precipi,
                                        RE,
                                        minute.ema.was.collected,
                                        PosMean,
                                        NegMean,
                                        distances,
                                        unique_locations,
                                        novel_locations,
                                        RE_two_decimal,
                                        UL_two_decimal,
                                        TOP3_PCA3_HL,
                                        RE_NoInterpolation))
write.csv(moves2a, file="/Users/aheller/Desktop/moves_daily_masterfile_share2.csv")
rm(list=ls())
library(afex)
library(lmerTest)
library(dplyr)
library(nlme)
A<-Sys.getenv("LOGNAME")
setwd(paste("/Users/", A, "/Dropbox/Miami/Research/MovesEMA/",sep=""))
movesEMA = read.table("moves_daily_masterfile_share2.csv", header=T, sep=",",as.is=TRUE,strip.white=TRUE,fill=TRUE)
movesEMA$DOW <- as.factor(weekdays(as.Date(movesEMA$date, "%m/%d/%y"))) #Convert date to day of week
movesEMA$cohort <- as.factor(movesEMA$cohort) # Convert cohort number to a factor
movesEMA$subject.id <- as.factor(movesEMA$subject.id) # Convert cohort number to a factor
movesEMA <- movesEMA %>% group_by(subject.id) %>% mutate(subj_RowNum = dplyr::row_number())

movesEMA <- data.frame(movesEMA)
movesEMA2 <- movesEMA
for (ch in levels(movesEMA$subject.id)) {
  currSub <- subset(movesEMA, subject.id==ch)
  if (sum(!is.na(currSub$PosMean))<20) {
    movesEMA2 <- movesEMA2[movesEMA2$subject.id != ch, ]
  }
}

movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(RE_z = scale(RE))

##### Positive Affect - RE association #####
zRE_M6 <- lme(PosMean ~ RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2, random =~1 + RE_z|subject.id, na.action = na.omit, correlation=corAR1())
summary(zRE_M6)
randeff2<-coef(zRE_M6)
randeff2$subID <- rownames(randeff2)
randeff2 <- subset(randeff2, select = c(subID, RE_z))

##### Negative Affect - RE association #####
zRE_M6.NA <- lme(NegMean ~ RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2, random =~1 + RE_z|subject.id, na.action = na.omit, correlation=corAR1())
summary(zRE_M6.NA)

##### Novel Locations x RE #####
xx<-movesEMA2 %>% group_by(subject.id) %>% summarize(COR=cor(RE,novel_locations, use = "pairwise.complete.obs"))
mean(xx$COR)
novelLocRE <- lme(RE ~ novel_locations, data=movesEMA2, random = ~1 | subject.id , na.action = na.omit, correlation=corAR1())
summary(novelLocRE)

##### Positive Affect - Novel Locations Association  #####
movesEMA3 <- movesEMA2[movesEMA2$subj_RowNum > 10,]
movesEMA3 <- movesEMA3 %>% group_by(subject.id) %>% mutate(novel_locations_z = scale(novel_locations))

novelLocPA <- lme(PosMean ~ novel_locations_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA3,
                  random = ~1 + novel_locations_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(novelLocPA)

novelLocPA_RE <- lme(PosMean ~ novel_locations_z + RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA3,
                     random = ~1 + novel_locations_z + RE_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(novelLocPA_RE)

##### Positive Affect - Sociodemographic space exploration #####
RE.PCA <- lme(TOP3_PCA3_HL ~ RE, data=movesEMA2, random = ~1 | subject.id, na.action = na.omit, correlation = corAR1())
summary(RE.PCA)

movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(TOP3_PCA3_HL_z = scale(TOP3_PCA3_HL))
TOP3_PCA3_HL_M6 <- lme(PosMean ~ TOP3_PCA3_HL_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2,
                       random = ~1 + TOP3_PCA3_HL_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(TOP3_PCA3_HL_M6)

TOP3_PCA3_HL_M6 <- lme(PosMean ~ TOP3_PCA3_HL_z + RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2,
                       random = ~1 + TOP3_PCA3_HL_z + RE_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(TOP3_PCA3_HL_M6)

##### Directionality analyses for RE ####
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>%  mutate(lag1.RE_z = dplyr::lag(RE_z, n = 1, default = NA))
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>%  mutate(lead1.RE_z = dplyr::lead(RE_z, n = 1, default = NA))

# Lag analysis
zRE_M6.lag <- lme(PosMean ~ RE_z + lag1.RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2, 
                  random = ~1 + RE_z + lag1.RE_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(zRE_M6.lag)

# Lead analysis #note - model won't fit with RE as a random slope #
zRE_M6.lead <- lme(lead1.RE_z ~ RE_z + PosMean + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2, 
                  random = ~1 + PosMean | subject.id, na.action = na.omit, correlation = corAR1())
summary(zRE_M6.lead)


##### Directionality analyses for Novelty #####
movesEMA3 <- movesEMA3 %>% group_by(subject.id) %>%  mutate(lag1.novel_locations_z = dplyr::lag(novel_locations_z, n = 1, default = NA))
movesEMA3 <- movesEMA3 %>% group_by(subject.id) %>%  mutate(lead1.novel_locations_z = dplyr::lead(novel_locations_z, n = 1, default = NA))

# Lag analysis
znovelLoc_M6.lag <- lme(PosMean ~ novel_locations_z + lag1.novel_locations_z  + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort,
                        data=movesEMA3, random = ~1 + novel_locations_z + lag1.novel_locations_z | subject.id, 
                        na.action = na.omit, correlation = corAR1())
summary(znovelLoc_M6.lag) 

# Lead analysis
znovelLoc_M6.lead <- lme(lead1.novel_locations_z ~ PosMean + novel_locations_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort,
                         data=movesEMA3, random = ~1 + PosMean  | subject.id, 
                         na.action = na.omit, correlation = corAR1())

summary(znovelLoc_M6.lead) 

##### Directionality analyses for Sociodemographic space #####
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>%  mutate(lag1.TOP3_PCA3_HL_z = dplyr::lag(TOP3_PCA3_HL_z, n = 1, default = NA))
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>%  mutate(lead1.TOP3_PCA3_HL_z = dplyr::lead(TOP3_PCA3_HL_z, n = 1, default = NA))

#Lag analysis
PCA.HL_M6.lag <- lme(PosMean ~  TOP3_PCA3_HL_z + lag1.TOP3_PCA3_HL_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, 
                     data=movesEMA2, random = ~1 + TOP3_PCA3_HL_z | subject.id,
                     na.action = na.omit, correlation = corAR1())

summary(PCA.HL_M6.lag)

#Lead analysis
PCA.HL_M6.lead <- lme(lead1.TOP3_PCA3_HL_z ~ PosMean + TOP3_PCA3_HL_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort,
                      data=movesEMA2, random = ~1 + PosMean | subject.id,
                      na.action = na.omit, correlation = corAR1())

summary(PCA.HL_M6.lead) 


##### SUPPLEMENTAL Analyses #######
PA_NA.lme  <- lme(PosMean ~ NegMean, data = movesEMA2,
                     random = ~1 + NegMean | subject.id, na.action = na.omit, correlation = corAR1())
summary(PA_NA.lme) 

##### Unique Location Analyses #######
UniqueLoc_RE  <- lme(RE ~ unique_locations, data = movesEMA2,
                     random = ~1 | subject.id, na.action = na.omit, correlation = corAR1())
summary(UniqueLoc_RE) 

movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(unique_locations_z = scale(unique_locations))
zuniqueLoc_M6 <- lme(PosMean ~ unique_locations_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data = movesEMA2,
                     random = ~1 + unique_locations_z | subject.id,
                     na.action = na.omit, correlation = corAR1())
summary(zuniqueLoc_M6)

zuniqueLoc_M6_RE <- lme(PosMean ~ unique_locations_z + RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2,
                      random = ~1 + RE_z + unique_locations_z | subject.id,
                      na.action = na.omit, correlation = corAR1())
summary(zuniqueLoc_M6_RE) 

##### Lower GPS spatial resolution Analyses #####
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(RE_two_decimal_z = scale(RE_two_decimal))
RE2Decimal_M6 <- lme(PosMean ~  RE_two_decimal_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data = movesEMA2,
                     random = ~1 + RE_two_decimal_z | subject.id,
                     na.action = na.omit, correlation = corAR1())
summary(RE2Decimal_M6)

movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(UL_two_decimal_z = scale(UL_two_decimal))
RE.UL.twoDecimal_M6 <- lme(PosMean ~ RE_two_decimal_z + UL_two_decimal_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA2,
                      random =~1 + RE_two_decimal_z + UL_two_decimal_z | subject.id, na.action = na.omit, correlation = corAR1())
summary(RE.UL.twoDecimal_M6)


##### RE - Positive affect analyses with no interpolation for location #####
movesEMA2 <- movesEMA2 %>% group_by(subject.id) %>% mutate(RE_NoInterpolation_z = scale(RE_NoInterpolation))
RE_NoInterpolation_M6 <- lme(PosMean ~  RE_NoInterpolation_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data = movesEMA2,
                     random = ~1 + RE_NoInterpolation_z | subject.id,
                     na.action = na.omit, correlation = corAR1())
summary(RE_NoInterpolation_M6)


##### Identifying best fitting model #####
#subset data frame to be able to "dredge it" - identify best fitting model for RE - Positive Affect association
library(MuMIn)

myvars<-c("PosMean",  "RE_z", "minute.ema.was.collected","mean.temp","precipi","distances", "DOW", "cohort","subject.id")
movesEMA4 <- movesEMA2[,myvars]
movesEMA4<-movesEMA4[complete.cases(movesEMA4),]

options(na.action = "na.fail")
zRE_M6 <- lme(PosMean ~ RE_z + minute.ema.was.collected + mean.temp + precipi + distances + DOW + cohort, data=movesEMA4, 
              random =~1 + RE_z|subject.id, correlation=corAR1())

A1<-dredge(zRE_M6)
options(na.action = "na.omit")

zRE.bestAIC <- lme(PosMean ~  RE_z + DOW + cohort, data=movesEMA2,
                   random = ~1 + RE_z | subject.id, na.action = na.omit, correlation=corAR1())
summary(zRE.bestAIC)

