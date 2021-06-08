rm(list = ls())
library(JM)
library(table1)
library(ggplot2)
library(tidyverse)
library(comorbidity)
library(zoo)
library(rcompanion)
library(gridExtra)



set.seed(30)

dx <- read.csv("dx.csv")
elixhauser_icd10 <- comorbidity(x = as.data.frame(filter(dx,icd_version==10)), id = "hadm_id", 
                          code = "icd_code", score = "elixhauser", icd = "icd10", assign0 = FALSE)
elixhauser_icd9 <- comorbidity(x = as.data.frame(filter(dx,icd_version==9)), id = "hadm_id", 
                                code = "icd_code", score = "elixhauser", icd = "icd9", assign0 = FALSE)
elix = rbind(elixhauser_icd10,elixhauser_icd9)
EthinictyList = c('WHITE', 'ASIAN','BLACK/AFRICAN AMERICAN','HISPANIC/LATINO')



cohort <- read.csv("cohort.csv", stringsAsFactors=TRUE)%>%
  left_join(elix, by = c('hadm_id'))%>%
  mutate(male=ifelse(gender=="M",1,0))%>%
  mutate(ethnicity=as.factor(ethnicity))%>%
  mutate(ethnicity=relevel(ethnicity, ref = 'WHITE'))%>%
  mutate(sepsis=ifelse(sepsis3=="true",1,0))%>%  
  mutate(U_ICUHR = as.integer(ceiling(difftime(as.POSIXct(as.character(deathtime), format="%Y-%m-%d %H:%M:%OS"),
                                              as.POSIXct(as.character(intime), format="%Y-%m-%d %H:%M:%OS"), units="hours"))))%>%
  mutate(icu_day = U_ICUHR/24) #icu outime might be later than death time.
  

table1(~anchor_age + U_ICUHR + sofa_adm+ male+oasis+factor(sepsis)
       +factor(vent)+wscore_vw+factor(chf)+factor(rf)+factor(metacanc)+factor(aids)| ethnicity, 
       data=cohort,overall="Total")



sofa_recent <- read.csv("~/Desktop/icu_elos/sofa_hr.csv")%>%
  arrange(stay_id,hr)%>%
  mutate(renal_recent = ifelse((hr==0)&is.na(renal),0,renal))%>%
  mutate(renal_recent = na.locf(renal_recent,na.rm=FALSE))%>%
  mutate(liver_recent = ifelse((hr==0)&is.na(liver),0,liver))%>%
  mutate(liver_recent = na.locf(liver_recent,na.rm=FALSE))%>%
  mutate(cns_recent = ifelse((hr==0)&is.na(cns),0,cns))%>%
  mutate(cns_recent = na.locf(cns,na.rm=FALSE))%>%
  mutate(cv_recent = ifelse((hr==0)&is.na(cv),0,cv))%>%
  mutate(cv_recent = na.locf(cv_recent,na.rm=FALSE))%>%
  mutate(coag_recent = ifelse((hr==0)&is.na(coag),0,coag))%>%
  mutate(coag_recent = na.locf(coag_recent,na.rm=FALSE))%>%
  mutate(resp_recent = ifelse((hr==0)&is.na(resp),0,resp))%>%
  mutate(resp_recent = na.locf(resp_recent,na.rm=FALSE))%>%
  mutate(sofa_recent = renal_recent+liver_recent+cns_recent+cv_recent+coag_recent+resp_recent)%>%
  mutate(sofa_diff_recent = ifelse(stay_id==lag(stay_id),sofa_recent-lag(sofa_recent),999))%>%
  mutate(hr_day = hr/24)
  
  






  
  
dat = sofa_recent %>%
  filter(hr>=24)%>%
  left_join(cohort, by = c('stay_id'))%>%
  left_join(sofa_recent[sofa_recent$hr==24,c('stay_id','sofa_recent')], by = c('stay_id'))%>%
  filter(ethnicity=="WHITE" | ethnicity=="ASIAN" | ethnicity=="BLACK/AFRICAN AMERICAN" | ethnicity == "HISPANIC/LATINO")%>%
  mutate(ethnicity = factor(ethnicity))%>%
  filter(hr<=U_ICUHR)%>%
  mutate(sofa_diff = sofa_recent.x-sofa_recent.y)


dat.id = unique(dat[,c("stay_id","ethnicity","anchor_age","sofa_adm","wscore_vw","U_ICUHR",'male',"oasis","sepsis","vent")])
table1(~anchor_age + U_ICUHR + sofa_adm+ factor(male)+oasis+factor(sepsis)+factor(vent)+wscore_vw| ethnicity, 
       data=dat.id,overall="Total")

dat_dif = dat%>%
  #filter(U_ICUHR>=240)%>%
  filter(hr%%24==0)%>%
  mutate(hr2death = 24*(ceiling((hr-U_ICUHR)/24)))
  
dat_dif_fw = dat_dif%>%
  #filter(U_ICUHR>120)%>%
  filter(hr<=240)
  



Sum = groupwiseMean(sofa_recent.x ~ hr + ethnicity,
                    data   = dat_dif_fw,
                    conf   = 0.95,
                    digits = 3)
Sum$hr = as.factor(Sum$hr)
p1 = ggplot(Sum, aes(x=hr, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average SOFA score in 120 hrs after admission for people who survived more than 120 hrs")
Sum = groupwiseMean(sofa_diff ~ hr + ethnicity,
                    data   = dat_dif_fw,
                    conf   = 0.95,
                    digits = 3)
Sum$hr = as.factor(Sum$hr)
p2 = ggplot(Sum, aes(x=hr, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average relative change of SOFA score in 120 hrs after admission for people who survived more than 120 hrs")

grid.arrange(p1, p2, nrow = 2)

dat_dif = dat%>%
  filter(U_ICUHR<=120)%>%
  filter(hr%%24==0)%>%
  mutate(hr2death = 24*(ceiling((hr-U_ICUHR)/24)))


dat_dif_bw = dat_dif%>%
  left_join(dat_dif[dat_dif$hr2death==0,c('stay_id','sofa_recent.x')], by = c('stay_id'))%>%
  filter(!is.na(sofa_recent.x.y))%>%
  mutate(sofa_diff_death = sofa_recent.x.x-sofa_recent.x.y)%>%
  filter(hr2death>=-72)%>%
  filter(hr2death<=0)

Sum = groupwiseMean(sofa_recent.x.x ~ hr2death + ethnicity,
                    data   = dat_dif_bw,
                    conf   = 0.95,
                    digits = 3)
Sum$hr2death = as.factor(Sum$hr2death)
p1 = ggplot(Sum, aes(x=hr2death, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr2death, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average SOFA score for peole died in first 240 hr")


Sum = groupwiseMean(sofa_diff_death ~ hr2death + ethnicity,
                    data   = dat_dif_bw,
                    conf   = 0.95,
                    digits = 3)
Sum$hr2death = as.factor(Sum$hr2death)
p2 = ggplot(Sum, aes(x=hr2death, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr2death, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average relatve SOFA score change  for peole died in first 240 hr")

grid.arrange(p1, p2, nrow = 2)

mdat = dat  %>%
  mutate(age10=floor(anchor_age/10))%>%
  mutate(sofa_change = sofa_recent.x-sofa_recent.y)%>%
  mutate(died = ifelse(icu_day>10,0,1))%>%
  mutate(hr_day = hr_day-1)%>%
  mutate(icu_day = icu_day-0.999)%>%
  #filter(sofa_diff_recent>0)%>%
  filter(hr %% 6 ==0)%>%
  filter(hr<=240)

mdat = dat  %>%
  filter(icu_day>25/24)%>%
  mutate(age10=floor(anchor_age/10))%>%
  mutate(sofa_change = sofa_recent.x-sofa_recent.y)%>%
  filter(ifelse(hr>24*10,hr %% 24 == 0,sofa_diff_recent!=0))%>%
  mutate(hr_day = ifelse(hr_day<icu_day,hr_day+0.001,hr_day))%>%
  mutate(hr_day = hr_day-1)%>%
  mutate(died = 1)%>%
  mutate(icu_day = icu_day-0.999)


length(unique(mdat$stay_id))
mdat.id = unique(mdat[,c("stay_id","ethnicity","anchor_age","sofa_adm","wscore_vw","U_ICUHR",'male',"oasis",
                         "sepsis","vent","icu_day","died","age10","sofa_recent.y")])
table1(~anchor_age + U_ICUHR + sofa_adm+ factor(male)+oasis+factor(sepsis)+factor(vent)+wscore_vw| ethnicity, 
       data=mdat.id,overall="Total")

rawfit <- survfit(Surv(icu_day+1,died)~  ethnicity, data=mdat.id)
plot(rawfit,  lty=1:4,col=c(3,4,6,9),xlim=c(0, 10),
     xlab="Days After ICU admission", ylab="Survival")
legend(4, 1, EthinictyList,lty=1:4,col=c(3,4,6,9), bty='n')



fitLME_ns <- lme(sofa_change ~ bs(hr_day, knots = c(1,2,3,5,7))+ethnicity+sofa_recent.y+age10+male,
                 random = ~ hr_day | stay_id, 
                 data = mdat,control=list(opt = "optim"))


fitCox_ns <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+sofa_recent.y,data=mdat.id,x=TRUE)

fitJM_pwph_ns_1 <- jointModel(fitLME_ns,fitCox_ns,timeVar="hr_day",
                              method="piecewise-PH-aGH",
                              control=list(verbose=TRUE,iter.EM=10,optimizer='optim'))
summary(fitJM_pwph_ns_1)



fitJM_pwph_ns_2 <- jointModel(fitLME_ns,fitCox_ns,timeVar="hr_day",
                              method="piecewise-PH-aGH",interFact = list(value = ~ ethnicity, data = mdat.id),
                              control=list(verbose=TRUE,iter.EM=10,optimizer='optim'))
summary(fitJM_pwph_ns_2)



fitLME_ns <- lme(sofa_change~ bs(hr_day, knots = c(1,2,3,5,7,9))*ethnicity+sofa_recent.y+age10+male,
                 random = ~ hr_day | stay_id, 
                 data = mdat,control=list(opt = "optim"))
fitCox_ns <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+sofa_recent.y,data=mdat.id,x=TRUE)

fitJM_pwph_ns_3 <- jointModel(fitLME_ns,fitCox_ns,timeVar="hr_day",
                              method="piecewise-PH-aGH",
                              control=list(verbose=TRUE,iter.EM=2,optimizer='optim'))
summary(fitJM_pwph_ns_3)



fitJM_pwph_ns_4 <- jointModel(fitLME_ns,fitCox_ns,timeVar="hr_day",
                              method="piecewise-PH-aGH",interFact = list(value = ~ ethnicity, data = mdat.id),
                              control=list(verbose=TRUE,iter.EM=10,optimizer='optim'))
summary(fitJM_pwph_ns_4)





DF <- with(mdat, expand.grid(ethnicity = levels(ethnicity),
                             hr_day = seq(0, 9, len = 10)))
DF$sofa_recent.y=4
DF$age10=7
DF$male=1
mdat_pred <- predict(fitJM_pwph_ns_3, DF, interval = "confidence", return = TRUE)
require(lattice)
xyplot(pred + low + upp ~ hr_day | ethnicity, data = mdat_pred,
       type = "l", col = c(3,4,6,9), lty = c(1:4), lwd = 2,
       ylab = "Average Sofa Score")

ggplot(mdat_pred,
       aes(x = hr_day,
           y = pred)) +
  # Add a ribbon with the confidence band
  geom_smooth(
    aes(
      # lower and upper bound of the ribbon
      ymin = low, ymax = upp,
      # Different colour for men/women
      fill = ethnicity, colour = ethnicity
    ),
    stat = "identity") +
  xlab("ICU days") +
  ylab("Sofa")


ss <- survfitJM(fitJM_pwph_ns_1, newdata = mdat[mdat$stay_id=='33516489',], idVar = "stay_id", M = 50)
plot(ss,xlim=c(0,60))
plot(ss, include.y = TRUE, add.last.time.axis.tick = TRUE, legend = TRUE,conf.int=TRUE,fill.area=TRUE)









#=============code gradeyard===========================
#=============code gradeyard===========================
#=============code gradeyard===========================
#=============code gradeyard===========================
#=============code gradeyard===========================
#=============code gradeyard===========================
#=============code gradeyard===========================
fitLME_ns <- lme(sofa_recent.x ~ bs(hr_day, knots = c(log(4)), degree = 1)*ethnicity+sofa_recent.y+age10+male,
                 random = list(stay_id = pdDiag(form = ~ bs(log(hr_day), knots = c(log(4)), degree = 1))), 
                 data = mdat,control=list(opt = "optim"))

###########full model with slope control for icuday ==================
fitLME_72h_elix <- lme(sofa_recent  ~  hr_day*ethnicity+U_sofa_adm+age10+male+elix_vws,na.action = na.omit,
                  random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))
###keep 2-3 degrees of poly to use splines.poly(hr_day,3)
fitCox_72h_elix <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+U_sofa_adm+elix_vws,na.action = na.omit,
                    data=mdat.id,x=TRUE)
dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 10,11,12), random = ~ 1, 
              indRandom = 2)
fitJM_pwph_72h_elix <- jointModel(fitLME_72h_elix,fitCox_72h_elix,timeVar="hr_day",
                             method="piecewise-PH-GH",
                             parameterization = "both", derivForm = dForm,
                             control=list(verbose=TRUE,iter.EM=100))
saveRDS(fitJM_pwph_72h_elix, "./fitJM_pwph_72h_elix.rds")


fitLME_72h <- lme(sofa_recent  ~  hr_day*ethnicity+U_sofa_adm+age10+male,na.action = na.omit,
                  random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))
###keep 2-3 degrees of poly to use splines.poly(hr_day,3)
fitCox_72h <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+U_sofa_adm,na.action = na.omit,
                    data=mdat.id,x=TRUE)
dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 9,10,11), random = ~ 1, 
              indRandom = 2)
fitJM_pwph_72h <- jointModel(fitLME_72h,fitCox_72h,timeVar="hr_day",
                             method="piecewise-PH-GH",
                             parameterization = "both", derivForm = dForm,
                             control=list(verbose=TRUE,iter.EM=100))

saveRDS(fitJM_pwph_72h, "./fitJM_pwph_72h.rds")

mdat = dat  %>%
  #filter(Time2LastICU<=48)%>% # from admission to icu in less than 24 hours
  #filter(ICU2Death<=72)%>% # from icu2death in less than 24 hours
  filter(ICU2Death==-1)%>% # from icu2death in less than 24 hours
  filter(ifelse(hr<24*10, sofa_diff_recent!=0, TRUE))%>% 
  filter(ifelse(hr>24*10,hr %% 24 == 0,TRUE))%>%
  arrange(stay_id,hr)

mdat.id = unique(mdat[,c("stay_id","ethnicity","anchor_age","U_sofa_adm","admsofagrp","elix_vws",
                         "U_ICUHR",'male',"Oasis","U_Sepsis","Vent","Lauguage","agegrp","icu_day",
                         "died","age10")])
table1(~anchor_age + U_ICUHR + U_sofa_adm+factor(male)+Oasis+factor(U_Sepsis)+factor(Vent)+elix_vws| ethnicity, 
       data=mdat.id,overall="Total")

rawfit <- survfit(Surv(icu_day,died)~  ethnicity, data=mdat.id)
plot(rawfit, mark.time=FALSE, lty=1:4,col=c(3,4,6,9),
     xlab="Days After ICU admission", ylab="Survival")
legend(1, 1, EthinictyList,lty=1:4,col=c(3,4,6,9), bty='n')

###########full model with slope control for icuday ==================
fitLME_elix <- lme(sofa_recent  ~  hr_day*ethnicity+U_sofa_adm+age10+male+elix_vws,na.action = na.omit,
                   random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))
###keep 2-3 degrees of poly to use splines.poly(hr_day,3)
fitCox_elix <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+U_sofa_adm+elix_vws,na.action = na.omit,
                     data=mdat.id,x=TRUE)
dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 10,11,12), random = ~ 1, 
              indRandom = 2)
fitJM_pwph_elix <- jointModel(fitLME_elix,fitCox_elix,timeVar="hr_day",
                              method="piecewise-PH-GH",
                              parameterization = "both", derivForm = dForm,
                              control=list(verbose=TRUE,iter.EM=100))
saveRDS(fitJM_pwph_elix, "./fitJM_pwph_elix.rds")

fitLME <- lme(sofa_recent  ~  hr_day*ethnicity+U_sofa_adm+age10+male,na.action = na.omit,
              random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))
###keep 2-3 degrees of poly to use splines.poly(hr_day,3)
fitCox <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+U_sofa_adm,na.action = na.omit,
                data=mdat.id,x=TRUE)
dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 9,10,11), random = ~ 1, 
              indRandom = 2)
fitJM_pwph <- jointModel(fitLME,fitCox,timeVar="hr_day",
                         method="piecewise-PH-GH",
                         parameterization = "both", derivForm = dForm,
                         control=list(verbose=TRUE,iter.EM=100))
saveRDS(fitJM_pwph, "./fitJM_pwph.rds")

summary(fitJM_pwph)

mdat$y = mdat$sofa_recent-mdat$U_sofa_adm
knots = c(1,5,133)
fitLME_ns <- lme(sofa_recent ~ ns(hr_day,3)*ethnicity+U_sofa_adm+age10+male,
              random = list(stay_id = pdDiag(form = ~ ns(hr_day,3))), data = mdat,control=list(opt = "optim"))
###keep 2-3 degrees of poly to use splines.poly(hr_day,3)
fitCox_ns <- coxph(Surv(icu_day,died) ~  ethnicity+age10+male+U_sofa_adm,data=mdat.id,x=TRUE)
fitJM_pwph_ns <- jointModel(fitLME_ns,fitCox_ns,timeVar="hr_day",
                             method="piecewise-PH-aGH",
                             control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph_ns)
saveRDS(fitJM_pwph_ns, "./fitJM_pwph_ns3.rds")

ss <- survfitJM(fitJM_pwph, newdata = mdat[mdat$stay_id=='33516489',], idVar = "stay_id", M = 50)
plot(ss,xlim=c(0,60))
plot(ss, include.y = TRUE, add.last.time.axis.tick = TRUE, legend = TRUE,conf.int=TRUE,fill.area=TRUE)

par(mfrow = c(1, 2))
plot(fitJM_pwph_ns, 3, add.KM = TRUE, col = 2, lwd = 2, main = "Spline-approximated Baseline Risk Function")
res.surv <- residuals(fitJM_pwph_ns, process = "Event", type = "Cox")
sfit <- survfit(Surv(res.surv, fitJM_pwph_ns$y$d) ~ 1)
plot(sfit, mark.time = FALSE, conf.int = TRUE, lty = 1:2, 
     xlab = "Cox-Snell Residuals", ylab = "Survival Probability", 
     main = "Piecewise-constant Baseline Risk Function")
curve(pexp(x, lower.tail = FALSE), from = min(res.surv), to = max(res.surv), 
      add = TRUE, col = "red", lwd = 2)
legend("topright", c("Survival function of unit Exponential", 
                     "Survival function of Cox-Snell Residuals"), lty = 1, lwd = 2, col = 2:1, bty = "n")
res.Ypw <- residuals(fitJM_pwph_ns, type = "stand-Subject")
qqnorm(res.Ypw, main = "Piecewise-constant Baseline Risk Function"); qqline(res.Ypw)


ggplot(Ps, aes(hr_day,pred)) + geom_line(aes(colour = ethnicity))
p1 <- ggplot(Ps,aes(x=hr_day,y=pred)) + 
  geom_line(aes(colour = ethnicity)) +
  geom_ribbon(aes(ymin=low,ymax=upp),alpha=0.2,fill="red") +
  scale_y_continuous("y")
p1



dForm <- list(fixed = ~  dns(hr_day, 3)*ethnicity, indFixed = c(2,3,4,11,12,13,14,15,16,17,18,19), 
              random = ~ dns(hr_day, 3) , indRandom = c(2,3,4))
mfX.deriv <- model.frame(terms(dForm$fixed), data = mdat)
Xderiv <- model.matrix(dForm$fixed, mfX.deriv)

bb = fixef(fitLME)[dForm$indFixed])

dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 9,10,11), random = ~ 1, 
              indRandom = 2)

fitJM_pwph_poly <- jointModel(fitLME_poly,fitCox_poly,timeVar="hr_day",
                              method="piecewise-PH-aGH",
                              parameterization = "both", derivForm = dForm,
                              control=list(verbose=TRUE,iter.EM=100))

fitJM_pwph_slp_int <- jointModel(fitLME_slp,fitCox_slp,timeVar="hr_day",
                                 method="piecewise-PH-GH",
                                 parameterization = "both", derivForm = dForm,
                                 interFact = list(value = ~ ethnicity, data = mdat.id),
                                 control=list(verbose=TRUE,iter.EM=200))
summary(fitJM_pwph_slp_int)
anova(fitJM_pwph_slp,fitJM_pwph_slp_int)

require(ggplot2)
library(reshape2)
p_day = seq(from = 1, to = 100, by = 1)
betas = fitJM_pwph_slp$coefficients$betas
sofa_w = betas[1] + p_day*betas[2]
sofa_a = betas[1] + p_day*betas[2]+betas[3]+betas[9]*p_day
sofa_b = betas[1] + p_day*betas[2]+betas[4]+betas[10]*p_day
sofa_h = betas[1] + p_day*betas[2]+betas[5]+betas[11]*p_day
sofa_pred = data.frame(p_day = p_day, white=sofa_w,asian=sofa_a,black=sofa_b,hispanic=sofa_h)
df <- melt(sofa_pred ,  id.vars = 'p_day', variable.name = 'Ethnicity')
ggplot(df, aes(p_day,value)) + geom_line(aes(colour = Ethnicity))





ni <- with(mdat, tapply(sofa_recent, stay_id, length))
summary(ni)

xyplot(stay_id ~ icu_day | ethnicity, data = mdat, groups = stay_id, subset = died == 1,
       type = "b", xlab = "Visit Time (Days)", ylab = "Patients", 
       main = "Visiting Process")

xyplot(sofa_recent ~ icu_day, group = stay_id, data = mdat[mdat$stay_id %in% sample(mdat$stay_id, 10), ], 
       xlab = "days", type = "l", ylim = range(mdat$sofa_recent))

xyplot(sofa_recent ~ icu_day | ethnicity, data = mdat, xlab = "Months", 
       type = c("p", "smooth"), lwd = 2)








###########full model with slope control for total day==================
fitLME_slp <- lme(sofa_recent  ~  hr_day*ethnicity+U_sofa_adm, 
              random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))
fitCox_slp <- coxph(Surv(icu_day,died) ~  ethnicity+age10,data=mdat.id,x=TRUE)
dForm <- list(fixed = ~ 1 + ethnicity, indFixed = c(2, 7,8,9), random = ~ 1, 
              indRandom = 2)
fitJM_pwph_slp <- jointModel(fitLME_slp,fitCox_slp,timeVar="hr_day",method="piecewise-PH-GH",
                         parameterization = "both", derivForm = dForm,
                         control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph_slp)
fitJM_pwph_slp_int <- jointModel(fitLME_slp,fitCox_slp,timeVar="hr_day",method="piecewise-PH-GH",
                          parameterization = "both", derivForm = dForm,
                          interFact = list(value = ~ ethnicity, data = mdat.id),
                          control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph_slp_int)
anova(fitJM_pwph_slp,fitJM_pwph_slp_int)
################model without slope control================
fitLME <- lme(sofa_recent  ~  hr_day+ethnicity+U_sofa_adm, 
              random = ~ 1 | stay_id, data = mdat,control=list(opt = "optim"))
fitCox <- coxph(Surv(icu_day,died) ~  ethnicity+age10,data=mdat.id,x=TRUE)
fitJM_pwph <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",
                         control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph)
fitJM_pwph_int <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",
                             interFact = list(value = ~ ethnicity, data = mdat.id),
                             control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph_int)
anova(fitJM_pwph,fitJM_pwph_int)










fitJM_Con <- jointModel(fitLME,fitCox,timeVar="hr_day",method="Cox-PH-GH",control=list(verbose=TRUE))
summary(fitJM_Con)





table1(~anchor_age + U_ICUHR + U_sofa_adm+male+Oasis+U_Sepsis+Vent| Lauguage, 
       data=dat.id,overall="Total")

table1(~anchor_age + U_ICUHR + U_sofa_adm+male+Oasis+U_Sepsis+Vent| Lauguage, 
       data=mdat.id,overall="Total")

rawfit <- survfit(Surv(icu_day,died)~  Lauguage, data=mdat)
plot(rawfit, mark.time=FALSE, lty=1:4,col=c(3,4,6,9),
     xlab="Days After ICU admission", ylab="Survival")
legend(50, 1, EthinictyList,
       lty=1:4,col=c(3,4,6,9), bty='n')


fitLME <- lme(sofa_recent  ~  hr_day*Lauguage+U_sofa_adm, 
              random = ~ hr_day | stay_id, data = mdat,control=list(opt = "optim"))

fitCox <- coxph(Surv(icu_day,died) ~  Lauguage+age10,
                data=mdat.id,x=TRUE)


fitJM_pwph <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",
                         control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph)

fitJM_pwph_int <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",
                             parameterization = "both", derivForm = dForm,
                             interFact = list(value = ~ ethnicity, data = mdat.id),
                             control=list(verbose=TRUE,iter.EM=100))
summary(fitJM_pwph_int)
anova(fitJM_pwph,fitJM_pwph_int)




fitJM_Con <- jointModel(fitLME,fitCox,timeVar="hr_day",method="Cox-PH-GH",control=list(verbose=TRUE))
summary(fitJM_Con)


basedata = unique(mdat[,c("stay_id","icu_day","died","ethnicity","age10","agegrp","admsofagrp",'U_sofa_adm')])
ldat = mdat[,c("stay_id","hr_day","sofarecgrp","sofa_recent","sofa24grp","sofa_24hours")]
newdata = tmerge(basedata, basedata, id=stay_id, death = event(icu_day, died))
tdata <- tmerge(newdata, ldat, id=stay_id, sofa_tdc = tdc(hr_day, sofarecgrp))
tdfit <- coxph(Surv(tstart, tstop, death) ~ ethnicity+sofa_tdc+agegrp,data=tdata)
summary(tdfit)
print(cox.zph(tdfit))


library(coxme)
tdfitme <- coxme(Surv(tstart, tstop, death) ~ ethnicity+sofa_tdc+agegrp +(1|stay_id),data=tdata)
summary(tdfitme)


library(joineR)  
ldat = mdat[,c("stay_id","hr_day","sofa_recent")]
sdata = UniqueVariables(mdat,var.col = c("icu_day", "died"),id.col = "stay_id")
cdata = UniqueVariables(mdat,var.col = c("ethnicity", "U_sofa_adm","age10"),id.col = "stay_id")
jdata = jointdata(
  longitudinal = ldat,
  survival = sdata,
  baseline = cdata,
  id.col = "stay_id",
  time.col = "hr_day"
)
fit <- joint(data = jdata,
             long.formula = sofa_recent ~ 1 + hr_day*ethnicity+U_sofa_adm+age10,
             surv.formula = Surv(U_LOS, died) ~ ethnicity,
             model = "intslope")
summary(fit)


heart.surv <- UniqueVariables(heart.valve,
                              var.col = c("fuyrs", "status"),
                              id.col = "num")
heart.long <- heart.valve[, c("num", "time", "log.lvmi")]
heart.cov <- UniqueVariables(heart.valve,
                             c("age", "hs", "sex"),
                             id.col = "num")
heart.valve.jd <- jointdata(longitudinal = heart.long,
                            baseline = heart.cov,
                            survival = heart.surv,
                            id.col = "num",
                            time.col = "time")
fit <- joint(data = heart.valve.jd,
             long.formula = log.lvmi ~ 1 + time + hs,
             surv.formula = Surv(fuyrs, status) ~ hs,
             model = "intslope")
summary(fit)

km= survfit(Surv(U_LOS, died) ~ ethnicity, data = mdat)
ggsurvplot(km, data = mdat, risk.table = FALSE)


fit2 <- mjoint(
  formLongFixed = list("sofa" = sofa_recent ~ hr_day+ethnicity+U_sofa_adm+male+age10+Vent),
  formLongRandom = list("sofa" = ~ hr_day | stay_id),
  formSurv = Surv(U_LOS,died)~  ethnicity + age10 + U_sofa_adm,
  data = mdat,
  timeVar = "hr_day",
  verbose = TRUE)


fitJMBayes <- jointModelBayes(fitLME_recent,fitCox,timeVar="hr_day",method="piecewise-PH-GH",control=list(verbose=TRUE))


fitLME_24h <- lme(sofa_24hours  ~ hr_day+ethnicity+U_sofa_adm+male+age10+Vent, 
                  random = ~ hr_day | stay_id, data = mdat, control = lmeControl(opt = "optim"))

fitJM_24h <- jointModel(fitLME_24h,fitCox,timeVar="hr_day",method="piecewise-PH-GH",control=list(verbose=TRUE))
summary(fitJM_24h)


fitCox <- coxph(Surv(U_LOS,died)~  ethnicity + anchor_age,
                data=arrange(unique(mdat[,c('stay_id','ethnicity','anchor_age','icu_day','died','U_LOS')]),stay_id),x=TRUE)

fitJM <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",control=list(verbose=TRUE))


fitJM <- jointModel(fitLME,fitCox,timeVar="hr_day",method="piecewise-PH-GH",
                    control=list(tol1=1e-7,tol2=1e-7,iter.EM=1,only.EM=FALSE,optimizer="optim",verbose=TRUE,GKk=15))
summary(fitJM)
plot(fitJM)

dat_person = dat
dat_person$hr = NULL
dat_person$sofa_24hours = NULL
dat_unique_person = unique(dat_person)

fitLME.null <- lme(sofa_24hours ~ hr+gender+ethnicity+Oasis+sepsis+Vent+male,
                   random = ~ 1 | stay_id, data = dat)


dat$adm_sofa_grp = cut(dat$adm_SOFA,breaks=c(0,6,9,12,Inf), include.lowest=TRUE,
                 labels=c("0-6","7-9","10-12",">12"))
dat$sofa_24hr_grp = cut(dat$sofa_hr24,breaks=c(0,6,9,12,Inf), include.lowest=TRUE,
                       labels=c("0-6","7-9","10-12",">12"))
dat$sofa_48hr_grp = cut(dat$sofa_hr48,breaks=c(0,6,9,12,Inf), include.lowest=TRUE,
                        labels=c("0-6","7-9","10-12",">12"))

WeibullReg(Surv(Los, death) ~ factor(ethnicity) + anchor_age+adm_SOFA+Oasis+Vent_Day1, data=filter(dat,Los>1))
WeibullDiag(Surv(Los, death) ~ factor(ethnicity), data = filter(dat,Los>1), labels=c("WHITE", "AFRICAN AMERICAN", "ASIAN", "HISPANIC/LATINO"))


ttest_adm_sofa = c((t.test(dat$adm_SOFA[dat$ethnicity=='WHITE'],dat$adm_SOFA[dat$ethnicity=='ASIAN']))$p.value,
                  (t.test(dat$adm_SOFA[dat$ethnicity=='WHITE'],dat$adm_SOFA[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                  (t.test(dat$adm_SOFA[dat$ethnicity=='WHITE'],dat$adm_SOFA[dat$ethnicity=='HISPANIC/LATINO']))$p.value)
ttest_age =  c((t.test(dat$anchor_age[dat$ethnicity=='WHITE'],dat$anchor_age[dat$ethnicity=='ASIAN']))$p.value,
               (t.test(dat$anchor_age[dat$ethnicity=='WHITE'],dat$anchor_age[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
               (t.test(dat$anchor_age[dat$ethnicity=='WHITE'],dat$anchor_age[dat$ethnicity=='HISPANIC/LATINO']))$p.value)

ttest_oasis =  c((t.test(dat$Oasis[dat$ethnicity=='WHITE'],dat$Oasis[dat$ethnicity=='ASIAN']))$p.value,
                   (t.test(dat$Oasis[dat$ethnicity=='WHITE'],dat$Oasis[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                   (t.test(dat$Oasis[dat$ethnicity=='WHITE'],dat$Oasis[dat$ethnicity=='HISPANIC/LATINO']))$p.value)
wilcox_los =  c((wilcox.test(dat$Los[dat$ethnicity=='WHITE'],dat$Los[dat$ethnicity=='ASIAN']))$p.value,
                 (wilcox.test(dat$Los[dat$ethnicity=='WHITE'],dat$Los[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                 (wilcox.test(dat$Los[dat$ethnicity=='WHITE'],dat$Los[dat$ethnicity=='HISPANIC/LATINO']))$p.value)


                 
ztest_died_in_icu = c(prop.test(c(sum(dat$Died_in_ICU[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Died_in_ICU[dat$ethnicity=='ASIAN'],na.rm = TRUE)),
                           c(length(na.omit(dat$Died_in_ICU[dat$ethnicity=='WHITE'])),length(na.omit(dat$Died_in_ICU[dat$ethnicity=='ASIAN']))))$p.value,
                 prop.test(c(sum(dat$Died_in_ICU[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Died_in_ICU[dat$ethnicity=='AFRICAN AMERICAN'],na.rm = TRUE)),
                           c(length(na.omit(dat$Died_in_ICU[dat$ethnicity=='WHITE'])),length(na.omit(dat$Died_in_ICU[dat$ethnicity=='AFRICAN AMERICAN']))))$p.value,
                 prop.test(c(sum(dat$Died_in_ICU[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Died_in_ICU[dat$ethnicity=='HISPANIC/LATINO'],na.rm = TRUE)),
                           c(length(na.omit(dat$Died_in_ICU[dat$ethnicity=='WHITE'])),length(na.omit(dat$Died_in_ICU[dat$ethnicity=='HISPANIC/LATINO']))))$p.value)
ztest_vent_firstday = c(prop.test(c(sum(dat$Vent_Day1[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Vent_Day1[dat$ethnicity=='ASIAN'],na.rm = TRUE)),
                                c(length(na.omit(dat$Vent_Day1[dat$ethnicity=='WHITE'])),length(na.omit(dat$Vent_Day1[dat$ethnicity=='ASIAN']))))$p.value,
                      prop.test(c(sum(dat$Vent_Day1[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Vent_Day1[dat$ethnicity=='AFRICAN AMERICAN'],na.rm = TRUE)),
                                c(length(na.omit(dat$Vent_Day1[dat$ethnicity=='WHITE'])),length(na.omit(dat$Vent_Day1[dat$ethnicity=='AFRICAN AMERICAN']))))$p.value,
                      prop.test(c(sum(dat$Vent_Day1[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Vent_Day1[dat$ethnicity=='HISPANIC/LATINO'],na.rm = TRUE)),
                                c(length(na.omit(dat$Vent_Day1[dat$ethnicity=='WHITE'])),length(na.omit(dat$Vent_Day1[dat$ethnicity=='HISPANIC/LATINO']))))$p.value)


table(dat$language)
sum(is.na(dat$Sepsis[dat$language=='ENGLISH']))
sum(is.na(dat$Sepsis[dat$language=='Non-Engl']))


tot_ethnicity = table(dat$ethnicity)
sepsis_missing_ethnicity = c(sum(is.na(dat$Sepsis[dat$ethnicity=='WHITE'])),
                             sum(is.na(dat$Sepsis[dat$ethnicity=='ASIAN'])),
                             sum(is.na(dat$Sepsis[dat$ethnicity=='AFRICAN AMERICAN'])),
                             sum(is.na(dat$Sepsis[dat$ethnicity=='HISPANIC/LATINO'])))
for (eth in EthinictyList) 
{
  print(eth)
  print('los--')
  print(quantile(dat$Los[dat$ethnicity==eth]))
  print('na in sepsis')
  print(sum(is.na(dat$Sepsis[dat$ethnicity==eth])))
  
}
ztest_sepsis = c(prop.test(c(sum(dat$Sepsis[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Sepsis[dat$ethnicity=='ASIAN'],na.rm = TRUE)),
                           c(length(na.omit(dat$Sepsis[dat$ethnicity=='WHITE'])),length(na.omit(dat$Sepsis[dat$ethnicity=='ASIAN']))))$p.value,
                 prop.test(c(sum(dat$Sepsis[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Sepsis[dat$ethnicity=='AFRICAN AMERICAN'],na.rm = TRUE)),
                           c(length(na.omit(dat$Sepsis[dat$ethnicity=='WHITE'])),length(na.omit(dat$Sepsis[dat$ethnicity=='AFRICAN AMERICAN']))))$p.value,
                 prop.test(c(sum(dat$Sepsis[dat$ethnicity=='WHITE'],na.rm = TRUE),sum(dat$Sepsis[dat$ethnicity=='HISPANIC/LATINO'],na.rm = TRUE)),
                           c(length(na.omit(dat$Sepsis[dat$ethnicity=='WHITE'])),length(na.omit(dat$Sepsis[dat$ethnicity=='HISPANIC/LATINO']))))$p.value)

wilcox_pctavgsofa = c((wilcox.test(dat$pct_beyond_avg_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_avg_sofa[dat$ethnicity=='ASIAN']))$p.value,
                      (wilcox.test(dat$pct_beyond_avg_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_avg_sofa[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                      (wilcox.test(dat$pct_beyond_avg_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_avg_sofa[dat$ethnicity=='HISPANIC/LATINO']))$p.value)
wilcox_pct75sofa = c((wilcox.test(dat$pct_beyond_p75_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_p75_sofa[dat$ethnicity=='ASIAN']))$p.value,
                     (wilcox.test(dat$pct_beyond_p75_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_p75_sofa[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                     (wilcox.test(dat$pct_beyond_p75_sofa[dat$ethnicity=='WHITE'],dat$pct_beyond_p75_sofa[dat$ethnicity=='HISPANIC/LATINO']))$p.value)
wilcox_avgsofa =  c((wilcox.test(dat$avg_sofa[dat$ethnicity=='WHITE'],dat$avg_sofa[dat$ethnicity=='ASIAN']))$p.value,
                    (wilcox.test(dat$avg_sofa[dat$ethnicity=='WHITE'],dat$avg_sofa[dat$ethnicity=='AFRICAN AMERICAN']))$p.value,
                    (wilcox.test(dat$avg_sofa[dat$ethnicity=='WHITE'],dat$avg_sofa[dat$ethnicity=='HISPANIC/LATINO']))$p.value)


