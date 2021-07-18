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
eicu <- read.csv("eicu_died.csv")%>%
  mutate(age_int = ifelse(age=="> 89",89,as.numeric(age)))%>%
  mutate(bmi = ifelse(admissionheight>10,10000*admissionweight/admissionheight^2,NA))%>%
  mutate(age_cat = ifelse(age=="> 89",">89","<=89"))


table1(~age_int  + sofa_24hour+ factor(gender)+apache_iv+icu_los_hours| ethnicity, 
       data=eicu,overall="Total")

table1(~age_int +age_cat + sofa_24hour+ factor(gender)+apache_iv+icu_los_hours
       +apachedxgroup+bmi
       | ethnicity, 
       data=filter(eicu,(unitdischargelocation=="Death")&(icu_los_hours>24)&((ethnicity=="African American")
                                             |(ethnicity=="Hispanic")
                                             |(ethnicity=="Caucasian")
                                             |(ethnicity=="Asian"))),overall="Total")

eicu_sofa_24hr <- read_csv("eicu_sofa_24hr.csv")%>%
  left_join(eicu,by = c('patientunitstayid'))%>%
  mutate(sofa_diff = sofa_24hours-sofa_24hour)%>%
  filter((hr<=120)&((ethnicity=="African American")
                    |(ethnicity=="Hispanic")
                    |(ethnicity=="Caucasian")
                    |(ethnicity=="Asian")))


Sum = groupwiseMean(sofa_24hours ~ hr + ethnicity,
                    data   = eicu_sofa_24hr,
                    conf   = 0.95,
                    digits = 3)
Sum$hr = as.factor(Sum$hr)
p1 = ggplot(Sum, aes(x=hr, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average SOFA score in 120 hrs after admission for people who survived more than 120 hrs")
Sum = groupwiseMean(sofa_diff ~ hr + ethnicity,
                    data   = eicu_sofa_24hr,
                    conf   = 0.95,
                    digits = 3)
Sum$hr = as.factor(Sum$hr)
p2 = ggplot(Sum, aes(x=hr, y=Mean,color = ethnicity))+
  geom_point(position = position_dodge(width = 0.6))+
  geom_errorbar(aes(hr, ymin = Trad.lower, ymax = Trad.upper), 
                position = position_dodge(width = 0.6))+
  theme_bw()+ggtitle("Average relative change of SOFA score in 120 hrs after admission for people who survived more than 120 hrs")

grid.arrange(p1, p2, nrow = 2)


