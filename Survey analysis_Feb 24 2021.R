library(ggplot2)
library(tidyr)
library(tidyverse)
library(psych) # psychometrics
library(caret) # ML workflow
library(MVN) # multivariate normality
library(sjPlot) #visualization
library(fastDummies) # do i use this?***
library(scales)  # for percentage scales
library(glmnet)
library(nnet)
library(leaps) # stepwise regression
library(car) # for VIF function, also recoding
library(RColorBrewer)
library(ggpubr)
library(questionr)
library(broom)


setwd("~/Desktop/Stuff/PhD/Research projects/Genome BC Perceptions of Gene Editing in Agriculture/SURVEY ANALYSIS/Survey analysis")
survey <- read.csv('Perceptions of gene editing in agriculture_August 3 2020 EXTRA CLEAN numeric.csv', sep= ',', stringsAsFactors=FALSE, header=TRUE)


##Demographics ETC####
#cleaningfor regression analysis
gender_full <-as.factor(survey$Q26_clean)
political_full <- survey$Q29
religious_full <-survey$Q32
uni_full <-survey$Q28
race_full <-survey$Q30_clean
age_full <-survey$Q27
familiar_gm <-dplyr::na_if(survey$Q1_1, 6)
familiar_ge <-dplyr::na_if(survey$Q1_2, 6)
familiar_gd <-dplyr::na_if(survey$Q1_3, 6)
similar <-survey$Q2
income <-survey$Q33
natural <- dplyr::na_if(survey$Q23_1, 6)


#######cleaning for graphs####
#gender (simplified in XLS)
survey$gender <- as.factor(survey$Q26_clean)
levels(survey$gender) <- c("1" = "Female",
                    "2" = "Male")

###run these for graphs
#political, simplified
survey <- survey %>%
  mutate(political = case_when(
    Q29 <=2 ~ "Liberal",
    Q29 ==3 ~ "Moderate",
    Q29 >=4 ~ "Conservative"
  ))
survey$political <- as.factor(survey$political)
political <- survey$political

#religious, simplified
survey <- survey %>% 
  mutate(religious = ifelse(Q32<=3, "religious", "non-religious"))
survey$religious <- as.factor(survey$religious)
religious <- survey$religious

#education, simplified (university degree) 
survey <- survey %>% 
  mutate(uni = ifelse(Q28<=3, "uni", "no uni"))
survey$uni <- as.factor(survey$uni)
uni <- survey$uni

#race, simplified
survey <- survey %>% 
  mutate(race = ifelse(Q30==6, "white", "non-white"))
survey$race <- as.factor(survey$race)
race <- survey$race

#age, simplified
survey <- survey %>%
  mutate(age = case_when(
    Q27_clean >=65 ~ "Older",
    Q27_clean <65 & Q27_clean>35 ~ "Middle",
    Q27_clean <=35 ~ "Young"
  ))
survey$age <- as.factor(survey$age)
age <-survey$age


#FOR GRAPHS


race_count <- as.factor(survey$Q30_clean)
levels(race_count) <- c("1" = "East Asian",
                        "2" = "Black",
                        "3" = "American Indian, First Nation, or Pacific Islander",
                        "4" = "Hispanic",
                        "5" = "Mixed race",
                        "6" = "White",
                        "7" = "Didn't say",
                        "8" = "South Asian",
                        "9" = "Middle Eastern")
summary(race_count)


#### Factor analyses for regressions *NOT full exploratory analyses, just  results for regressions* ##########

# CLIMATE SCALE
#remove NA's
survey$Q21_1.na <-na_if(survey$Q21_1, 6)
survey$Q21_2.na <-na_if(survey$Q21_2, 6)
survey$Q21_3.na <-na_if(survey$Q21_3, 6)
survey$Q21_4.na <-na_if(survey$Q21_4, 6)
survey$Q21_5.na <-na_if(survey$Q21_5, 6)
survey$Q21_6.na <-na_if(survey$Q21_6, 6)
survey$Q21_7.na <-na_if(survey$Q21_7, 6)

# Load data
climate <- data.frame(Q21_1.na, Q21_2.na, Q21_3.na, Q21_4.na, Q21_5.na, Q21_6.na, Q21_7.na)
keys.climate <- c(1,-1,-1, 1, 1, 1, -1)
reverse.code(keys.climate, climate)

# Factor analysis
climate.fa <- fa(climate, nfactors=2, scores=TRUE, rotate="varimax", fm="pa")
print.psych(climate.fa, cut = 0.3, sort = TRUE)
fa.diagram(climate.fa)


# Cronbach alpha scores
climate.a <- data.frame(Q21_6.na, Q21_4.na, Q21_1.na, Q21_5.na)
climate.b <- data.frame(Q21_2.na, Q21_3.na, Q21_7.na)
psych::alpha(climate.a, check.keys=TRUE)
psych::alpha(climate.b, check.keys=TRUE)


# split climate factors and plot responses to DVs

#add scores to original dataset
survey$ccscores1 <- data.frame(climate.fa$scores[,1])
survey$ccscores2 <- data.frame(climate.fa$scores[,2])

survey <- survey %>%
  mutate(cc1 = ifelse(ccscores1 >0, 1, -1))
survey$cc1 <- as.factor(survey$cc1)

survey <- survey %>%
  mutate(cc2 = ifelse(ccscores2 >0, 1, -1))
survey$cc2 <- as.factor(survey$cc2)

##ALTERNATIVE APPROACH--mean scores.

climate1mean=rowMeans(cbind(survey$Q21_6.na, survey$Q21_4.na, survey$Q21_1.na, survey$Q21_5.na),na.rm=TRUE)
climate2mean=rowMeans(cbind(survey$Q21_2.na, survey$Q21_3.na, survey$Q21_7.na),na.rm=TRUE)

climatebothmean = rowMeans(cbind(survey$Q21_6.na, survey$Q21_4.na, survey$Q21_1.na, survey$Q21_5.na, survey$Q21_2.na, survey$Q21_3.na, survey$Q21_7.na),na.rm=TRUE)


# GLOBALIZATION SCALE

#remove NA's
survey$Q24_1.na <-na_if(survey$Q24_1, 6)
survey$Q24_3.na <-na_if(survey$Q24_3, 6)
survey$Q24_4.na <-na_if(survey$Q24_4, 6)

# Load data
global <- data.frame(Q24_1.na, Q24_3.na, Q24_4.na)

# Factor analysis using Psych package
global.fa <- fa(global, nfactors=1, scores=TRUE,rotate="varimax", fm="pa")
print.psych(global.fa, cut = 0.3, sort = TRUE)
fa.diagram(global.fa)

global.a <- data.frame(Q24_1.na, Q24_4.na, Q24_3.na)
psych::alpha(global.a)

#add scores to original dataset
survey$globalscores1 <- data.frame(global.fa$scores[,1])

survey <- survey %>%
  mutate(global1 = ifelse(globalscores1 >0, 1, -1))
survey$global1 <- as.factor(survey$global1)

#alternative means approach--
globalmean=rowMeans(cbind(survey$Q24_1.na, survey$Q24_4.na, survey$Q24_3.na),na.rm=TRUE)


# TRUST SCALE

#remove NA's
survey$Q20_1.na <-na_if(survey$Q20_1, 6)
survey$Q20_2.na <-na_if(survey$Q20_2, 6)
survey$Q20_3.na <-na_if(survey$Q20_3, 6)

trust <- data.frame(Q20_1.na, Q20_2.na, Q20_3.na)


trust.fa <- fa(trust, nfactors=1, scores=TRUE,rotate="varimax", fm="pa")
print.psych(trust.fa, cut = 0.3, sort = TRUE)


trust.a <- data.frame(Q20_1.na, Q20_2.na, Q20_3.na)
psych::alpha(trust.a)


#add scores to original dataset
survey$trustscores1 <- data.frame(trust.fa$scores[,1])

survey <- survey %>%
  mutate(trust1 = ifelse(trustscores1 >0, 1, -1))
survey$trust1 <- as.factor(survey$trust1)

#alternative approach
trustmean=rowMeans(cbind(survey$Q20_1.na, survey$Q20_2.na, survey$Q20_3.na),na.rm=TRUE)


# TECH SCALE

#remove NA's
survey$Q22_2.na <-na_if(survey$Q22_2, 6)
survey$Q22_5.na <-na_if(survey$Q22_5, 6)

tech <- data.frame(Q22_2.na, Q22_5.na)

tech.fa <- fa(tech, nfactors=1, scores=TRUE,rotate="varimax", fm="pa")
print.psych(tech.fa, cut = 0.3, sort = TRUE)

tech.a <- data.frame(Q22_5.na, Q22_2.na)
psych::alpha(tech.a)

#add scores to original dataset
survey$techscores1 <- data.frame(tech.fa$scores[,1])

survey <- survey %>%
  mutate(tech1 = ifelse(techscores1 >0, 1, -1))
survey$tech1 <- as.factor(survey$tech1)

#tech mean
techmean=rowMeans(cbind(survey$Q22_5.na, survey$Q22_2.na),na.rm=TRUE)


# TOXICS SCALE

#remove NA's
survey$Q23_2.na <-na_if(survey$Q23_2, 6)
survey$Q23_3.na <-na_if(survey$Q23_3, 6)
survey$Q23_5.na <-na_if(survey$Q23_5, 6)

toxics <- data.frame(Q23_2.na, Q23_3.na, Q23_5.na)

scree(toxics)
toxics.fa <- fa(toxics, nfactors=1, scores=TRUE,rotate="varimax", fm="pa")
print.psych(toxics.fa, cut = 0.3, sort = TRUE)

toxics.a <- data.frame(Q23_3.na, Q23_2.na, Q23_5.na)
psych::alpha(toxics.a)

survey$toxicsscores1 <- data.frame(toxics.fa$scores[,1])

survey <- survey %>%
  mutate(toxics1 = ifelse(toxicsscores1 >0, 1, -1))
survey$toxics1 <- as.factor(survey$toxics1)

#mean score
toxicmean=rowMeans(cbind(survey$Q23_2.na, survey$Q23_3.na, survey$Q23_5.na),na.rm=TRUE)


# GREEN REVOLUTION SCALE

#remove NA's

#remove observations with unfamiliarity with GR
# survey[!survey$Q25_1.na=='1',]
# survey[!survey$Q25_1.na=='2',]

survey$Q25_2.na <-na_if(survey$Q25_2, 6)
survey$Q25_3.na <-na_if(survey$Q25_3, 6)
survey$Q25_4.na <-na_if(survey$Q25_4, 6)
survey$Q25_5.na <-na_if(survey$Q25_5, 6)
survey$Q25_6.na <-na_if(survey$Q25_6, 6)
survey$Q25_7.na <-na_if(survey$Q25_7, 6)
survey$Q25_8.na <-na_if(survey$Q25_8, 6)

green <- data.frame(survey$Q25_2.na, survey$Q25_3.na, survey$Q25_4.na, survey$Q25_5.na, survey$Q25_6.na, survey$Q25_7.na, survey$Q25_8.na)
keys.green <- c(1, 1, -1, 1, -1, -1, -1)
reverse.code(keys.green, green)

green.fa <- fa(green, nfactors=2, scores=TRUE,rotate="varimax", fm="pa")
print.psych(green.fa, cut = 0.3, sort = TRUE)

green.a <- data.frame(survey$Q25_4.na, survey$Q25_6.na, survey$Q25_8.na, survey$Q25_7.na)
green.b <- data.frame(survey$Q25_5.na, survey$Q25_2.na, survey$Q25_3.na)

psych::alpha(green.a, check.keys=TRUE)
psych::alpha(green.b, check.keys=TRUE)


#add scores to original dataset (ASK GUILLAUME ABOUT THIS?)
survey$greenscores1 <- data.frame(green.fa$scores[,1])
survey$greenscores2 <- data.frame(green.fa$scores[,2])

survey <- survey %>%
  mutate(green1 = ifelse(greenscores1 >0, 1, -1))
survey$green1 <- as.factor(survey$green1)

survey <- survey %>%
  mutate(green2 = ifelse(greenscores2 >0, 1, -1))
survey$green2 <- as.factor(survey$green2)

green1mean=rowMeans(cbind(survey$Q25_4.na, survey$Q25_6.na, survey$Q25_8.na, survey$Q25_7.na),na.rm=TRUE)
green2mean=rowMeans(cbind(survey$Q25_5.na, survey$Q25_2.na, survey$Q25_3.na),na.rm=TRUE)

greenbothmean = rowMeans(cbind(survey$Q25_5.na, survey$Q25_2.na, survey$Q25_3.na, survey$Q25_4.na, survey$Q25_6.na, survey$Q25_8.na, survey$Q25_7.na),na.rm=TRUE)

survey$green1mean <-green1mean
survey$green2mean <-green2mean
survey$greenbothmean <-greenbothmean


###### DV1 data cleaning #####
# Combine and label variables, convert "don't know/not sure" responses to NA's


# Tomatoes
# pre-nudge
survey$Q4.na <- na_if(survey$Q4, 6)
# post-nudge
Q6.na <- na_if(survey$Q6, 6)
# change from pre- to post-nudge
tomato.change <- Q6.na - Q4.na

# Cattle
# pre-nudge
survey$Q8.na <- na_if(survey$Q8, 6)
# post-nudge
Q10.na <- na_if(survey$Q10, 6)
#change from pre- to post-nudge
cattle.change <- Q10.na - Q8.na

# Wheat
# pre-nudge
survey$Q12.na <- na_if(survey$Q12, 6)
# post-nudge
Q14.na <- na_if(survey$Q14, 6)
#change from pre- to post-nudge
wheat.change <- Q14.na - Q12.na

#RECODE???
#library(car)
#Q4.na = recode(Q4.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q6.na = recode(Q6.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q8.na = recode(Q8.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q10.na = recode(Q10.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q12.na = recode(Q12.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q14.na = recode(Q14.na, "1=5; 2=4; 3=3; 4=2; 5=1") 

comfort=rowMeans(cbind(survey$Q4.na, survey$Q8.na, survey$Q12.na),na.rm=TRUE)


#### DV1 bar charts #############

###examine plots of those who are both in agreement with both GR1 and GR2
tomato.greenboth <- survey %>% 
  filter(green1mean>'3' & green2mean>'3') %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) 

tomato.greenneither <- survey %>%
  filter(green1mean<'3' & green2mean<'3') %>%
         count(Q4.na)

tomato.greenskeptic <- survey %>%
  filter(green1mean>'3' & green2mean<'3') %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) 

tomato.greenoptimist <- survey %>%
  filter(green1mean<'3' & green2mean>'3') %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) 


tomatoboth.plot<-ggplot(data=tomato.greenboth, aes(x=Q4.na, y=percent)) +
  geom_bar(stat="identity")
tomatoboth.plot

tomatoneither.plot<-ggplot(data=tomato.greenneither, aes(x=Q4.na, y=n)) +
  geom_bar(stat="identity")

tomatogreenskeptic.plot<-ggplot(data=tomato.greenskeptic, aes(x=Q4.na, y=percent)) +
  geom_bar(stat="identity")
tomatogreenskeptic.plot
  

tomatogreenoptimist.plot<-ggplot(data=tomato.greenoptimist, aes(x=Q4.na, y=percent)) +
  geom_bar(stat="identity")
tomatogreenoptimist.plot




cattle.greenboth <- survey %>% 
  filter(green1mean>'3' & green2mean>'3') %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) 


cattle.greenoptimist <- survey %>%
  filter(green1mean<'3' & green2mean>'3') %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) 

cattleboth.plot<-ggplot(data=cattle.greenboth, aes(x=Q8.na, y=percent)) +
  geom_bar(stat="identity")
cattleboth.plot

cattlegreenoptimist.plot<-ggplot(data=cattle.greenoptimist, aes(x=Q8.na, y=percent)) +
  geom_bar(stat="identity")
cattlegreenoptimist.plot



wheat.greenboth <- survey %>% 
  filter(green1mean>'3' & green2mean>'3') %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) 

wheat.greenoptimist <- survey %>%
  filter(green1mean<'3' & green2mean>'3') %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) 

wheatboth.plot<-ggplot(data=wheat.greenboth, aes(x=Q12.na, y=percent)) +
  geom_bar(stat="identity")
wheatboth.plot

wheatgreenoptimist.plot<-ggplot(data=cattle.greenoptimist, aes(x=Q12.na, y=percent)) +
  geom_bar(stat="identity")
wheatgreenoptimist.plot






cattle.greenboth <- survey %>% 
  filter(green1mean>'3' & green2mean>'3') %>%
  count(Q8.na)

cattle.greenneither <- survey %>%
  filter(green1mean<'3' & green2mean<'3') %>%
  count(Q8.na)

cattle.greenskeptic <- survey %>%
  filter(green1mean>'3' & green2mean<'3') %>%
  count(Q8.na)

cattle.greenoptimist <- survey %>%
  filter(green1mean<'3' & green2mean>'3') %>%
  count(Q8.na)



#of those who are comfortable, people were not in agreement about 
# GR skepticism.
comfort <- survey %>% 
  filter(as.numeric(Q4.na)<'3') %>%
  count(green1mean)
comfort.plot<-ggplot(data=comfort, aes(x=green1mean, y=n)) +
  geom_bar(stat="identity")

# Of those who were comfortable, people were in agreement 
# and largely optimistic about the GR
comfort <- survey %>% 
  filter(as.numeric(Q4.na)<'3') %>%
  count(green2mean)
comfort2.plot<-ggplot(data=comfort, aes(x=green2mean, y=n)) +
  geom_bar(stat="identity")





survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")

tomato_plot<- ggplot(data=tomato.greenboth, aes(Q4.na, fill=as.factor(Q4.na)))+
  geom_histogram()
tomato_plot

########GENDER
# Tomatoes subset by gender, pre-nudge
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                "2" = "Comfortable",
                "3" = "Neither",
                "4" = "Uncomfortable",
                "5" = "Very uncomfortable")
tomato.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.gender, aes(x = Q4.na, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomatoes subset by gender, post-nudge
survey$Q6.na <- as.factor(Q6.na)
levels(Q6.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q6.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.gender, aes(x = Q4.na, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by gender
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.change.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.gender, aes(x = tomato.change, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by gender
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.gender, aes(x = Q8.na, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by gender
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.gender, aes(x = cattle.change, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by gender
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
wheat.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.gender, aes(x = Q12.na, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by gender
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
wheat.change.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.gender, aes(x = wheat.change, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar




########POLITICAL
# Tomatoes subset by political
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.political, aes(x = Q4.na, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by political
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
tomato.change.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.political, aes(x = tomato.change, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by political
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.political, aes(x = Q8.na, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Conservative", "Moderate", "Liberal"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by political
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.political, aes(x = cattle.change, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by political
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                    "2" = "Comfortable",
                    "3" = "Neither",
                    "4" = "Uncomfortable",
                    "5" = "Very uncomfortable")
wheat.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.political, aes(x = Q12.na, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by gender
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                          "2" = "Comfortable",
                          "3" = "Neither",
                          "4" = "Uncomfortable",
                          "5" = "Very uncomfortable")
wheat.change.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.political, aes(x = wheat.change, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar






########RELIGIOUS
# Tomatoes subset by religious
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.religious, aes(x = Q4.na, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by religious
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
tomato.change.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.religious, aes(x = tomato.change, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by religious
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.religious, aes(x = Q8.na, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by religious
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.religious, aes(x = cattle.change, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by religious
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                    "2" = "Comfortable",
                    "3" = "Neither",
                    "4" = "Uncomfortable",
                    "5" = "Very uncomfortable")
wheat.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.religious, aes(x = Q12.na, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by religious
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                          "2" = "Comfortable",
                          "3" = "Neither",
                          "4" = "Uncomfortable",
                          "5" = "Very uncomfortable")
wheat.change.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.religious, aes(x = wheat.change, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar





########UNI
# Tomatoes subset by uni
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.uni, aes(x = Q4.na, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by uni
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
tomato.change.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.uni, aes(x = tomato.change, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by uni
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.uni, aes(x = Q8.na, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by uni
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.uni, aes(x = cattle.change, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by uni
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                    "2" = "Comfortable",
                    "3" = "Neither",
                    "4" = "Uncomfortable",
                    "5" = "Very uncomfortable")
wheat.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.uni, aes(x = Q12.na, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by uni
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                          "2" = "Comfortable",
                          "3" = "Neither",
                          "4" = "Uncomfortable",
                          "5" = "Very uncomfortable")
wheat.change.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.uni, aes(x = wheat.change, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar







########RACE
# Tomatoes subset by race
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.race, aes(x = Q4.na, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by race
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
tomato.change.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.race, aes(x = tomato.change, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by race
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.race, aes(x = Q8.na, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by race
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.race, aes(x = cattle.change, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by race
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                    "2" = "Comfortable",
                    "3" = "Neither",
                    "4" = "Uncomfortable",
                    "5" = "Very uncomfortable")
wheat.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.race, aes(x = Q12.na, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by race
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                          "2" = "Comfortable",
                          "3" = "Neither",
                          "4" = "Uncomfortable",
                          "5" = "Very uncomfortable")
wheat.change.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.race, aes(x = wheat.change, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar







########AGE
# Tomatoes subset by age
survey$Q4.na <- as.factor(Q4.na)
levels(Q4.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
tomato.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.bar <-ggplot(tomato.age, aes(x = Q4.na, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.bar

# Tomato CHANGE (post-nudge) subset by age
survey$tomato.change <- as.factor(tomato.change)
levels(tomato.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
tomato.change.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.bar <-ggplot(tomato.change.age, aes(x = tomato.change, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.bar




# Cattle subset by age
survey$Q8.na <- as.factor(Q8.na)
levels(Q8.na) <- c("1" = "Very comfortable",
                   "2" = "Comfortable",
                   "3" = "Neither",
                   "4" = "Uncomfortable",
                   "5" = "Very uncomfortable")
cattle.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.bar <-ggplot(cattle.age, aes(x = Q8.na, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.bar

# Cattle CHANGE (post-nudge) subset by age
survey$cattle.change <- as.factor(cattle.change)
levels(cattle.change) <- c("1" = "Very comfortable",
                           "2" = "Comfortable",
                           "3" = "Neither",
                           "4" = "Uncomfortable",
                           "5" = "Very uncomfortable")
cattle.change.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.bar <-ggplot(cattle.change.age, aes(x = cattle.change, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.bar



# Wheat subset by age
survey$Q12.na <- as.factor(Q12.na)
levels(Q12.na) <- c("1" = "Very comfortable",
                    "2" = "Comfortable",
                    "3" = "Neither",
                    "4" = "Uncomfortable",
                    "5" = "Very uncomfortable")
wheat.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.bar <-ggplot(wheat.age, aes(x = Q12.na, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.bar

# Wheat CHANGE (post-nudge) subset by age
survey$wheat.change <- as.factor(wheat.change)
levels(wheat.change) <- c("1" = "Very comfortable",
                          "2" = "Comfortable",
                          "3" = "Neither",
                          "4" = "Uncomfortable",
                          "5" = "Very uncomfortable")
wheat.change.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.bar <-ggplot(wheat.change.age, aes(x = wheat.change, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.bar


###SCALES

#CLIMATE SCALE 1: TOMATO
tomato.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.cc1 <- tomato.cc1[1:12,]
tomato.cc1.bar <-ggplot(tomato.cc1, aes(x = Q4.na, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.cc1.bar

#CLIMATE SCALE 2: TOMATO
tomato.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.cc2 <- tomato.cc2[1:12,]
tomato.cc2.bar <-ggplot(tomato.cc2, aes(x = Q4.na, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.cc2.bar

#CLIMATE SCALE 1: TOMATO CHANGE
tomato.cc1.change <- survey %>% 
  group_by(cc1) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.cc1.bar <-ggplot(na.omit(tomato.cc1.change), aes(x = tomato.change, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.cc1.bar


#CLIMATE SCALE 2: TOMATO CHANGE
tomato.cc2.change <- survey %>% 
  group_by(cc2) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.cc2.bar <-ggplot(na.omit(tomato.cc2.change), aes(x = tomato.change, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.cc2.bar





#CLIMATE SCALE 1: CATTLE
cattle.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.cc1 <- cattle.cc1[1:12,]
cattle.cc1.bar <-ggplot(cattle.cc1, aes(x = Q8.na, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.cc1.bar

#CLIMATE SCALE 2: TOMATO
cattle.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.cc2 <- cattle.cc2[1:12,]
cattle.cc2.bar <-ggplot(cattle.cc2, aes(x = Q8.na, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.cc2.bar

#CLIMATE SCALE 1: CATTLE CHANGE
cattle.cc1.change <- survey %>% 
  group_by(cc1) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.cc1.bar <-ggplot(na.omit(cattle.cc1.change), aes(x = cattle.change, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.cc1.bar


#CLIMATE SCALE 2: CATTLE CHANGE
cattle.cc2.change <- survey %>% 
  group_by(cc2) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.cc2.bar <-ggplot(na.omit(cattle.cc2.change), aes(x = cattle.change, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.cc2.bar





#CLIMATE SCALE 1: WHEAT
wheat.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.cc1 <- wheat.cc1[1:12,]
wheat.cc1.bar <-ggplot(wheat.cc1, aes(x = Q12.na, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.cc1.bar

#CLIMATE SCALE 2: WHEAT
wheat.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.cc2 <- wheat.cc2[1:12,]
wheat.cc2.bar <-ggplot(wheat.cc2, aes(x = Q12.na, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.cc2.bar

#CLIMATE SCALE 1: WHEAT CHANGE
wheat.cc1.change <- survey %>% 
  group_by(cc1) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.cc1.bar <-ggplot(na.omit(wheat.cc1.change), aes(x = wheat.change, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.cc1.bar


#CLIMATE SCALE 2: WHEAT CHANGE
wheat.cc2.change <- survey %>% 
  group_by(cc2) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.cc2.bar <-ggplot(na.omit(wheat.cc2.change), aes(x = wheat.change, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.cc2.bar









###GLOBALIZATION SCALE

#GLOBALIZATION SCALE 1: TOMATO
tomato.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.global1 <- tomato.global1[1:12,]
tomato.global1.bar <-ggplot(tomato.global1, aes(x = Q4.na, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.global1.bar



#GLOBAL SCALE 1: TOMATO CHANGE
tomato.global1.change <- survey %>% 
  group_by(global1) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.global1.bar <-ggplot(na.omit(tomato.global1.change), aes(x = tomato.change, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.global1.bar





#GLOBAL SCALE 1: CATTLE
cattle.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.global1 <- cattle.global1[1:12,]
cattle.global1.bar <-ggplot(cattle.global1, aes(x = Q8.na, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.global1.bar


#GLOBAL SCALE 1: CATTLE CHANGE
cattle.global1.change <- survey %>% 
  group_by(global1) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.global1.bar <-ggplot(na.omit(cattle.global1.change), aes(x = cattle.change, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.global1.bar







#CLIMATE SCALE 1: WHEAT
wheat.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.global1 <- wheat.global1[1:12,]
wheat.global1.bar <-ggplot(wheat.global1, aes(x = Q12.na, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.global1.bar


#CLIMATE SCALE 1: WHEAT CHANGE
wheat.global1.change <- survey %>% 
  group_by(global1) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.global1.bar <-ggplot(na.omit(wheat.global1.change), aes(x = wheat.change, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.global1.bar





###TECH SCALE

#TECH SCALE 1: TOMATO
tomato.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.tech1 <- tomato.tech1[1:12,]
tomato.tech1.bar <-ggplot(tomato.tech1, aes(x = Q4.na, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.tech1.bar

#TECH SCALE 2: TOMATO
tomato.tech2 <- survey %>% 
  group_by(tech2) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.tech2 <- tomato.tech2[1:12,]
tomato.tech2.bar <-ggplot(tomato.tech2, aes(x = Q4.na, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.tech2.bar

#TECH SCALE 1: TOMATO CHANGE
tomato.tech1.change <- survey %>% 
  group_by(tech1) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.tech1.bar <-ggplot(na.omit(tomato.tech1.change), aes(x = tomato.change, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.tech1.bar


#TECH SCALE 2: TOMATO CHANGE
tomato.tech2.change <- survey %>% 
  group_by(tech2) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.tech2.bar <-ggplot(na.omit(tomato.tech2.change), aes(x = tomato.change, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.tech2.bar





#TECH SCALE 1: CATTLE
cattle.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.tech1 <- cattle.tech1[1:12,]
cattle.tech1.bar <-ggplot(cattle.tech1, aes(x = Q8.na, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.tech1.bar

#TECH SCALE 2: Cattle
cattle.tech2 <- survey %>% 
  group_by(tech2) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.tech2 <- cattle.tech2[1:12,]
cattle.tech2.bar <-ggplot(cattle.tech2, aes(x = Q8.na, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.tech2.bar

#TECH SCALE 1: CATTLE CHANGE
cattle.tech1.change <- survey %>% 
  group_by(tech1) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.tech1.bar <-ggplot(na.omit(cattle.tech1.change), aes(x = cattle.change, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.tech1.bar


#TECH SCALE 2: CATTLE CHANGE
cattle.tech2.change <- survey %>% 
  group_by(tech2) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.tech2.bar <-ggplot(na.omit(cattle.tech2.change), aes(x = cattle.change, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.tech2.bar





#TECH SCALE 1: WHEAT
wheat.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.tech1 <- wheat.tech1[1:12,]
wheat.tech1.bar <-ggplot(wheat.tech1, aes(x = Q12.na, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.tech1.bar

#TECH SCALE 2: WHEAT
wheat.tech2 <- survey %>% 
  group_by(tech2) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.tech2 <- wheat.tech2[1:12,]
wheat.tech2.bar <-ggplot(wheat.tech2, aes(x = Q12.na, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.tech2.bar

#Tech SCALE 1: WHEAT CHANGE
wheat.tech1.change <- survey %>% 
  group_by(tech1) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.tech1.bar <-ggplot(na.omit(wheat.tech1.change), aes(x = wheat.change, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.tech1.bar


#CLIMATE SCALE 2: WHEAT CHANGE
wheat.tech2.change <- survey %>% 
  group_by(tech2) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.tech2.bar <-ggplot(na.omit(wheat.tech2.change), aes(x = wheat.change, y = percent, fill=tech2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.tech2.bar






###GREEN REVOLUTION SCALE

#GREEN SCALE 1: TOMATO
tomato.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.green1 <- tomato.green1[1:12,]
tomato.green1.bar <-ggplot(tomato.green1, aes(x = Q4.na, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.green1.bar

#GREEN SCALE 2: TOMATO
tomato.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q4.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.green2 <- tomato.green2[1:12,]
tomato.green2.bar <-ggplot(tomato.green2, aes(x = Q4.na, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Tomato case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.green2.bar

#GREEN SCALE 1: TOMATO CHANGE
tomato.green1.change <- survey %>% 
  group_by(green1) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.green1.bar <-ggplot(na.omit(tomato.green1.change), aes(x = tomato.change, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.green1.bar


#GREEN SCALE 2: TOMATO CHANGE
tomato.green2.change <- survey %>% 
  group_by(green2) %>%
  count(tomato.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

tomato.change.green2.bar <-ggplot(na.omit(tomato.green2.change), aes(x = tomato.change, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Tomato case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
tomato.change.green2.bar





#GREEN SCALE 1: CATTLE
cattle.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.green1 <- cattle.green1[1:12,]
cattle.green1.bar <-ggplot(cattle.green1, aes(x = Q8.na, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.green1.bar

#GREEN SCALE 2: TOMATO
cattle.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q8.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.green2 <- cattle.green2[1:12,]
cattle.green2.bar <-ggplot(cattle.green2, aes(x = Q8.na, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Cattle case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.green2.bar

#GREEN SCALE 1: CATTLE CHANGE
cattle.green1.change <- survey %>% 
  group_by(green1) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.green1.bar <-ggplot(na.omit(cattle.green1.change), aes(x = cattle.change, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.green1.bar


#Green SCALE 2: CATTLE CHANGE
cattle.green2.change <- survey %>% 
  group_by(green2) %>%
  count(cattle.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

cattle.change.green2.bar <-ggplot(na.omit(cattle.green2.change), aes(x = cattle.change, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Cattle case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
cattle.change.green2.bar





#GREEN SCALE 1: WHEAT
wheat.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.green1 <- wheat.green1[1:12,]
wheat.green1.bar <-ggplot(wheat.green1, aes(x = Q12.na, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.green1.bar

#GREEN SCALE 2: WHEAT
wheat.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q12.na) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.green2 <- wheat.green2[1:12,]
wheat.green2.bar <-ggplot(wheat.green2, aes(x = Q12.na, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(labels = c("Very comfortable",
                              "Comfortable",
                              "Neither comfortable nor uncomfortable",
                              "Uncomfortable",
                              "Very uncomfortable",
                              "Don't know/not sure"))+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Wheat case: Comfort PRE nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.green2.bar

#GREEN SCALE 1: WHEAT CHANGE
wheat.green1.change <- survey %>% 
  group_by(green1) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.green1.bar <-ggplot(na.omit(wheat.green1.change), aes(x = wheat.change, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.green1.bar


#GREEN SCALE 2: WHEAT CHANGE
wheat.green2.change <- survey %>% 
  group_by(green2) %>%
  count(wheat.change) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

wheat.change.green2.bar <-ggplot(na.omit(wheat.green2.change), aes(x = wheat.change, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree", "NA"))+
  labs(title="Wheat case: index of change, from pre to POST nudge", x="")+
  theme(axis.ticks.x = element_blank())
wheat.change.green2.bar









### DV2 data cleaning ########

#Recode
Q16A.na <- na_if(survey$Q16A, 6)
Q16B.na <- na_if(survey$Q16B, 6)
Q18A.na <- na_if(survey$Q18A, 6)
Q18B.na <- na_if(survey$Q18B, 6)

#REVERSE CODING ALL RESPONSES so that they are easier to understand/interpret!
#Q16A.na = recode(Q16A.na, "1=5; 2=4; 3=3; 4=2; 5=1; 7=7; 8=8") 
#Q16B.na = recode(Q16B.na, "1=5; 2=4; 3=3; 4=2; 5=1; 7=7; 8=8")
#Q18A.na = recode(Q18A.na, "1=5; 2=4; 3=3; 4=2; 5=1; 7=7; 8=8") 
#Q18B.na = recode(Q18B.na, "1=5; 2=4; 3=3; 4=2; 5=1; 7=7; 8=8") 


#Clean ordinal responses
Q16A_na2 <- dplyr::na_if(Q16A.na, 7)
Q16A_na3 <- dplyr::na_if(Q16A_na2, 8)

Q16B_na2 <- dplyr::na_if(Q16B.na, 7)
Q16B_na3 <- dplyr::na_if(Q16B_na2, 8)

Q18A_na2 <- dplyr::na_if(Q18A.na, 7)
Q18A_na3 <- dplyr::na_if(Q18A_na2, 8)

Q18B_na2 <- dplyr::na_if(Q18B.na, 7)
Q18B_na3 <- dplyr::na_if(Q18B_na2, 8)

# Combine split tradeoff questions into single columns to make cleaned variables
survey$Q16=rowSums(cbind(Q16A_na3,Q16B_na3),na.rm=TRUE)
survey$Q18=rowSums(cbind(Q18A_na3,Q18B_na3),na.rm=TRUE)
survey$Q16[survey$Q16==0] <-NA
survey$Q16
survey$Q18[survey$Q18==0] <-NA
survey$Q18

## Create nominal groups for responses for DV2
# Combine split tradeoff questions into single columns to make cleaned variables

#first, exclude don't know's
survey$Q16A_na <- dplyr::na_if(survey$Q16A, 6)
survey$Q16B_na <- dplyr::na_if(survey$Q16B, 6)
survey$Q18A_na <- dplyr::na_if(survey$Q18A, 6)
survey$Q18B_na <- dplyr::na_if(survey$Q18B, 6)


#then, compile columns for turning into categoricals
survey$Q16_cat=rowSums(cbind(survey$Q16A_na,survey$Q16B_na), na.rm=TRUE)
survey$Q18_cat=rowSums(cbind(survey$Q18A_na,survey$Q18B_na),na.rm=TRUE)


#Turn variables into categorical/factors
survey$Q16_dk=cut(survey$Q16_cat, br=c(0,2,3,5,8), labels= c("pro-GE", "neutral", "anti-GE", "opt-out"))
summary(survey$Q16_dk)

survey$Q18_dk=cut(survey$Q18_cat, br=c(0,2,3,5,8), labels= c("pro-GE", "neutral", "anti-GE", "opt-out"))
summary(survey$Q18_dk)


#binomial data prep

#first, exclude don't know's
Q16A__na <- dplyr::na_if(Q16A.na, 3)
Q16A__na2 <- dplyr::na_if(Q16A__na, 4)
Q16A__na3 <- dplyr::na_if(Q16A__na2, 5)

Q16B__na <- dplyr::na_if(Q16B.na, 3)
Q16B__na2 <- dplyr::na_if(Q16B__na, 4)
Q16B__na3 <- dplyr::na_if(Q16B__na2, 5)

Q18A__na <- dplyr::na_if(Q18A.na, 3)
Q18A__na2 <- dplyr::na_if(Q18A__na, 4)
Q18A__na3 <- dplyr::na_if(Q18A__na2, 5)

Q18B__na <- dplyr::na_if(Q18B.na, 3)
Q18B__na2 <- dplyr::na_if(Q18B__na, 4)
Q18B__na3 <- dplyr::na_if(Q18B__na2, 5)

#then, compile columns for turning into categoricals
Q16_bi=rowSums(cbind(Q16A__na3,Q16B__na3), na.rm=TRUE)
Q18_bi=rowSums(cbind(Q18A__na3,Q18B__na3), na.rm=TRUE)
Q16_bi[Q16_bi==0] <-NA
Q18_bi[Q18_bi==0] <-NA



#Split sample--prep for regressions
Q16_bi_split=Q16B__na3
Q18_bi_split=Q18B__na3



#Turn variables into binary variable
Q16_all_bi=cut(Q16_bi, br=c(0,2,8), labels= c("pro-GE", "opt-out"))
summary(Q16_all_bi)

Q18_all_bi=cut(Q18_bi, br=c(0,2,8), labels= c("pro-GE", "opt-out"))
summary(Q18_all_bi)




######DV2 ordered regression bar chart example #####
survey$Q16 <- as.factor(survey$Q16)
levels(survey$Q16) <- c("Large increase in gene-edited crops",
                             "Small increase in gene-edited crops",
                             "No change",
                             "Small decrease in gene-edited crops",
                             "Large decrease in gene-edited crops",
                             "Don't know/not sure")
pest.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.gender.bar <-ggplot(na.omit(pest.gender), aes(x = Q16, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Pesticide tradeoff, by gender", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")
pest.gender.bar


#SPLIT SAMPLE
survey$Q16B_na3 <- Q16B_na3
Q16B_na3 <- as.factor(Q16B_na3)
survey$Q18B_na3 <- Q18B_na3
Q18B_na3 <- as.factor(Q18B_na3)

levels(Q16B_na3) <- c("Large increase in gene-edited crops",
                      "Small increase in gene-edited crops",
                      "No change",
                      "Small decrease in gene-edited crops",
                      "Large decrease in gene-edited crops",
                      "Don't know/not sure")
pest.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q16B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.gender.bar <-ggplot(na.omit(pest.gender), aes(x = Q16B_na3, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Pesticide tradeoff, by gender", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")
pest.gender.bar




#Pesticide tradeoff--TRUST scale
pest.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.trust1.bar <-ggplot(na.omit(pest.trust1), 
                         aes(x = Q16, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")

pest.trust1.bar

#SPLIT SAMPLE
pest.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q16B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.trust1.bar <-ggplot(na.omit(pest.trust1), 
                         aes(x = Q16B_na3, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")
pest.trust1.bar


biodiv.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.trust1.bar <-ggplot(na.omit(biodiv.trust1), 
                         aes(x = Q16, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")
biodiv.trust1.bar

#SPLIT SAMPLE:

biodiv.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q18B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.trust1.bar <-ggplot(na.omit(biodiv.trust1), 
                           aes(x = Q18B_na3, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")

biodiv.trust1.bar



biodiv.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc2.bar <-ggplot(na.omit(biodiv.cc2), 
                           aes(x = Q16, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Climate ambivalence")+
  labs(title="", x="")

biodiv.cc2.bar

#Split sample

biodiv.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q18B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc2.bar <-ggplot(na.omit(biodiv.cc2), 
                        aes(x = Q18B_na3, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Climate ambivalence")+
  labs(title="", x="")

biodiv.cc2.bar



pest.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green1.bar <-ggplot(na.omit(pest.green1), aes(x = Q16, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
pest.green1.bar

#SPLIT SAMPLE
pest.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q16B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green1.bar <-ggplot(na.omit(pest.green1), aes(x = Q16B_na3, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
pest.green1.bar

#Pesticide tradeoff--green scale 2
pest.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green2.bar <-ggplot(na.omit(pest.green2), aes(x = Q16, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
pest.green2.bar



#Pesticide tradeoff--green scale 2
pest.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q16B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green2.bar <-ggplot(na.omit(pest.green2), aes(x = Q16B_na3, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
pest.green2.bar

#Biodiversity tradeoff--green scale
biodiv.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q18) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green1.bar <-ggplot(na.omit(biodiv.green1), aes(x = Q18, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
biodiv.green1.bar

#Biodiv tradeoff--green scale 2
biodiv.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q18) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green2.bar <-ggplot(na.omit(biodiv.green2), aes(x = Q18, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
biodiv.green2.bar

#SPLIT SAMPLE
#Pesticide tradeoff--green scale 2
biodiv.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q18B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green2.bar <-ggplot(na.omit(biodiv.green2), aes(x = Q18B_na3, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
biodiv.green2.bar




#Pesticide tradeoff--Tech scale
pest.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q16) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.tech1.bar <-ggplot(na.omit(pest.tech1), aes(x = Q16, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Techno-skepticism")+
  labs(title="", x="")
pest.tech1.bar


#Pesticide tradeoff--Tech scale
pest.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q16B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.tech1.bar <-ggplot(na.omit(pest.tech1), aes(x = Q16B_na3, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Techno-skepticism")+
  labs(title="", x="")
pest.tech1.bar




#Biodiversity tradeoff--tech scale
biodiv.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q18) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.tech1.bar <-ggplot(na.omit(biodiv.tech1), aes(x = Q18, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="Techno-skepticism")+
  labs(title="", x="")
biodiv.tech1.bar




#SPLIT SAMPLE
#Biodiversity tradeoff--tech scale
biodiv.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q18B_na3) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.tech1.bar <-ggplot(na.omit(biodiv.tech1), aes(x = Q18B_na3, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="Techno-skepticism")+
  labs(title="", x="")
biodiv.tech1.bar




#Biodiversity tradeoff--TRUST scale
biodiv.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q18) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.trust1.bar <-ggplot(na.omit(biodiv.trust1), 
                           aes(x = Q18, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")

biodiv.trust1.bar



###### DV2 bar charts #########

# Pesticide tradeoff, subset by gender
pest.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.gender.bar <-ggplot(na.omit(pest.gender), aes(x = Q16_dk, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Pesticide tradeoff, by gender", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")
pest.gender.bar


# Biodiversity tradeoff, subset by gender
biodiv.gender <- survey %>% 
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.gender.bar <-ggplot(na.omit(biodiv.gender), aes(x = Q18_dk, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  labs(title="Biodiversity tradeoff, by gender", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.gender.bar










# Pesticide tradeoff, subset by political
pest.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.political.bar <-ggplot(na.omit(pest.political), aes(x = Q16_dk, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Pesticide tradeoff, by political orientation", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

pest.political.bar


# Biodiversity tradeoff, subset by political
biodiv.political <- survey %>% 
  filter(!is.na(political)) %>%
  group_by(political) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.political.bar <-ggplot(na.omit(biodiv.political), aes(x = Q18_dk, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  labs(title="Biodiversity tradeoff, by political orientation", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.political.bar





# Pesticide tradeoff, subset by religious
pest.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.religious.bar <-ggplot(na.omit(pest.religious), aes(x = Q16_dk, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Pesticide tradeoff, by religiosity", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

pest.religious.bar


# Biodiversity tradeoff, subset by religious
biodiv.religious <- survey %>% 
  filter(!is.na(religious)) %>%
  group_by(religious) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.religious.bar <-ggplot(na.omit(biodiv.religious), aes(x = Q18_dk, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-religious", "religious"))+
  labs(title="Biodiversity tradeoff, by religiosity", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.religious.bar




# Pesticide tradeoff, subset by uni
pest.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.uni.bar <-ggplot(na.omit(pest.uni), aes(x = Q16_dk, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Pesticide tradeoff, by education", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

pest.uni.bar


# Biodiversity tradeoff, subset by uni
biodiv.uni <- survey %>% 
  filter(!is.na(uni)) %>%
  group_by(uni) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.uni.bar <-ggplot(na.omit(biodiv.uni), aes(x = Q18_dk, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("no uni", "uni"))+
  labs(title="Biodiversity tradeoff, by education", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.uni.bar


# Pesticide tradeoff, subset by age
pest.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.age.bar <-ggplot(na.omit(pest.age), aes(x = Q16_dk, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Pesticide tradeoff, by age", x="")+
  theme(axis.ticks.x = element_blank())
pest.age.bar


# Biodiversity tradeoff, subset by age
biodiv.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.age.bar <-ggplot(na.omit(biodiv.age), aes(x = Q18_dk, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Biodiversity tradeoff, by age", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.age.bar




# Pesticide tradeoff, subset by race
pest.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.race.bar <-ggplot(na.omit(pest.race), aes(x = Q16_dk, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Pesticide tradeoff, by race", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

pest.race.bar


# # Pesticide tradeoff, subset by Familiarity with GM
pest.gm <- survey %>% 
  filter(!is.na(familiar_gm)) %>%
  group_by(familiar_gm) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.gm.bar <-ggplot(na.omit(familiar_gm), aes(x = Q16_dk, y = percent, fill=familiar_gm)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Strongly disagree", "Disagree", "Neutral", "Agree", "Strongly agree"))+
  labs(title="Pesticide tradeoff, by familiarity w GM", x="")+
  theme(axis.ticks.x = element_blank())
pest.gm.bar


# Biodiversity tradeoff, subset by age
biodiv.age <- survey %>% 
  filter(!is.na(age)) %>%
  group_by(age) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.age.bar <-ggplot(na.omit(biodiv.age), aes(x = Q18_dk, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  labs(title="Biodiversity tradeoff, by age", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.age.bar


# Biodiversity tradeoff, subset by race
biodiv.race <- survey %>% 
  filter(!is.na(race)) %>%
  group_by(race) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.race.bar <-ggplot(na.omit(biodiv.race), aes(x = Q18_dk, y = percent, fill=race)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("non-white", "white"))+
  labs(title="Biodiversity tradeoff, age", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")

biodiv.race.bar





#Pesticide tradeoff--Climate scale
pest.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.cc1.bar <-ggplot(na.omit(pest.cc1), aes(x = Q16_dk, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Pesticide tradeoff, by climate scale 1", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels= c("disagree", "agree"))

pest.cc1.bar

#Pesticide tradeoff--climate scale 2
pest.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.cc2.bar <-ggplot(na.omit(pest.cc2), aes(x = Q16_dk, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Pesticide tradeoff, by climate scale 2", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

pest.cc2.bar


#Biodiversity tradeoff--Climate scale
biodiv.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc1.bar <-ggplot(na.omit(biodiv.cc1), aes(x = Q18_dk, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Biodiv tradeoff, by climate scale 1", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree, agree"))

biodiv.cc1.bar

#Pesticide tradeoff--climate scale 2
biodiv.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc2.bar <-ggplot(na.omit(biodiv.cc2), aes(x = Q18_dk, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Biodiv tradeoff, by climate scale 2", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

biodiv.cc2.bar







#Pesticide tradeoff--Globalization scale
pest.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.global1.bar <-ggplot(na.omit(pest.global1), aes(x = Q16_dk, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Pesticide tradeoff, by globalization scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

pest.global1.bar


#Biodiversity tradeoff--Globalization scale
biodiv.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.global1.bar <-ggplot(na.omit(biodiv.global1), aes(x = Q18_dk, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  labs(title="Biodiv tradeoff, by globalization scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))


biodiv.global1.bar






#Pesticide tradeoff--Tech scale
pest.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.tech1.bar <-ggplot(na.omit(pest.tech1), aes(x = Q16_dk, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Pesticide tradeoff, by tech scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))
pest.tech1.bar




#Biodiversity tradeoff--tech scale
biodiv.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.tech1.bar <-ggplot(na.omit(biodiv.tech1), aes(x = Q18_dk, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, tech scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

biodiv.tech1.bar







#Pesticide tradeoff--Green revolution scale
pest.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green1.bar <-ggplot(na.omit(pest.green1), aes(x = Q16_dk, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Pesticide tradeoff, by Green Revolution scale 1", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

pest.green1.bar

#Pesticide tradeoff--green scale 2
pest.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green2.bar <-ggplot(na.omit(pest.green2), aes(x = Q16_dk, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Pesticide tradeoff, by Green Revolution scale 2", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

pest.green2.bar


#Biodiversity tradeoff--green scale
biodiv.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green1.bar <-ggplot(na.omit(biodiv.green1), aes(x = Q18_dk, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, by Green Revolution scale 1", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

biodiv.green1.bar

#Pesticide tradeoff--green scale 2
biodiv.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green2.bar <-ggplot(na.omit(biodiv.green2), aes(x = Q18_dk, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, by Green Revolution scale 2", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

biodiv.green2.bar







#Pesticide tradeoff--Toxics scale
pest.toxics1 <- survey %>% 
  group_by(toxics1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.toxics1.bar <-ggplot(na.omit(pest.toxics1), aes(x = Q16_dk, y = percent, fill=toxics1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Pesticide tradeoff, by toxics scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

pest.toxics1.bar




#Biodiversity tradeoff--toxics scale
biodiv.toxics1 <- survey %>% 
  group_by(toxics1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.toxics1.bar <-ggplot(na.omit(biodiv.toxics1), aes(x = Q18_dk, y = percent, fill=toxics1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, by toxics scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"))

biodiv.toxics1.bar







#####DV2 bar charts COMBINED --Demographics#####


# Pesticide tradeoff, subset by gender
pest.gender.bar <-ggplot(na.omit(pest.gender), aes(x = Q16_dk, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
pest.gender.bar


# Biodiversity tradeoff, subset by gender
biodiv.gender.bar <-ggplot(na.omit(biodiv.gender), aes(x = Q18_dk, y = percent, fill=gender)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Female", "Male", "NA"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
biodiv.gender.bar


# Pesticide tradeoff, subset by political
pest.political.bar <-ggplot(na.omit(pest.political), aes(x = Q16_dk, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
pest.political.bar


# Biodiversity tradeoff, subset by political
biodiv.political.bar <-ggplot(na.omit(biodiv.political), aes(x = Q18_dk, y = percent, fill=political)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Conservative", "Liberal", "Moderate"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
biodiv.political.bar


# Pesticide tradeoff, subset by religious
pest.religious.bar <-ggplot(na.omit(pest.religious), aes(x = Q16_dk, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Non-religious", "Religious"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
pest.religious.bar


# Biodiversity tradeoff, subset by religious
biodiv.religious.bar <-ggplot(na.omit(biodiv.religious), aes(x = Q18_dk, y = percent, fill=religious)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Non-religious", "Religious"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
biodiv.religious.bar


# Pesticide tradeoff, subset by uni
pest.uni.bar <-ggplot(na.omit(pest.uni), aes(x = Q16_dk, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("No uni", "Uni"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
pest.uni.bar

# Biodiversity tradeoff, subset by uni
biodiv.uni.bar <-ggplot(na.omit(biodiv.uni), aes(x = Q18_dk, y = percent, fill=uni)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("No uni", "Uni"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
biodiv.uni.bar


# Pesticide tradeoff, subset by age
pest.age.bar <-ggplot(na.omit(pest.age), aes(x = Q16_dk, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
pest.age.bar


# Biodiversity tradeoff, subset by age
biodiv.age.bar <-ggplot(na.omit(biodiv.age), aes(x = Q18_dk, y = percent, fill=age)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("Middle", "Older", "Young"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2")+
  labs(title="", x="")
biodiv.age.bar

  

#All demographics plots together (pesticide tradeoff)
pest.dem.all <-ggarrange(pest.gender.bar, pest.age.bar, pest.uni.bar, 
                    pest.political.bar, pest.religious.bar,
                    labels = c("gender", "age", "education", 
                               "politics", "religiosity"),
                    ncol=3, nrow=2)
pest.dem.all

#All demographics plots together (biodiversity tradeoff)
biodiv.dem.all <-ggarrange(biodiv.gender.bar, biodiv.age.bar, biodiv.uni.bar, 
                         biodiv.political.bar, biodiv.religious.bar,
                         labels = c("gender", "age", "education", 
                                    "politics", "religiosity"),
                         ncol=3, nrow=2)
biodiv.dem.all


#####DV2 bar charts COMBINED --Scales#####

#Pesticide tradeoff--Climate scale
survey$Q16_all_bi <- Q16_all_bi
survey$Q18_all_bi <- Q18_all_bi


pest.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.cc1.bar <-ggplot(na.omit(pest.cc1), aes(x = Q16_all_bi, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels= c("disagree", "agree"), 
                    name = "Climate urgency")+
  labs(title="", x="")
pest.cc1.bar

#Pesticide tradeoff--climate scale 2
pest.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.cc2.bar <-ggplot(na.omit(pest.cc2), aes(x = Q16_all_bi, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  scale_fill_discrete(labels=c("disagree", "agree"))+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Climate ambivalence")+
  labs(title="", x="")
pest.cc2.bar


#Biodiversity tradeoff--Climate scale
biodiv.cc1 <- survey %>% 
  group_by(cc1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc1.bar <-ggplot(na.omit(biodiv.cc1), aes(x = Q18_dk, y = percent, fill=cc1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Climate urgency")+
  labs(title="", x="")

biodiv.cc1.bar

#Pesticide tradeoff--climate scale 2
biodiv.cc2 <- survey %>% 
  group_by(cc2) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.cc2.bar <-ggplot(na.omit(biodiv.cc2), aes(x = Q18_dk, y = percent, fill=cc2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, by climate scale 2", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Climate ambivalence")+
  labs(title="", x="")
biodiv.cc2.bar


#CLIMATE COMBINED (pesticide)
pest.cc.all <-ggarrange(pest.cc1.bar, pest.cc2.bar,
                        labels = c("Climate urgency", "Climate ambivalence"),
                        ncol=2, nrow=1)
pest.cc.all

#CLIMATE COMBINED (biodiversity)
biodiv.cc.all <-ggarrange(biodiv.cc1.bar, biodiv.cc2.bar,
                        labels = c("Climate urgency", "Climate ambivalence"),
                        ncol=2, nrow=1)
biodiv.cc.all



#Pesticide tradeoff--Globalization scale
pest.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q16_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.global1.bar <-ggplot(na.omit(pest.global1), 
                          aes(x = Q16_dk, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Corporate skepticism")+
  labs(title="", x="")

pest.global1.bar


#Biodiversity tradeoff--Globalization scale
biodiv.global1 <- survey %>% 
  group_by(global1) %>%
  count(Q18_dk) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.global1.bar <-ggplot(na.omit(biodiv.global1), 
                            aes(x = Q18_dk, y = percent, fill=global1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Corporate skepticism")+
  labs(title="", x="")
biodiv.global1.bar


#Pesticide tradeoff--TRUST scale
pest.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.trust1.bar <-ggplot(na.omit(pest.trust1), 
                         aes(x = Q16_all_bi, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")

pest.trust1.bar




#Biodiversity tradeoff--TRUST scale
biodiv.trust1 <- survey %>% 
  group_by(trust1) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.trust1.bar <-ggplot(na.omit(biodiv.trust1), 
                         aes(x = Q18_all_bi, y = percent, fill=trust1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), name="Trust")+
  labs(title="", x="")

biodiv.trust1.bar


#COMBINED GLOBAL AND TRUST
#(pesticide)
pest.gt.all <-ggarrange(pest.global1.bar, pest.trust1.bar,
                        labels = c("Corporate skepticism", "Trust"),
                        ncol=2, nrow=1)
pest.gt.all

# COMBINED GLOBAL AND TRUST (biodiversity)
biodiv.gt.all <-ggarrange(biodiv.global1.bar, biodiv.trust1.bar,
                          labels = c("Corporate skepticism", "Trust"),
                          ncol=2, nrow=1)
biodiv.gt.all



#Pesticide tradeoff--Tech scale
pest.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.tech1.bar <-ggplot(na.omit(pest.tech1), aes(x = Q16_all_bi, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Techno-skepticism")+
  labs(title="", x="")
pest.tech1.bar




#Biodiversity tradeoff--tech scale
biodiv.tech1 <- survey %>% 
  group_by(tech1) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.tech1.bar <-ggplot(na.omit(biodiv.tech1), aes(x = Q18_all_bi, y = percent, fill=tech1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="Techno-skepticism")+
  labs(title="", x="")
biodiv.tech1.bar



#Pesticide tradeoff--Toxics scale
pest.toxics1 <- survey %>% 
  group_by(toxics1) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.toxics1.bar <-ggplot(na.omit(pest.toxics1), aes(x = Q16_all_bi, y = percent, fill=toxics1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"), 
                    name="Bodily resilience")+
  labs(title="", x="")
pest.toxics1.bar




#Biodiversity tradeoff--toxics scale
biodiv.toxics1 <- survey %>% 
  group_by(toxics1) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.toxics1.bar <-ggplot(na.omit(biodiv.toxics1), aes(x = Q18_all_bi, y = percent, fill=toxics1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  labs(title="Biodiv tradeoff, by toxics scale", x="")+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="Bodily resilience")+
  labs(title="", x="")
biodiv.toxics1.bar





#COMBINED Tech and Toxics (Pest)
pest.tt.all <-ggarrange(pest.tech1.bar, pest.toxics1.bar,
                        labels = c("Techno-skepticism", 
                                   "Bodily resilience"),
                        ncol=2, nrow=1)
pest.tt.all

#COMBINED Tech and toxics (biodiversity)
biodiv.tt.all <-ggarrange(biodiv.tech1.bar, biodiv.toxics1.bar,
                          labels = c("Techno-skepticism", 
                                     "Bodily resilience"),
                          ncol=2, nrow=1)
biodiv.tt.all




#Pesticide tradeoff--Green revolution scale
pest.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green1.bar <-ggplot(na.omit(pest.green1), aes(x = Q16_all_bi, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
pest.green1.bar

#Pesticide tradeoff--green scale 2
pest.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green2.bar <-ggplot(na.omit(pest.green2), aes(x = Q16_all_bi, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
pest.green2.bar


#Biodiversity tradeoff--green scale
biodiv.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green1.bar <-ggplot(na.omit(biodiv.green1), aes(x = Q18_all_bi, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
biodiv.green1.bar

#Pesticide tradeoff--green scale 2
biodiv.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green2.bar <-ggplot(na.omit(biodiv.green2), aes(x = Q18_all_bi, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
biodiv.green2.bar


#COMBINED Green Rev --pest
pest.gr.all <-ggarrange(pest.green1.bar, pest.green2.bar,
                        labels = c("GR skepticism", 
                                   "GR optimism"),
                        ncol=2, nrow=1)
pest.gr.all

#COMBINED Green Rev --biodiv
biodiv.gr.all <-ggarrange(biodiv.tech1.bar, biodiv.toxics1.bar,
                          labels = c("GR skepticism", 
                                     "GR optimism"),
                          ncol=2, nrow=1)
biodiv.gr.all


### Cross tab of DV questions #####

dv <- table(comfort, survey$Q16, survey$Q18)

library(gmodels)
CrossTable(survey$Q12.na, survey$Q16)



########## Some wilcoxon tests ####################

# tomato test
wilcox.test(Q4.na, Q6.na, paired = TRUE)

# cattle test
wilcox.test(Q8.na, Q10.na, paired = TRUE)

# wheat test
wilcox.test(Q12.na, Q14.na, paired = TRUE)

#comparison of tomatoes and cattle
wilcox.test(Q4.na, Q8.na, paired=TRUE)

#comparison of tomatoes and wheat
wilcox.test(Q4.na, Q12.na, paired=TRUE)

#comparison of cattle and wheat
wilcox.test(Q8.na, Q12.na, paired=TRUE)


# Wilcoxon test comparing preferences between pesticide and biodiversity tradeoffs
wilcox.test(survey$Q16, survey$Q18, paired=TRUE)



##### DV1 regression on change index ####

tomato.model2 = glm(survey$tomato.change ~ gender + climate.fa$scores + tech.fa$scores 
                  + green.fa$scores + global.fa$scores, 
                  data=regdata, family="binomial")
summary(tomato.model2)

cattle.model2 = glm(survey$cattle.change ~ gender + climate.fa$scores + tech.fa$scores 
                    + green.fa$scores + global.fa$scores, 
                    data=regdata, family="binomial")
summary(cattle.model2)

wheat.model2 = glm(survey$wheat.change ~ gender + climate.fa$scores + tech.fa$scores 
                    + green.fa$scores + global.fa$scores, 
                    data=regdata, family="binomial")
summary(wheat.model2)






######Correlation matrices####
survey2 <-cbind(age_full, gender_full, political_full, religious_full, 
                uni_full, race_full, income,
                trustmean, familiar_ge, globalmean, 
                climate2mean, green1mean, green2mean)

library(expss)
colnames(survey2) <- c("Age","Gender", "Political (conservatism)", 
                       "Religiosity", "Formal education", "Race/ethnicity",
                       "Income", "Trust", "Familiarity with gene editing",
                       "Corporate criticism", "Climate ambivalence",
                       "Green Revolution criticism", "Green Revolution optimism")

survey2 <- char2numeric(survey2)
library(corrplot)
corrplot(cor(survey2, method="pearson", use="complete.obs"), 
         tl.col = "black", tl.srt = 45, method = "number",  order = "FPC")

correlation <- corrplot(cor(survey2, method="pearson", use="complete.obs"), 
         tl.col = "black", tl.srt = 45, method = "circle",  order = "FPC")

library(corrr)
survey2 %>% 
  correlate() %>% 
  network_plot()

#correlations of gr skepticism with everything else
s2 <- correlate(survey2)
s2 %>%
  focus(green1mean) %>%
  mutate(rowname = reorder(term, green1mean)) %>%
  ggplot(aes(term, green1mean)) +
  geom_col() + coord_flip()

#correlations of gr optimism with everything else
s2 %>%
  focus(green2mean) %>%
  mutate(rowname = reorder(term, green2mean)) %>%
  ggplot(aes(term, green2mean)) +
  geom_col() + coord_flip()

gr <-cbind(green1mean, green2mean)
library(corrplot)
corrplot(cor(gr, method="pearson", use="complete.obs"))

gr_2 <- survey$Q25_2.na
gr_3 <- survey$Q25_3.na
gr_4 <- survey$Q25_4.na
gr_5 <- survey$Q25_5.na
gr_6 <- survey$Q25_6.na
gr_7 <- survey$Q25_7.na
gr_8 <- survey$Q25_8.na
gr2 <-cbind(gr_2, gr_3, gr_4, gr_5, gr_6, gr_7, gr_8)
library(corrplot)
corrplot(cor(gr2, method="pearson", use="complete.obs"))

####Ordinal regression on DV2####


###Regressions on ordinal variables for TRADEOFF questions
pestdv <- survey$Q16 #Q16B_na3 # 
pestdv <- factor(pestdv, ordered=TRUE)
biodivdv <- survey$Q18 #Q18B_na3 #
biodivdv <- factor(biodivdv, ordered=TRUE)
pestdata <- data.frame(pestdv, gender_full, uni_full, age_full, 
                       religious_full, political_full, income, trustmean, familiar_ge, 
                       globalmean, climate2mean, green1mean, green2mean)
biodivdata <- data.frame(biodivdv, gender_full, uni_full, age_full, 
                         religious_full, political_full, income, trustmean, familiar_ge,
                         globalmean, climate2mean, green1mean, 
                         green2mean)

#remove missing data so that step function runs
pestdata<-na.omit(pestdata)
biodivdata<-na.omit(biodivdata)

##Pesticide tradeoff
pest_step <- glm(formula = pestdv ~ trustmean + familiar_ge +
                   globalmean + 
                   climate2mean + green1mean + green2mean +
                   gender_full + uni_full + age_full +  
                   religious_full + political_full + income,
                 family = binomial(link="logit"), data = pestdata)
summary(pest_step)
pest_or <- odds.ratio(pest_step, level=0.95)
pest_or
write.csv(pest_or, file="pest_ordinal_or.csv")

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(pest_step)

##Biodiversity tradeoff 
biodiv_step <- glm(formula = biodivdv ~ trustmean + familiar_ge + 
                     globalmean + 
                     climate2mean + green1mean + green2mean + 
                     gender_full + uni_full + age_full +  
                     religious_full + political_full + income,
                   family = binomial(link="logit"), data = biodivdata)
summary(biodiv_step)

biodiv_or <- odds.ratio(biodiv_step, level=0.95)
biodiv_or
write.csv(biodiv_or, file="biodiv_ordinal_or.csv")

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(biodiv_step)



### Plots

pest_ordered <- plot_model(pest_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey") +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 12, base_family = "")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 0.3, 1, 0.7), "cm"))

biodiv_ordered <- plot_model(biodiv_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey") +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 12, base_family = "")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 1, 1, 0.3), "cm"))+
  theme(axis.title.y=element_blank())

library(cowplot) 
ordered_combined <- plot_grid(pest_ordered, 
          biodiv_ordered, 
          labels = c('Pesticide tradeoff', 'Biodiversity tradeoff'), 
          label_size = 12)
save_plot("./ordered combined.png", ordered_combined, base_height=8, ncol=1, nrow=1)




#######Regressions on categoricals (Multinomial logistic)####
pestdv_dk <- survey$Q16_dk
biodivdv_dk <- survey$Q18_dk
pestdata_dk <- data.frame(pestdv_dk, gender_full, uni_full, race_full, age_full, 
                          religious_full, political_full, income, familiar_gm, 
                          familiar_ge, familiar_gd, similar, trustmean, toxicmean, 
                          techmean, globalmean, climate1mean, climate2mean, green1mean, green2mean)
biodivdata_dk <- data.frame(biodivdv_dk, gender_full, uni_full, race_full, age_full, 
                            religious_full, political_full, income, familiar_gm, 
                            familiar_ge, familiar_gd, similar, trustmean, toxicmean, 
                            techmean, globalmean, climate1mean, climate2mean, green1mean, green2mean)

#remove missing data so that step function runs
pestdata_dk<-na.omit(pestdata_dk)
biodivdata_dk<-na.omit(biodivdata_dk)

##Pesticide tradeoff
#backward selection (pesticide data)
require(nnet)
pest_all_dk <- multinom(pestdv_dk ~., data=pestdata_dk)
step(pest_all_dk, direction = "backward")

#forward selection (pesticide)
pest_start_dk <- multinom(pestdv_dk ~ 1, data=pestdata_dk)
step(pest_start_dk, direction="forward", scope=formula(pest_all_dk))

#stepwise selection (pesticide)
pest_start_dk <- multinom(pestdv_dk ~ 1, data=pestdata_dk)
step(pest_start_dk, direction="both", scope=formula(pest_all_dk))

#based on stepwise selection, what is best model?
pest_step_dk <- multinom(formula = pestdv_dk ~ age_full + political_full + gender_full,
                 data = pestdata_dk)
summary(pest_step_dk)

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(pest_step_dk)

pest_multi_or <- odds.ratio(pest_step_dk, level=0.95)
write.csv(pest_multi_or, file="pest multinomial or.csv")




### p value
z<- summary(pest_step_dk)$coefficients/summary(pest_step_dk)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1))*2
p
coef(pest_step_dk)

# extract coefficients and exponentiate
exp(coef(pest_step_dk))


library(broom)

pest.df <- tidy(pest_step_dk, conf.int = T)  # Convert model to dataframe for easy manipulation
pest.df
write.csv(pest.df, file="pest multinomial.csv")

z <- pest.df %>% mutate(or = exp(estimate),  # Odds ratio/gradient
                         var.diag = diag(vcov(pest_step_dk)),  # Variance of each coefficient
                         OR_adjusted = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 
z



###Odds ratio adjusted (need to use this)
z <- pest.df %>% mutate(or = exp(estimate),  # Odds ratio/gradient
                          var.diag = diag(vcov(pest_step_dk)),  # Variance of each coefficient
                          OR_adjusted = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

print(as_tibble(z), n=24)


#Confidence intervals
get.or.se <- function(biodiv_step_dk) {
  broom::tidy(biodiv_step_dk) %>% 
    mutate(or = exp(estimate),
           var.diag = diag(vcov(biodiv_step_dk)),
           or.se = sqrt(or^2 * var.diag)) %>%
    select(or.se) %>% unlist %>% unname
}

get.or.se(biodiv_step_dk)



##Biodiversity tradeoff
#backward selection (biodiv data)
biodiv_all_dk <- multinom(biodivdv_dk ~ ., data=biodivdata_dk)
step(biodiv_all_dk, direction = "backward")


#forward selection (biodiv)
biodiv_start_dk <- multinom(biodivdv_dk ~ 1, data=biodivdata_dk)
step(biodiv_start_dk, direction="forward", scope=formula(biodiv_all_dk))

#stepwise selection (biodiv)
biodiv_start_dk <- multinom(biodivdv_dk ~ 1, data=biodivdata_dk)
step(biodiv_start_dk, direction="both", scope=formula(biodiv_all_dk))

#based on stepwise selection, what is best model?
biodiv_step_dk <- multinom(formula = biodivdv_dk ~ trustmean + techmean +
                             green2mean + green1mean + age_full
                           + gender_full, 
                  data = biodivdata_dk)
summary(biodiv_step_dk)

biodiv_multi_or <- odds.ratio(biodiv_step_dk, level=0.95)
write.csv(biodiv_multi_or, file="biodiv multinomial or.csv")


#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(biodiv_step_dk)


###NOTE: VIFs are really high. SO, need to remove some items from this model...
##namely demographic stuff that was NOT in the pest model, e.g. political views, income
#that wasn't enough, so also removed familiarity with gm
##the VIFs were STILL high! what now???




### p value, standard errors
z<- summary(biodiv_step_dk)$coefficients/summary(pest_step_dk)$standard.errors
z
p <- (1 - pnorm(abs(z), 0, 1))*2
p

#What is this?
coef(biodiv_step_dk)
exp(coef(biodiv_step_dk))


library(broom)

biodiv.df <- tidy(biodiv_step_dk, conf.int = T)  # Convert model to dataframe for easy manipulation
biodiv.df
write.csv(biodiv.df, file="biodiv multinomial.csv")

z <- biodiv.df %>% mutate(or = exp(estimate),  # Odds ratio/gradient
                        var.diag = diag(vcov(biodiv_step_dk)),  # Variance of each coefficient
                        OR_adjusted = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 
z

###Odds ratio adjusted (need to use this)
z <- biodiv.df %>% mutate(or = exp(estimate),  # Odds ratio/gradient
                         var.diag = diag(vcov(biodiv_step_dk)),  # Variance of each coefficient
                         OR_adjusted = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 

print(as_tibble(z), n=24)


#Confidence intervals
get.or.se <- function(biodiv_step_dk) {
  broom::tidy(biodiv_step_dk) %>% 
    mutate(or = exp(estimate),
           var.diag = diag(vcov(biodiv_step_dk)),
           or.se = sqrt(or^2 * var.diag)) %>%
    select(or.se) %>% unlist %>% unname
}

get.or.se(biodiv_step_dk)




####### Binomial regressions #########

##RUN THIS OR NOT?
#survey[!survey$Q25_1.na=='1',]
#survey[!survey$Q25_1.na=='2',]


###examine plots of those who are both in agreement with both GR1 and GR2
#survey[gr1mean>'3',]
#survey[gr2mean>'3',]

pestdv_bi <- Q16_all_bi #Q16_bi_split #
pestdv_bi <- factor(pestdv_bi, ordered=TRUE)
biodivdv_bi <- Q18_all_bi #Q18_bi_split #
biodivdv_bi <- factor(biodivdv_bi, ordered=TRUE)
pestdata_bi <- data.frame(pestdv_bi, gender_full, uni_full, age_full, 
                          religious_full, political_full, income, trustmean, familiar_ge,
                          globalmean, climate1mean, climate2mean, green1mean, 
                          green2mean)
biodivdata_bi <- data.frame(biodivdv_bi, gender_full, uni_full, age_full, 
                            religious_full, political_full, income, trustmean, familiar_ge,
                            globalmean, climate1mean, climate2mean, green1mean, 
                            green2mean)



#remove missing data so that step function runs
pestdata_bi<-na.omit(pestdata_bi)
biodivdata_bi<-na.omit(biodivdata_bi)


##Pesticide tradeoff
pest_step_bi <- glm(pestdv_bi ~ trustmean +
                      familiar_ge + globalmean +
                      climate2mean +green1mean +green2mean+
                      gender_full + uni_full + age_full + 
                      religious_full + political_full +income,
                         family = binomial, data = pestdata_bi)
summary(pest_step_bi)

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(pest_step_bi)

pest_bi_or <- odds.ratio(pest_step_bi, level=0.95)
pest_bi_or 
write.csv(pest_bi_or, file="pest binomial or.csv")



##Biodiv tradeoff
biodiv_step_bi <- glm(biodivdv_bi ~ trustmean +
                        familiar_ge + globalmean +
                        climate2mean +green1mean +green2mean + 
                        gender_full + uni_full + age_full + 
                        religious_full + political_full +income,
                    family = binomial, data = biodivdata_bi)
summary(biodiv_step_bi)


#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(biodiv_step_bi)

biodiv_bi_or <- odds.ratio(biodiv_step_bi, level=0.95)
biodiv_bi_or
write.csv(biodiv_bi_or, file="biodiv binomial or.csv")





pest_bi <- plot_model(pest_step_bi, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
                           vline.color = "grey") +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 12, base_family = "Times") +
  theme(axis.title.x=element_blank())

biodiv_bi <- plot_model(biodiv_step_bi, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
                             vline.color = "grey") +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 12, base_family = "Times") +
  theme(axis.title.x=element_blank())

library(cowplot)
theme_set(theme_cowplot(font_size=12, font_family = "Times"))
bi_combined <- plot_grid(pest_bi, biodiv_bi, 
                         labels = c('Pesticide tradeoff', 'Biodiversity tradeoff'), 
                         label_size = 12) +
  theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())
save_plot("./bi combined.png", ordered_combined, base_height=8, ncol=1, nrow=1)





####Regressions on first DVs####

#Ordered logistic regressions
tomatodv <- survey$Q4.na
tomatodv <- factor(tomatodv, ordered=TRUE)
cattledv <- survey$Q8.na
cattledv <- factor(cattledv, ordered=TRUE)
wheatdv <- survey$Q12.na
wheatdv <- factor(wheatdv, ordered=TRUE)
tomatodata <- data.frame(tomatodv, trustmean, 
                         familiar_ge, globalmean, 
                         climate2mean, green1mean, green2mean,
                         gender_full, uni_full, age_full, 
                         religious_full, political_full, income)
cattledata <- data.frame(cattledv, trustmean, 
                         familiar_ge, globalmean, 
                         climate2mean, green1mean, green2mean,
                         gender_full, uni_full, age_full, 
                         religious_full, political_full, income)
wheatdata <- data.frame(wheatdv, trustmean, 
                        familiar_ge, globalmean, 
                        climate2mean, green1mean, green2mean,
                        gender_full, uni_full, age_full, 
                        religious_full, political_full, income)

#remove missing data so that step function runs
tomatodata <- na.omit(tomatodata)
cattledata <- na.omit(cattledata)
wheatdata <- na.omit(wheatdata)

##Tomato dv
tomato_step <- glm(formula = tomatodv ~ trustmean + familiar_ge + globalmean + 
                     climate2mean + green1mean + green2mean + gender_full + uni_full + 
                     age_full + religious_full + political_full + income,
                   data = tomatodata, binomial(link="logit"))
summary(tomato_step)

tomato_or <- odds.ratio(tomato_step, level=0.95)
tomato_or
write.csv(tomato_or, file="tomato_ordinal_or.csv")

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(tomato_step)


##cattle dv
cattle_step <- glm(formula = cattledv ~ trustmean + familiar_ge + 
                     globalmean + climate2mean + green1mean + green2mean + 
                     gender_full + uni_full + age_full + 
                     religious_full + political_full + income,
                   data = cattledata, binomial(link="logit"))
summary(cattle_step)
cattle_or <- odds.ratio(cattle_step, level=0.95)
cattle_or
write.csv(cattle_or, file="cattle_ordinal_or.csv")

vif(cattle_step)

##wheat dv
wheat_step <- glm(formula = wheatdv ~  
                     trustmean + familiar_ge + globalmean + 
                     climate2mean + green1mean + green2mean + gender_full + uni_full + 
                     age_full + religious_full + political_full + income,
                  data = wheatdata, binomial(link="logit"))
summary(wheat_step)

wheat_or <- odds.ratio(wheat_step, level=0.95)
wheat_or
write.csv(wheat_or, file="wheat_ordinal_or.csv")

vif(wheat_step)


tomato_model <- plot_model(tomato_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey") +ylim(0.4, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())
  

cattle_model <- plot_model(cattle_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey") +ylim(0.4, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())

wheat_model <- plot_model(wheat_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey") +ylim(0.4, 2) +
  scale_x_discrete(labels=list(
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    green1mean = "Green Revolution criticism",
    green2mean = "Green Revolution optimism",
    gender_full2 = "Gender (female)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())
  



cases_combined <- plot_grid(tomato_model, cattle_model, wheat_model, 
                         labels = c('Tomato', 'Cattle', 'Wheat'), 
                         label_size = 12, ncol=1)
save_plot("./cases combined.png", cases_combined, base_height=8)









#######Graphs on Green revolution scale items... #######

survey$gr2=cut(survey$Q25_2.na, br=c(0,2,3,5), labels= c("non-optimist", "neutral", "optimist"))
survey$gr2 <-as.factor(survey$gr2)
summary(survey$gr2)

Q25_2.na
Q25_3.na
Q25_4.na
Q25_5.na
Q25_6.na
Q25_7.na
Q25_8.na

#Pesticide tradeoff--Green revolution scale
pest.green1 <- survey %>% 
  group_by(gr2) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green1.bar <-ggplot(na.omit(pest.green1), aes(x = Q16_all_bi, y = percent, fill=gr2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
pest.green1.bar


#Pesticide tradeoff--green scale 2
pest.green2 <- survey %>% 
  group_by(Q25_2.na) %>%
  count(Q16_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

pest.green2.bar <-ggplot(na.omit(pest.green2), aes(x = Q16_all_bi, y = percent, fill=Q25_2.na)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
pest.green2.bar


#Biodiversity tradeoff--green scale
biodiv.green1 <- survey %>% 
  group_by(green1) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green1.bar <-ggplot(na.omit(biodiv.green1), aes(x = Q18_all_bi, y = percent, fill=green1)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR skepticism")+
  labs(title="", x="")
biodiv.green1.bar

#Pesticide tradeoff--green scale 2
biodiv.green2 <- survey %>% 
  group_by(green2) %>%
  count(Q18_all_bi) %>%
  mutate(percent = n / nrow(survey)*100) ### HOW DO I GET THIS TO BE as % per gender? not over all?

biodiv.green2.bar <-ggplot(na.omit(biodiv.green2), aes(x = Q18_all_bi, y = percent, fill=green2)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  theme_minimal()+
  coord_flip()+
  theme(axis.ticks.x = element_blank())+
  scale_fill_brewer(palette = "Dark2", labels=c("disagree", "agree"),
                    name="GR optimism")+
  labs(title="", x="")
biodiv.green2.bar





