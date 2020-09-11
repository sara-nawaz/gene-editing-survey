library(ggplot2)
library(tidyr)
library(tidyverse)
library(psych) # psychometrics
library(MVN) # multivariate normality # do i use this?
library(sjPlot)
library(fastDummies) # do i use this?
library(glmnet)
library(caret)
library(cowplot)


setwd("~/Desktop/Stuff/PhD/Research projects/Genome BC Perceptions of Gene Editing in Agriculture/SURVEY ANALYSIS/Survey analysis")
survey <- read.csv('Perceptions of gene editing in agriculture_August 3 2020 EXTRA CLEAN numeric.csv', sep= ',', stringsAsFactors=FALSE, header=TRUE)

###### Pre-questions #######

# Familiarity with GM, gene-editing, & drives
Q1_1.na <-na_if(survey$Q1_1, 6)
Q1_2.na <-na_if(survey$Q1_2, 6)
Q1_3.na <-na_if(survey$Q1_3, 6)
# combined familiarity
familiar <- data.frame(Q1_1.na, Q1_2.na, Q1_3.na)

# graphs of familiarity
familiar1<- ggplot(data=survey, aes(Q1_1, fill=as.factor(Q1_1)))+
geom_bar()+
  coord_cartesian(ylim = c(0, 1000))+
  labs(title="Familiarity with GM", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Familiarity", 
                      labels = c("Strongly disagree",
                                 "Disagree",
                                 "Neutral",
                                 "Agree",
                                 "Strongly agree",
                                 "Don't know/not sure"))
familiar1

tomato.comfort <-as.factor(survey$Q4)
levels(tomato.comfort) <-c("1" ="Very comfortable", 
                      "2" = "Comfortable",
                      "3" = "Neither comfortable nor uncomfortable",
                      "4" = "Uncomfortable",
                      "5" = "Very uncomfortable",
                      "6" = "Don't know/not sure")
familiar1 + facet_grid(tomato.comfort)


familiar2<- ggplot(data=survey, aes(Q1_2, fill=as.factor(Q1_2)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 1000))+
  labs(title="Familiarity with gene editing", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Familiarity", 
                      labels = c("Strongly disagree",
                                 "Disagree",
                                 "Neutral",
                                 "Agree",
                                 "Strongly agree",
                                 "Don't know/not sure"))
familiar2
familiar2 + facet_grid(tomato.comfort)

familiar3<- ggplot(data=survey, aes(Q1_3, fill=as.factor(Q1_3)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 1000))+
  labs(title="Familiarity with gene drives", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Familiarity", 
                      labels = c("Strongly disagree",
                                 "Disagree",
                                 "Neutral",
                                 "Agree",
                                 "Strongly agree",
                                 "Don't know/not sure"))
familiar3
familiar3 + facet_grid(tomato.comfort)


familiarity <- plot_grid(familiar1, familiar2, familiar3, 
                         align="hv", labels = NULL)
familiarity


# Perceived similarity between GM, gene editing, gene drives
Q2.na <- na_if(survey$Q2, 6)

###### DV1 data cleaning #####
# Combine and label variables, convert "don't know/not sure" responses to NA's


# Tomatoes
# pre-nudge
Q4.na <- na_if(survey$Q4, 6)
# post-nudge
Q6.na <- na_if(survey$Q6, 6)
# change from pre- to post-nudge
tomato.change <- Q6.na - Q4.na

# Cattle
# pre-nudge
Q8.na <- na_if(survey$Q8, 6)
# post-nudge
Q10.na <- na_if(survey$Q10, 6)
#change from pre- to post-nudge
cattle.change <- Q10.na - Q8.na

# Wheat
# pre-nudge
Q12.na <- na_if(survey$Q12, 6)
# post-nudge
Q14.na <- na_if(survey$Q14, 6)
#change from pre- to post-nudge
wheat.change <- Q14.na - Q12.na



# Histograms for cases
#### DV1 bar charts #############

# Histogram of tomato comfort
tomato.pre <- ggplot(data=survey, aes(Q4, fill=as.factor(Q4)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Tomato case: Comfort PRE nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))
tomato.pre

# Print tomato histogram POST nudge
tomato.post <- ggplot(data=survey, aes(Q6, fill=as.factor(Q6)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Tomato case: Comfort POST nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
      axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))
tomato.post

tomato.prepost <- plot_grid(tomato.pre, tomato.post, nrow=1,
                    hjust = -1, labels = NULL,
                    label_size=12)
tomato.prepost

# Tomato faceting 

# POLITICAL faceting
political <-as.factor(survey$Q29)
levels(political) <-c("1" ="Very liberal", 
                      "2" = "Liberal",
                      "3" = "Moderate",
                      "4" = "Conservative",
                      "5" = "Very conservative")
tomato.pre + facet_grid(political)
tomato.post + facet_grid(political)

# Gender faceting
gender <- as.factor(survey$Q26_clean)
levels(gender) <- c("1" = "Female",
                    "2" = "Male",
                    "3" = "Nonbinary",
                    "4" = "Transitioning")
tomato.pre + facet_grid(gender)
tomato.post + facet_grid(gender)

# Race faceting
race <-as.factor(survey$Q30_clean)
levels(race) <-c("1"="East Asian",
                 "2"="Black",
                 "3"= "American Indian, First Nation, or Pacific Islander",
                 "4"="Hispanic",
                 "5"="Non-white mixed race",
                 "6"="White",
                 "7"="Did not say",
                 "9"="Middle Eastern",
                 "8"="South Asian")

tomato.pre + facet_grid(race)
tomato.post + facet_grid(race) 



# Histogram of cattle comfort


# Histogram of tomato comfort
cattle.pre <- ggplot(data=survey, aes(Q8, fill=as.factor(Q8)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Cattle case: Comfort PRE nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))
cattle.pre

# Print cattle histogram POST nudge
cattle.post <- ggplot(data=survey, aes(Q10, fill=as.factor(Q10)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Cattle case: Comfort POST nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))

cattle.prepost <- plot_grid(cattle.pre, cattle.post, nrow=1,
                            hjust = -1, labels = NULL,
                            label_size=12)
cattle.prepost





# Cattle faceting


cattle.pre + facet_grid(political)
cattle.post + facet_grid(political) 

cattle.pre + facet_grid(gender)
cattle.post + facet_grid(gender) 

cattle.pre + facet_grid(race)
cattle.post + facet_grid(race) 


# Histogram of wheat comfort

# Histogram of wheat comfort
wheat.pre <- ggplot(data=survey, aes(Q12, fill=as.factor(Q12)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Wheat case: Comfort PRE nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))
wheat.pre

# Print Wheat histogram POST nudge
wheat.post <- ggplot(data=survey, aes(Q14, fill=as.factor(Q14)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 600))+
  labs(title="Wheat case: Comfort POST nudge", 
       y="Respondents", x="")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Very comfortable",
                                 "Comfortable",
                                 "Neither comfortable nor uncomfortable",
                                 "Uncomfortable",
                                 "Very uncomfortable",
                                 "Don't know/not sure"))
wheat.post

wheat.prepost <- plot_grid(wheat.pre, wheat.post, nrow=1,
                           align='h',
                            hjust = -1, labels = NULL,
                            label_size=12)
wheat.prepost





# wheat faceting

wheat.pre + facet_grid(political)
wheat.post + facet_grid(political) 

wheat.pre + facet_grid(gender)
wheat.post + facet_grid(gender) 

wheat.pre + facet_grid(race)
wheat.post + facet_grid(race)


### DV2 data cleaning ########
Q16A_na1 <- na_if(survey$Q16A, 6)
Q16A_na2 <- na_if(Q16A_na1, 7)
Q16A_na3 <- na_if(Q16A_na2, 8)

Q16B_na1 <- na_if(survey$Q16B, 6)
Q16B_na2 <- na_if(Q16B_na1, 7)
Q16B_na3 <- na_if(Q16B_na2, 8)

Q18A_na1 <- na_if(survey$Q18A, 6)
Q18A_na2 <- na_if(Q18A_na1, 7)
Q18A_na3 <- na_if(Q18A_na2, 8)

Q18B_na1 <- na_if(survey$Q18B, 6)
Q18B_na2 <- na_if(Q18B_na1, 7)
Q18B_na3 <- na_if(Q18B_na2, 8)

# Combine split tradeoff questions into single columns
survey$Q16=rowSums(cbind(Q16A_na3,Q16B_na3),na.rm=TRUE)
survey$Q18=rowSums(cbind(Q18A_na3,Q18B_na3),na.rm=TRUE)


survey$Q16_plus=rowSums(cbind(survey$Q16A,survey$Q16B),na.rm=TRUE)
survey$Q18_plus=rowSums(cbind(survey$Q18A,survey$Q18B),na.rm=TRUE)

###### DV2 bar charts #########
# Histogram of pesticide tradeoff preferences, *BEFORE* answers 6, 7, 8 removed
pest.hist <- ggplot(data=survey, aes(Q16_plus, fill=as.factor(Q16_plus)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 500))+
  labs(title="Histogram for pesticide tradeoff, including NA's", 
       y="Respondents", x="") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Large increase in gene-edited crops",
                                 "A small increase in gene-edited crops",
                                 "No change",
                                 "A small decrease in gene-edited crops",
                                 "A large decrease in gene-edited crops",
                                 "Don't know/not sure",
                                 "Prefer not to answer b/c no information on ownership",
                                 "Prefer not to answer b/c other ways to avoid pesticides"))

pest.hist

# Histogram of biodiversity tradeoff preferences, *BEFORE* answers 6, 7, 8 removed
biodiv.hist <- ggplot(data=survey, aes(Q18_plus, fill=as.factor(Q18_plus)))+
  geom_bar()+
  coord_cartesian(ylim = c(0, 500))+
  labs(title="Histogram for biodiversity tradeoff, including NA's", 
       y="Respondents", x="") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_discrete(name = "Comfort level", 
                      labels = c("Large increase in gene-edited crops",
                                 "A small increase in gene-edited crops",
                                 "No change",
                                 "A small decrease in gene-edited crops",
                                 "A large decrease in gene-edited crops",
                                 "Don't know/not sure",
                                 "Prefer not to answer b/c no information on ownership",
                                 "Prefer not to answer b/c other ways to conserve biodiversity"))
      
biodiv.hist
# Playing around with facets
survey$Q29 <- factor(survey$Q29, levels = c("1", "2", "3", "4", "5"), 
                  labels = c("1" ="Very liberal", 
                             "2" = "Liberal",
                             "3" = "Moderate",
                             "4" = "Conservative",
                             "5" = "Very conservative"))

# Visualizing the two tradeoff questions together
tradeoff.both <- plot_grid(pest.hist, biodiv.hist, nrow=1,
                           align='h',
                           hjust = -1, labels = NULL,
                           label_size=12)
tradeoff.both



# POLITICAL faceting
political <-as.factor(survey$Q29)
levels(political) <-c("1" ="Very liberal", 
                     "2" = "Liberal",
                     "3" = "Moderate",
                     "4" = "Conservative",
                     "5" = "Very conservative")
pest.hist + facet_grid(political)
biodiv.hist + facet_grid(political)


# Gender faceting
gender <- as.factor(survey$Q26_clean)
levels(gender) <- c("1" = "Female",
                    "2" = "Male",
                    "3" = "Nonbinary",
                    "4" = "Transitioning")
pest.hist + facet_grid(gender)
biodiv.hist + facet_grid(gender)


# Race faceting
race <-as.factor(survey$Q30_clean)
levels(race) <-c("1"="East Asian",
                 "2"="Black",
                 "3"= "American Indian, First Nation, or Pacific Islander",
                 "4"="Hispanic",
                 "5"="Non-white mixed race",
                 "6"="White",
                 "7"="Did not say",
                 "9"="Middle Eastern",
                 "8"="South Asian")

pest.hist + facet_grid(race)
biodiv.hist + facet_grid(race)               


# Technology faceting

# Tech1, pesticides
tech1 <-as.factor(survey$Q22_1)
levels(tech1) <-c("1"="Strongly disagree",
                 "2"="Disagree",
                 "3"= "Neutral",
                 "4"="Agree",
                 "5"="Strongly agree",
                 "6"="Don't know/not sure")

pest.hist + facet_grid(tech1)

#Tech2 pesticides
tech2 <-as.factor(survey$Q22_2)
levels(tech2) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

pest.hist + facet_grid(tech2)

#Tech3 pesticides
tech3 <-as.factor(survey$Q22_3)
levels(tech3) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

pest.hist + facet_grid(tech3)

#Tech4 pesticides
tech4 <-as.factor(survey$Q22_4)
levels(tech4) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

pest.hist + facet_grid(tech4)

#Tech5, pesticides
tech5 <-as.factor(survey$Q22_5)
levels(tech5) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

pest.hist + facet_grid(tech5)

#Tech and biodiv

#Tech1, biodiv
tech1 <-as.factor(survey$Q22_1)
levels(tech1) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

biodiv.hist + facet_grid(tech1)

#Tech2, biodiv
tech2 <-as.factor(survey$Q22_2)
levels(tech2) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

biodiv.hist + facet_grid(tech2)


#Tech3, biodiv
tech3 <-as.factor(survey$Q22_3)
levels(tech3) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

biodiv.hist + facet_grid(tech3)

#Tech4, biodiv
tech4 <-as.factor(survey$Q22_4)
levels(tech4) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

biodiv.hist + facet_grid(tech4)


#Tech5, biodiv
tech5 <-as.factor(survey$Q22_5)
levels(tech5) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

biodiv.hist + facet_grid(tech5)




# Climate faceting

#climate 1, pesticide
climate1 <-as.factor(survey$Q21_1)
levels(climate1) <-c("1"="Strongly disagree",
                  "2"="Disagree",
                  "3"= "Neutral",
                  "4"="Agree",
                  "5"="Strongly agree",
                  "6"="Don't know/not sure")

pest.hist + facet_grid(climate1)

#climate 2, pesticide
climate2 <-as.factor(survey$Q21_2)
levels(climate2) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate2)

#climate 3, pesticide
climate3 <-as.factor(survey$Q21_3)
levels(climate3) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate3)


#climate 4, pesticide
climate4 <-as.factor(survey$Q21_4)
levels(climate4) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate4)


#climate 5, pesticide
climate5 <-as.factor(survey$Q21_5)
levels(climate5) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate5)


climate6 <-as.factor(survey$Q21_6)
levels(climate6) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate6)

climate7 <-as.factor(survey$Q21_7)
levels(climate7) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(climate7)


#Climate 1, biodiv
climate1 <-as.factor(survey$Q21_1)
levels(climate1) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate1)


#climate 2, biodiv
climate2 <-as.factor(survey$Q21_2)
levels(climate2) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate2)

#climate 3, biodiv
climate3 <-as.factor(survey$Q21_3)
levels(climate3) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate3)


#climate 4, biodiv
climate4 <-as.factor(survey$Q21_4)
levels(climate4) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate4)


#climate 5, biodiv
climate5 <-as.factor(survey$Q21_5)
levels(climate5) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate5)

#climate 6, biodiv
climate6 <-as.factor(survey$Q21_6)
levels(climate6) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate6)

#climate 7, biodiv
climate7 <-as.factor(survey$Q21_7)
levels(climate7) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

biodiv.hist + facet_grid(climate7)



# Trust faceting

#Trust, pesticide
trust1 <-as.factor(survey$Q20_1)
levels(trust1) <-c("1"="Strongly disagree",
                     "2"="Disagree",
                     "3"= "Neutral",
                     "4"="Agree",
                     "5"="Strongly agree",
                     "6"="Don't know/not sure")

pest.hist + facet_grid(trust1)

#Trust, biodiv
trust1 <-as.factor(survey$Q20_1)
levels(trust1) <-c("1"="Strongly disagree",
                   "2"="Disagree",
                   "3"= "Neutral",
                   "4"="Agree",
                   "5"="Strongly agree",
                   "6"="Don't know/not sure")

biodiv.hist + facet_grid(trust1)


##### Explore who opted out ####




# Filter by those who opted out of pesticide tradeoff, and gender











#### Factor analyses ##########

# CLIMATE SCALE
#remove NA's
Q21_1.na <-na_if(survey$Q21_1, 6)
Q21_2.na <-na_if(survey$Q21_2, 6)
Q21_3.na <-na_if(survey$Q21_3, 6)
Q21_4.na <-na_if(survey$Q21_4, 6)
Q21_5.na <-na_if(survey$Q21_5, 6)
Q21_6.na <-na_if(survey$Q21_6, 6)
Q21_7.na <-na_if(survey$Q21_7, 6)

# Load data
climate <- data.frame(Q21_1.na, Q21_2.na, Q21_3.na, Q21_4.na, Q21_5.na, Q21_6.na, Q21_7.na)
keys.climate <- c(1,-1,-1, 1, 1, 1, -1)
reverse.code(keys.climate, climate)

# histogram to check visually if normally distributed
hist(climate)

# shapiro test to check if normally distributed
# if p-values are less than 0.5, then NOT normally distributed
shapiro.test(survey$Q21_1)
shapiro.test(survey$Q21_2)
shapiro.test(survey$Q21_3)
shapiro.test(survey$Q21_4)
shapiro.test(survey$Q21_5)
shapiro.test(survey$Q21_6)
shapiro.test(survey$Q21_7)

# KMO test--to check sampling adequacy, to see if worth looking at a correlation matrix
# 0.5 is unacceptable, 0.6-.69 is mediocre, 0.70-.79 middling, 0.8-.89 meritorious, etc.
KMO(climate)

# Bartlett's sphericity test
# tests if there are worthwhile correlations between items... 
# if P value is below significance level, data is appropriate for factor analysis
cortest.bartlett(climate)

# Check correlation matrix
cor(climate)
# Bartlett test to see if factor analysis is appropriate
cortest.bartlett(climate)


# PCA
climate.pca <- prcomp(climate, center=TRUE, scale.=TRUE)
summary(climate.pca)

# Scree test to see where eigenvalues fall for both PC and FA 
scree(climate)

# Parallel analysis
parallel = fa.parallel(climate, fm = "pa", fa = "fa")
print(parallel)

# Factor analysis
climate.fa <- fa(climate, nfactors=2, scores=TRUE, rotate="oblimin", fm="pa")
print.psych(climate.fa, cut = 0.3, sort = TRUE)
fa.diagram(climate.fa)

# Cronbach alpha scores
climate.a <- data.frame(Q21_6.na, Q21_4.na, Q21_1.na, Q21_5.na)
climate.b <- data.frame(Q21_2.na, Q21_3.na, Q21_7.na)
alpha(climate.a)
alpha(climate.b)


# GLOBALIZATION SCALE

#remove NA's
Q24_1.na <-na_if(survey$Q24_1, 6)
Q24_2.na <-na_if(survey$Q24_2, 6)
Q24_3.na <-na_if(survey$Q24_3, 6)
Q24_4.na <-na_if(survey$Q24_4, 6)


# Load data
global <- data.frame(Q24_1.na, Q24_2.na, Q24_3.na, Q24_4.na)
keys.global <- c(1,-1, 1, 1)
reverse.code(keys.global, global)


# shapiro test to check if normally distributed
# if p-values are less than 0.5, then NOT normally distributed
shapiro.test(survey$Q24_1)
shapiro.test(survey$Q24_2)
shapiro.test(survey$Q24_3)
shapiro.test(survey$Q24_4)

# KMO test--to check sampling adequacy, to see if worth looking at a correlation matrix
# 0.5 is unacceptable, 0.6-.69 is mediocre, 0.70-.79 middling, 0.8-.89 meritorious, etc.
KMO(global)

# Bartlett's sphericity test
# tests if there are worthwhile correlations between items... 
# if P value is below significance level, data is appropriate for factor analysis
cortest.bartlett(global)

# Check correlation matrix
cor(global)
# Bartlett test to see if factor analysis is appropriate
cortest.bartlett(global)

# PCA
global.pca <- prcomp(global, center=TRUE, scale.=TRUE)
summary(global.pca)

# Scree test to see where eigenvalues fall for both PC and FA 
scree(global) 

# Parallel analysis
parallel = fa.parallel(global, fm = "pa", fa = "fa")
print(parallel)

# Factor analysis using Psych package
global.fa <- fa(global, nfactors=2, scores=TRUE,rotate="oblimin", fm="pa")
print.psych(global.fa, cut = 0.3, sort = TRUE)
fa.diagram(global.fa)

global.a <- data.frame(Q24_1.na, Q24_4.na, Q24_3.na)
alpha(global.a)

# TRUST SCALE

#remove NA's
Q20_1.na <-na_if(survey$Q20_1, 6)
Q20_2.na <-na_if(survey$Q20_2, 6)
Q20_3.na <-na_if(survey$Q20_3, 6)
Q20_4.na <-na_if(survey$Q20_4, 6)

trust <- data.frame(Q20_1.na, Q20_2.na, Q20_3.na, Q20_4.na)



scree(trust) 
trust.fa <- fa(trust, nfactors=1, scores=TRUE,rotate="oblimin", fm="pa")
print.psych(trust.fa, cut = 0.3, sort = TRUE)

trust.a <- data.frame(Q20_1.na, Q20_2.na, Q20_3.na, Q20_4.na)
alpha(trust.a)

# TECH SCALE

#remove NA's
Q22_1.na <-na_if(survey$Q22_1, 6)
Q22_2.na <-na_if(survey$Q22_2, 6)
Q22_3.na <-na_if(survey$Q22_3, 6)
Q22_4.na <-na_if(survey$Q22_4, 6)
Q22_5.na <-na_if(survey$Q22_5, 6)

tech <- data.frame(Q22_1.na, Q22_2.na, Q22_3.na, Q22_4.na, Q22_5.na)
keys.tech <- c(1, 1, -1, -1, 1)
reverse.code(keys.tech, tech)

scree(tech)
tech.fa <- fa(tech, nfactors=2, scores=TRUE,rotate="oblimin", fm="pa")
print.psych(tech.fa, cut = 0.3, sort = TRUE)

tech.a <- data.frame(Q22_5.na, Q22_2.na, Q22_1.na)
tech.b <- data.frame(Q22_4.na, Q22_3.na)

alpha(tech.a)
alpha(tech.b)

# TOXICS SCALE

#remove NA's
Q23_1.na <-na_if(survey$Q23_1, 6)
Q23_2.na <-na_if(survey$Q23_2, 6)
Q23_3.na <-na_if(survey$Q23_3, 6)
Q23_4.na <-na_if(survey$Q23_4, 6)
Q23_5.na <-na_if(survey$Q23_5, 6)
Q23_6.na <-na_if(survey$Q23_6, 6)

toxics <- data.frame(Q23_1.na, Q23_2.na, Q23_3.na, Q23_4.na, Q23_5.na, Q23_6.na)
keys.toxics <- c(-1, 1, 1, 1, 1, -1)
reverse.code(keys.toxics, toxics)

scree(toxics)
toxics.fa <- fa(toxics, nfactors=2, scores=TRUE,rotate="oblimin", fm="pa")
print.psych(toxics.fa, cut = 0.3, sort = TRUE)

toxics.a <- data.frame(Q23_3.na, Q23_2.na, Q23_5.na)
toxics.b <- data.frame(Q23_6.na, Q23_4.na)
alpha(toxics.a)
alpha(toxics.b)

# GREEN REVOLUTION SCALE

#remove NA's
Q25_1.na <-na_if(survey$Q25_1, 6)
Q25_2.na <-na_if(survey$Q25_2, 6)
Q25_3.na <-na_if(survey$Q25_3, 6)
Q25_4.na <-na_if(survey$Q25_4, 6)
Q25_5.na <-na_if(survey$Q25_5, 6)
Q25_6.na <-na_if(survey$Q25_6, 6)
Q25_7.na <-na_if(survey$Q25_7, 6)
Q25_8.na <-na_if(survey$Q25_8, 6)

green <- data.frame(Q25_2.na, Q25_3.na, Q25_4.na, Q25_5.na, Q25_6.na, Q25_7.na, Q25_8.na)
keys.green <- c(1, 1, -1, 1, -1, -1, -1)
reverse.code(keys.green, green)

scree(green)
green.fa <- fa(green, nfactors=2, scores=TRUE,rotate="oblimin", fm="pa")
print.psych(green.fa, cut = 0.3, sort = TRUE)

green.a <- data.frame(Q25_4.na, Q25_6.na, Q25_8.na, Q25_7.na)
green.b <- data.frame(Q25_5.na, Q25_2.na)

alpha(green.a)
alpha(green.b)


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

##### Prepare data for regression analysis #########################


# add factor scores to dataset for regression analysis
regdata <- cbind(survey, climate.fa$scores, tech.fa$scores, toxics.fa$scores, green.fa$scores, global.fa$scores, trust.fa$scores) # ADD remaining/final!

#Labeling the data
# names(regdata) <- c() #FINISH CODE

# Turn dependent variable 1 into dummy variables for binomial logistic reg
survey$Q16_1 <- as.factor(survey$Q16 == "1")
survey$Q16_2 <- as.factor(survey$Q16 == "2")
survey$Q16_3 <- as.factor(survey$Q16 == "3")
survey$Q16_4 <- as.factor(survey$Q16 == "4")
survey$Q16_5 <- as.factor(survey$Q16 == "5")
survey$Q16_6 <- as.factor(survey$Q16 == "6")
survey$Q16_7 <- as.factor(survey$Q16 == "7")
survey$Q16_8 <- as.factor(survey$Q16 == "8")


####### Regression models ################
  

# Binomial logistic regression all responses to Q16
pest.model0 = glm(survey$Q16 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="gaussian")
summary(pest.model0)
rsq(pest.model, adj = T)




# Binomial logistic regression for those who answered 1 for Q16
pest.model1 = glm(survey$Q16_1 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model1)
rsq(pest.model, adj = T)

# Binomial logistic regression for those who answered 2 for Q16
pest.model2 = glm(survey$Q16_2 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model1)

# Binomial logistic regression for those who answered 3 for Q16
pest.model3 = glm(survey$Q16_3 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model3)

# Binomial logistic regression for those who answered 4 for Q16
pest.model4 = glm(survey$Q16_4 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model4)

# Binomial logistic regression for those who answered 2 for Q16
pest.model5 = glm(survey$Q16_5 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model5)

# Binomial logistic regression for those who answered 2 for Q16
pest.model6 = glm(survey$Q16_6 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model6)

# Binomial logistic regression for those who answered 2 for Q16
pest.model7 = glm(survey$Q16_7 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model7)

# Binomial logistic regression for those who answered 2 for Q16
pest.model8 = glm(survey$Q16_8 ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="binomial")
summary(pest.model8)


#calculate VIF scores to see where there is multicollinearity
vif(pest.model0)



# Regressions to explore what explains  
tomato.model1 = glm(Q4.na ~ survey$Q26_clean + Q27 + climate.fa$scores + tech.fa$scores 
                  +toxics.fa$scores + green.fa$scores + global.fa$scores+trust.fa$scores, 
                  data=regdata, family="gaussian")
summary(tomato.model1)



# PENALIZED REGRESSION

# Getting the independent variable
x_vars2 <- data.matrix(survey[, c("Q26", "Q27", "Q28", "Q29", "Q30", "Q31", "Q32", "Q33")])
y_var2 <- survey$Q16
lambda_seq <- 10^seq(2, -2, by = -.1)

# Splitting the data into test and train
set.seed(86)
train = sample(1:nrow(x_vars2), nrow(x_vars2)/2)
x_test = (-train)
y_test = y_var2[x_test]

cv_output <- cv.glmnet(x_vars2[train,], y_var2[train],
                       alpha = 1, lambda = lambda_seq, 
                       nfolds = 5)

# identifying best lamda
best_lam <- cv_output$lambda.min
best_lam


# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(x_vars2[train,], y_var2[train], alpha = 1, lambda = best_lam)

##PROBLEM HERE--TOO MANY NA's!!!!!

pred <- predict(lasso_best, s = best_lam, newx = x_vars2[x_test,])

final <- cbind(y_var2[y_test], pred)
# Checking the first six obs
head(final)

actual <- test$actual
preds <- test$predicted
rss <- sum((preds - actual) ^ 2)
tss <- sum((actual - mean(actual)) ^ 2)
rsq <- 1 - rss/tss
rsq

coef(lasso_best)











