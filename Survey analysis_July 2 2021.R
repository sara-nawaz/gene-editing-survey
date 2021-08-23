library(ggplot2)
library(ggcorrplot)
library(tidyr)
library(tidyverse)
library(psych) # psychometrics
library(MVN) # multivariate normality
library(sjPlot) #visualization
library(cowplot)
library(scales)  # for percentage scales
library(glmnet)
library(leaps) # stepwise regression
library(car) # for VIF function, also recoding
library(RColorBrewer)
library(ggpubr)
library(questionr)
library(broom)
library(expss)


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
income <-survey$Q33

#### Mean scores of factors for regressions (NOT full factor analysis) ##########

# CLIMATE SCALE
#remove NA's
survey$Q21_1.na <-na_if(survey$Q21_1, 6)
survey$Q21_2.na <-na_if(survey$Q21_2, 6)
survey$Q21_3.na <-na_if(survey$Q21_3, 6)
survey$Q21_4.na <-na_if(survey$Q21_4, 6)
survey$Q21_5.na <-na_if(survey$Q21_5, 6)
survey$Q21_6.na <-na_if(survey$Q21_6, 6)
survey$Q21_7.na <-na_if(survey$Q21_7, 6)

climate1mean=rowMeans(cbind(survey$Q21_6.na, survey$Q21_4.na, survey$Q21_1.na, survey$Q21_5.na),na.rm=TRUE)
climate2mean=rowMeans(cbind(survey$Q21_2.na, survey$Q21_3.na, survey$Q21_7.na),na.rm=TRUE)

# GLOBALIZATION SCALE

#remove NA's
survey$Q24_1.na <-na_if(survey$Q24_1, 6)
survey$Q24_3.na <-na_if(survey$Q24_3, 6)
survey$Q24_4.na <-na_if(survey$Q24_4, 6)

globalmean=rowMeans(cbind(survey$Q24_1.na, survey$Q24_4.na, survey$Q24_3.na),na.rm=TRUE)


# TRUST SCALE

#remove NA's
survey$Q20_1.na <-na_if(survey$Q20_1, 6)
survey$Q20_2.na <-na_if(survey$Q20_2, 6)
survey$Q20_3.na <-na_if(survey$Q20_3, 6)

trustmean=rowMeans(cbind(survey$Q20_1.na, survey$Q20_2.na, survey$Q20_3.na),na.rm=TRUE)


# TECH SCALE

#remove NA's
survey$Q22_2.na <-na_if(survey$Q22_2, 6)
survey$Q22_5.na <-na_if(survey$Q22_5, 6)

techmean=rowMeans(cbind(survey$Q22_5.na, survey$Q22_2.na),na.rm=TRUE)


# TOXICS SCALE

#remove NA's
survey$Q23_2.na <-na_if(survey$Q23_2, 6)
survey$Q23_3.na <-na_if(survey$Q23_3, 6)
survey$Q23_5.na <-na_if(survey$Q23_5, 6)

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

green1mean=rowMeans(cbind(survey$Q25_4.na, survey$Q25_6.na, survey$Q25_8.na, survey$Q25_7.na),na.rm=TRUE)
green2mean=rowMeans(cbind(survey$Q25_5.na, survey$Q25_2.na, survey$Q25_3.na),na.rm=TRUE)



###### DV1 data cleaning #####
# Combine and label variables, convert "don't know/not sure" responses to NA's

# Tomatoes
# pre-nudge
survey$Q4.na <- na_if(survey$Q4, 6)

# Cattle
# pre-nudge
survey$Q8.na <- na_if(survey$Q8, 6)

# Wheat
# pre-nudge
survey$Q12.na <- na_if(survey$Q12, 6)


#RECODE???
#library(car)
#Q4.na = recode(Q4.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q6.na = recode(Q6.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q8.na = recode(Q8.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q10.na = recode(Q10.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q12.na = recode(Q12.na, "1=5; 2=4; 3=3; 4=2; 5=1") 
#Q14.na = recode(Q14.na, "1=5; 2=4; 3=3; 4=2; 5=1") 

comfort=rowMeans(cbind(survey$Q4.na, survey$Q8.na, survey$Q12.na),na.rm=TRUE)


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

######Correlation matrices####
survey2 <-cbind(age_full, gender_full, political_full, religious_full, 
                uni_full, income,
                trustmean, familiar_ge, globalmean, 
                climate2mean, green1mean, green2mean)

colnames(survey2) <- c("Age","Gender (female)", "Political (conservatism)", 
                       "Religiosity", "Formal education",
                       "Income", "Trust", "Familiarity with gene editing",
                       "Corporate criticism", "Climate ambivalence",
                       "Green Revolution criticism", "Green Revolution optimism")
corr <- cor(survey2, use="complete.obs")
corrplot <- ggcorrplot(corr, hc.order = TRUE, type ="lower")+
  theme(text=element_text(size=16,  family="Times"))


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
pest_step <- glm(formula = pestdv ~ green2mean + green1mean +
                   globalmean + climate2mean + trustmean + 
                   familiar_ge +
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
biodiv_step <- glm(formula = biodivdv ~ green2mean + green1mean +
                     globalmean + climate2mean + trustmean + 
                     familiar_ge +
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

pest_ordered <- plot_model(pest_step, title = "Pesticide Use vs. GE", show.values = TRUE, dot.size = 3, value.offset=0.3,
           vline.color = "grey", colors= c("blue", "red"), geom.label.color = "black") +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 16, base_family = "Times")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 0.3, 1, 0.7), "cm"))

biodiv_ordered <- plot_model(biodiv_step, title = "Biodiversity loss vs. GE", show.values = TRUE, dot.size = 3, value.offset=0.3,
           vline.color = "grey", colors= c("blue", "red")) +ylim(0.5, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 16, base_family = "Times")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 1, 1, 0.3), "cm"))

ordered_combined <- plot_grid(pest_ordered, 
          biodiv_ordered)


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
pest_step_bi <- glm(pestdv_bi ~ green2mean + green1mean +
                      globalmean + climate2mean + trustmean + 
                      familiar_ge +
                      gender_full + uni_full + age_full +  
                      religious_full + political_full + income,
                         family = binomial, data = pestdata_bi)
summary(pest_step_bi)

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(pest_step_bi)

pest_bi_or <- odds.ratio(pest_step_bi, level=0.95)
pest_bi_or 
write.csv(pest_bi_or, file="pest binomial or.csv")



##Biodiv tradeoff
biodiv_step_bi <- glm(biodivdv_bi ~ green2mean + green1mean +
                        globalmean + climate2mean + trustmean + 
                        familiar_ge +
                        gender_full + uni_full + age_full +  
                        religious_full + political_full + income,
                    family = binomial, data = biodivdata_bi)
summary(biodiv_step_bi)


#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(biodiv_step_bi)

biodiv_bi_or <- odds.ratio(biodiv_step_bi, level=0.95)
biodiv_bi_or
write.csv(biodiv_bi_or, file="biodiv binomial or.csv")





pest_bi <- plot_model(pest_step_bi, title = "Opting out of Pesticide Tradeoff", show.values = TRUE, dot.size = 3, value.offset=0.3,
                           vline.color = "grey", colors= c("red", "blue")) +ylim(0.25, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 16, base_family = "Times")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 0.3, 1, 0.7), "cm"))

biodiv_bi <- plot_model(biodiv_step_bi, title = "Opting out of Biodiversity Tradeoff", show.values = TRUE, dot.size = 3, value.offset=0.3,
                             vline.color = "grey", colors= c("red", "blue")) +ylim(0.25, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_sjplot(base_size = 16, base_family = "Times")+
  theme(axis.title.x=element_blank())+
  theme(plot.margin = unit(c(1, 1, 1, 0.3), "cm"))

bi_combined <- plot_grid(pest_bi, 
                              biodiv_bi, 
                              label_size = 12)


save_plot("./bi combined.png", bi_combined, base_height=8, ncol=1, nrow=1)





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
tomato_step <- glm(formula = tomatodv ~ green2mean + green1mean +
                     globalmean + climate2mean + trustmean + 
                     familiar_ge +
                     gender_full + uni_full + age_full +  
                     religious_full + political_full + income,
                   data = tomatodata, binomial(link="logit"))
summary(tomato_step)

tomato_or <- odds.ratio(tomato_step, level=0.95)
tomato_or
write.csv(tomato_or, file="tomato_ordinal_or.csv")

#to check multicollinearity--calculate VIF scores for regressors from best stepwise outcome model
vif(tomato_step)


##cattle dv
cattle_step <- glm(formula = cattledv ~ green2mean + green1mean +
                     globalmean + climate2mean + trustmean + 
                     familiar_ge +
                     gender_full + uni_full + age_full +  
                     religious_full + political_full + income,
                   data = cattledata, binomial(link="logit"))
summary(cattle_step)
cattle_or <- odds.ratio(cattle_step, level=0.95)
cattle_or
write.csv(cattle_or, file="cattle_ordinal_or.csv")

vif(cattle_step)

##wheat dv
wheat_step <- glm(formula = wheatdv ~ green2mean + green1mean +
                    globalmean + climate2mean + trustmean + 
                    familiar_ge +
                    gender_full + uni_full + age_full +  
                    religious_full + political_full + income,
                  data = wheatdata, binomial(link="logit"))
summary(wheat_step)

wheat_or <- odds.ratio(wheat_step, level=0.95)
wheat_or
write.csv(wheat_or, file="wheat_ordinal_or.csv")

vif(wheat_step)


tomato_model <- plot_model(tomato_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey",colors= c("blue", "red")) +ylim(0.4, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())
  

cattle_model <- plot_model(cattle_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey",colors= c("blue", "red")) +ylim(0.4, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())

wheat_model <- plot_model(wheat_step, title = "", show.values = TRUE, dot.size = 1.5, value.offset=0.3,
           vline.color = "grey", colors= c("blue", "red")) +ylim(0.3, 2) +
  scale_x_discrete(labels=list(
    green2mean = "Green Revolution optimism",
    green1mean = "Green Revolution criticism",
    globalmean = "Corporate criticism",
    climate2mean = "Climate ambivalence",
    trustmean = "Trust",
    familiar_ge = "Familiarity with gene editing",
    gender_full2 = "Gender (male)",
    uni_full = "Education",
    age_full = "Age",
    religious_full = "Religiosity",
    political_full = "Political conservatism",
    income = "Income"
  ))+
  theme_gray(base_size = 18, base_family = "Times") +
  theme(axis.title.x=element_blank())
save_plot("./wheatplot.png", wheat_model, base_height=8)

  



cases_combined <- plot_grid(tomato_model, cattle_model, wheat_model, 
                         labels = c('Tomato', 'Cattle', 'Wheat'), 
                         label_size = 12, ncol=1)
save_plot("./cases combined.png", cases_combined, base_height=8)




