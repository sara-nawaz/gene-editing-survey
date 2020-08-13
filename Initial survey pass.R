library(ggplot2)
library(psych)
survey <- read.csv('Perceptions of gene editing in agriculture_August 3 2020 EXTRA CLEAN numeric.csv', sep= ',', stringsAsFactors=FALSE, header=TRUE)

# Dependent variable 1: Comfort level with three applications

# Tomatoes

# Average for tomatoes
mean(survey$Q4)
# Average for tomatoes with nudge (sugar beet gene)
mean(survey$Q6)

# Histogram for tomatoes
ggplot(data=survey, aes(Q4))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for Tomatoes", x="Comfort level", y="Respondents")

# Histogram for tomatoes + nudge
ggplot(data=survey, aes(Q6))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for Tomatoes+nudge", x="Comfort level", y="Respondents")

# Cattle

# Average for cattle
mean(survey$Q8)
# Average for cattle with nudge (more productive)
mean(survey$Q10)

# Histogram for cattle
ggplot(data=survey, aes(Q8))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for cattle", x="Comfort level", y="Respondents")

# Histogram for cattle+nudge
ggplot(data=survey, aes(Q10))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for cattle+nudge", x="Comfort level", y="Respondents")

# Wheat

# Average for wheat
mean(survey$Q12)
# Average for wheat with nudge (blueberries not wheat)
mean(survey$Q14)

# Histogram for wheat
ggplot(data=survey, aes(Q12))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for wheat", x="Comfort level", y="Respondents")

# Histogram for wheat plus nudge
ggplot(data=survey, aes(Q14))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for wheat+nudge", x="Comfort level", y="Respondents")



# Dependent variable 2: Tradeoff preferences

# Combine split tradeoff questions into single columns
survey$Q16=rowSums(cbind(survey$Q16A,survey$Q16B),na.rm=TRUE)
survey$Q18=rowSums(cbind(survey$Q18A,survey$Q18B),na.rm=TRUE)

# Histogram of pesticide tradeoff preferences
ggplot(data=survey, aes(Q16))+
geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for pesticide tradeoff", x="Comfort level", y="Respondents")

# Histogram of biodiversity tradeoff preferences
ggplot(data=survey, aes(Q18))+
  geom_histogram()+
  stat_bin(binwidth=1)+
  labs(title="Histogram for biodiv tradeoff", x="Comfort level", y="Respondents")


# Factor analyses

# Climate scale
# Load data
climate <- data.frame(survey$Q21_1, survey$Q21_2, survey$Q21_3, survey$Q21_4, survey$Q21_5, survey$Q21_6, survey$Q21_7)
# PCA
climate.pca <- prcomp(climate, center=TRUE, scale.=TRUE)
summary(climate.pca)
# Note where slope of plotted EVs flattens out
screeplot(climate.pca, type = "l", npcs = 15, main = "Screeplot of all PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# Factor analysis
climate.fa <- factanal(climate, factors= 2, rotation="varimax")
climate.fa
climate.fa$scores

# Globalization scale
# Load data
global <- data.frame(survey$Q24_1, survey$Q24_2, survey$Q24_3, survey$Q24_4)
# PCA
global.pca <- prcomp(global, center=TRUE, scale.=TRUE)
summary(global.pca)
#Note where slope of plotted EVs flattens out
screeplot(global.pca, type = "l", npcs = 15, main = "Screeplot of all PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# Factor analysis
global.fa <- factanal(global, factors= 1, rotation="varimax")
global.fa

# Green revolution scale
# Load data
green <- data.frame(survey$Q25_2, survey$Q25_3, survey$Q25_4, survey$Q25_5, survey$Q25_6, survey$Q25_7)
# PCA
green.pca <- prcomp(green, center=TRUE, scale.=TRUE)
summary(green.pca)
# Note where slope of plotted EVs flattens out
screeplot(green.pca, type = "l", npcs = 15, main = "Screeplot of all PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# Factor analysis
green.fa <- factanal(green, factors= 1)
green.fa

# Tech, globalization, green rev scale ("views on agricultural technology" scale?)
# Load data
tgg <- data.frame(survey$Q22_1, survey$Q22_2, survey$Q22_3, survey$Q22_4, survey$Q22_5, 
                  survey$Q25_2, survey$Q25_3, survey$Q25_4, survey$Q25_5, survey$Q25_6, survey$Q25_7,
                  survey$Q24_1, survey$Q24_2, survey$Q24_3, survey$Q24_4)
# PCA
green.pca <- prcomp(green, center=TRUE, scale.=TRUE)
summary(green.pca)
# Note where slope of plotted EVs flattens out
screeplot(green.pca, type = "l", npcs = 15, main = "Screeplot of all PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# Factor analysis
green.fa <- factanal(green, factors= 1)
green.fa


