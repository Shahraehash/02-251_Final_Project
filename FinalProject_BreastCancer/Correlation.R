#Cancer comparison file
data_Cancer <- read.csv("CancerComparison.csv")
data_Cancer <- data_Cancer[1:60483,2:9]

#File Comparison File
data_File <- read.csv("DataFileComparison.csv")
data_File <- data_File[1:60483,2:9]

#Sequencing Comparison File
data_Sequencing <- read.csv("SequencingComparison.csv")
data_Sequencing <- data_Sequencing[1:60483, 2:9]

#Correlation between Cancers
cor_Cancer <- cor(data_Cancer, use = "complete.obs")
round(cor_Cancer, 2)
mean(cor_Cancer)
sd(cor_Cancer)
library(corrplot)
corrplot(cor_Cancer, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5)
library(ggpubr)
t.test(data_Cancer$InvasiveDuctalBreastCancer_Avg, data_Cancer$Angiosarcoma_Avg, alternative = "two.sided", var.equal = FALSE)

#Correlation between Files
cor_File <- cor(data_File, use = "complete.obs")
round(cor_File, 2)
mean(cor_File)
sd(cor_File)
library(corrplot)
corrplot(cor_File, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5)


#Correlation between Sequencing
cor_Sequencing <- cor(data_Sequencing, use = "complete.obs")
round(cor_Sequencing, 2)
mean(cor_Sequencing)
sd(cor_Sequencing)
library(corrplot)
corrplot(cor_Sequencing, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5)

