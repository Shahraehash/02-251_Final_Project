filedata <- read.csv("DataFilePearsonComparison.csv")
filedata <- filedata[1:60483,2:3]
library("ggpubr")
ggscatter(filedata, x = "InvasiveDuctalBreastAvg_FPKM", y = "InvasiveDuctalBreastAvg_Counts",
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall",
          xlab = "FPKM Data", ylab = "Count Data")
ggqqplot(filedata$InvasiveDuctalBreastAvg_Counts, ylab = "Count Data")
ggqqplot(filedata$InvasiveDuctalBreastAvg_FPKM, ylab = "FPKM Data")

seqdata <- read.csv("SequencingTypePearsonComparison.csv")
seqdata <- seqdata[1:60483,2:3]
library("ggpubr")
ggscatter(seqdata, x = "AngiosarcomaAvg_HTSeq", y = "AngiosarcomaAvg_STAR",
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
          xlab = "HTSeq", ylab = "STAR")