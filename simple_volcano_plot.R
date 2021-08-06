library(compare)

#set your working directory
setwd("/path/to/your/working/directory/")

# read the csv file into a dataframe
data <- read.csv(file = 'your_file.csv')
head(data)

#analysing attributes of the data
colnames(data)
dim(data)

# keep only the fields needed for the plot (change the name accordingly)
# adjusted p value = significance
data <- data[,c('id','pval','foldChange')]
head(data)

#plotting a simple volcano-plot
#change the limits of the x & y-axis accordingly
with(data, plot(foldChange, -log10(pval),pch=20, main="Insert Title of the plot here", xlim=c(-8,8),ylim=c(0,8)))

#highlighting points which satisfy the condition
#this will highlight UP & DOWN-regulated genes in your plot
#change the thresholds for fold change & p-value accordingly
with(subset(data, pval<.05 & abs(foldChange)>2), points(log2(foldChange), -log10(pval), pch=20, col="red"))
with(subset(data, pval<.05 & abs(foldChange)<0.5), points(log2(foldChange), -log10(pval), pch=20, col="blue"))
dim(data)

#save differentially expressed gene list to a file
a1 <- subset(data, pval<.05 & abs(foldChange) > 2)
head(a1)
dim(a1)
write.csv(a1,file = "UPREGULATED.csv")
a2 <- subset(data, pval<.05 & abs(foldChange) < 0.5)
head(a2)
dim(a2)
write.csv(a2,file = "DOWNREGULATED.csv")
