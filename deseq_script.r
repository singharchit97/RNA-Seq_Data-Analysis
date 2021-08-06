#loading library DESeq
library( "DESeq" )

#setting up the working directory
setwd("/path/to/your/working/directory")

#reading the .csv file into a dataframe
#should contain readcounts (raw) from all the samples (experimental & control groups)
seq_data <- read.csv(file = 'deseq.csv',header = TRUE)
#to view the first six rows of the dataframe
head(seq_data)

#writing a new dataframe from the previous dataframe
#new dataframe contains the info on the condition(could be Wild type vs mutant/ could be untreated vs treated sample)
#club the duplicates together under the same condition name
#also mention if data is single end/ paired end sequenced
SeqdataDesign = data.frame(row.names = colnames(seq_data),condition = c( "control", "experiment" ),libType = c( "paired-end", "paired-end" ) )
#view the dataframe
SeqdataDesign

#for sequencing data with single end and paired end reads
#Process the paired end data seperately from the single end data
#For now if your data contains both type of reads, then use the following snippet of code to seperate them out into two datasets
pairedSamples = SeqdataDesign$libType == "paired-end" #mention "single-end for single end sequencing"
#creating new dataframe which contains only paired end sample data
countTable = seq_data[ , pairedSamples ]
condition = SeqdataDesign$condition[ pairedSamples ]
head(countTable)
condition

#creating new dataset containing the dataframe (countTable) and a factor(condition)
cds = newCountDataSet( countTable, condition )
#Normalization step
#the function estimateSizeFactors estimates the size factors from the count data
#we divide each column of the count table by the size factor for this column, the count values are brought
#to a common scale, making them comparable
cds = estimateSizeFactors(cds)
head(counts(cds, normalized=TRUE))
head(cds)
write.csv(counts(cds, normalized=TRUE),'normalized_deseq_counts.csv')

#Variance estimation
#DESeq relies on an estimation of the typical relationship between the data’s variance and their mean,
#or, equivalently, between the data’s dispersion and their mean
#estimating dispersion
cds = estimateDispersions(cds)
head(cds)

#to see the differential expression between conditions“WT_Rv”and“M_Rv”,we simply call the function nbinomTest
#returns a data frame with the pvalues and other useful information
res = nbinomTest(cds,"control", "experiment")
head(res)
#plotting the log 2 fold changes against the mean normalised counts, colouring in red those genes that are significant at 10% FDR
plotMA(res)
#plotting a frequency histogram of all the p-values
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
#filtering the most significantly differentially expressed genes
resSig = res[ res$padj < 0.05, ] #cutoff used here is 0.05 
head( resSig[ order(resSig$pval), ] )
#filter out the most underexpressed genes
head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
#filter out the most overexpressed genes
head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] )
#writimg the dataframe into a .csv file for reference
write.csv(res,file = "control-vs-experiment_result_file.csv")
