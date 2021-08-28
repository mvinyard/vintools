#! R script
# https://lashlock.github.io/compbio/R_presentation.html
# load DESeq2
library(DESeq2)

# arguments
args = commandArgs(trailingOnly=TRUE)
outsdir = args[1] # "./"
run_name = args[2] # "something"

print(c(outsdir, run_name))

# load data
countData <- read.csv('./tmp/gex_counts.csv', header = TRUE, sep = ",")
metaData <- read.csv('./tmp/gex_meta.csv', header=TRUE, sep=",")
DESeq2_Data <- DESeqDataSetFromMatrix(countData=countData,
                                      colData=metaData,
                                      design=~condition, tidy = TRUE)

# run DESeq2
DESeq2_Data <- DESeq(DESeq2_Data)
DESeq2_Data_results <- results(DESeq2_Data) # by default from tutorial; these vars are res and dds
DESeq2_vsData <- vst(DESeq2_Data, blind=FALSE)
pcs = plotPCA(DESeq2_vsData, intgroup="condition",  returnData=TRUE)

# write output
write.table(data.frame(DESeq2_Data_results), file=paste(outsdir, "/", run_name, "_", "DESeq2_results_df.csv", sep=""))
write.table(assay(DESeq2_vsData), paste(outsdir, "/", run_name, "_", "DESeq2_vsData.csv", sep=""))
write.table(pcs, paste(outsdir, "/", run_name, "_", "DESeq2_pca.csv", sep=""))