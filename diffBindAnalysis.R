library(DiffBind)
library(dplyr)

setwd("/path/to/pipeline")
diffBindDir <- file.path("diffBind")
dir.create(diffBindDir, recursive = TRUE, showWarnings = FALSE)

# Create a sample sheet
setwd
sampleSheet <- data.frame(
  SampleID = c("ALL1", "ALL2", "WT1", "WT2"),
  Tissue = rep("proB", 4),
  Factor = rep("ATAC", 4),
  Condition = c("Leukemia", "Leukemia", "WildType", "WildType"),
  Replicate = c(1, 2, 1, 2),
  bamReads = c("Data/B-ALL1_filtered.bam", "Data/B-ALL2_filtered.bam", "Data/WT1_filtered.bam", "data/WT2_filtered.bam"),
  ControlID = c(NA, NA, NA, NA),
  bamControl = c(NA, NA, NA, NA),
  Peaks = c("macs2_peaks/filtered_p10/B-ALL1_peaks.bed", "macs2_peaks/filtered_p10/B-ALL2_peaks.bed", "macs2_peaks/filtered_p10/WT1_peaks.bed", "macs2_peaks/filtered_p10/WT2_peaks.bed"),
  PeakCaller = rep("macs", 4)
)
write.csv(sampleSheet, file = "diffBind/sampleSheet.csv", row.names = FALSE)

setwd(diffBindDir)
# DiffBind Analysis
ALLdba <-dba(minOverlap=1,sampleSheet="sampleSheet.csv",config=data.frame(RunParallel=TRUE, reportInit="ATAC",DataType=DBA_DATA_GRANGES,AnalysisMethod=DBA_DESEQ2, minQCth=20, bCorPlot=FALSE,th=0.05, bUsePval=FALSE),peakCaller="macs", peakFormat="macs", skipLines=1,bAddCallerConsensus=FALSE,bRemoveM=TRUE, bRemoveRandom=TRUE,bSummarizedExperiment=FALSE)
ALLdba <- dba.count(ALLdba, bParallel=TRUE)
ALLdba <- dba.contrast(ALLdba,design = FALSE, group1 = ALLdba$masks$Leukemia)
ALLdba<-dba.analyze(ALLdba)

# Plot 
pdf("ALLvsWT_scatter_plot.pdf")
dba.plotMA(ALLdba,method=ALLdba$config$AnalysisMethod,contrast=1,th=0.05,bXY=TRUE, bNormalized=TRUE,bSmooth=FALSE)
dev.off()

pdf("ma_plot.pdf")
dba.plotMA(ALLdba,method=ALLdba$config$AnalysisMethod,contrast=1,th=0.05, bNormalized=TRUE,bSmooth=FALSE)
dev.off() 

# Save differentially accessible peaks
ALLdba.db3<-dba.report(ALLdba, contrast=1, method=ALLdba$config$AnalysisMethod, bNormalized=TRUE, bCalled=TRUE, bCounts=TRUE, bCalledDetail=TRUE, bAll=TRUE, initString='LEUvsWT',ext='csv',th=0.05, bUsePval=FALSE,DataType=DBA_DATA_FRAME)
write.table(ALLdba.db3,file="LEUvsWT.xls",sep="\t",row.names = FALSE)

# Identify lost and gained peaks by foldchange
leuGain <- ALLdba.db3 %>% filter(Fold >= 0)
leuLoss <- ALLdba.db3 %>% filter(Fold <= 0)

write.table(leuGain,file="B-ALL_gain_diffBind.bed",sep="\t", row.names = FALSE)
write.table(leuLoss,file="B-ALL_loss_diffBind.bed",sep="\t", row.names = FALSE)
