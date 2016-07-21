######
#rm(list=ls())
library("ggplot2")
library("reshape")
library("Cairo")

args = commandArgs(trailingOnly = TRUE)
pathToTable = args[1]
outputOfPlots = args[2]
qcdf = read.csv(file = pathToTable)


### Total Alignment metrics

CairoPNG(file=paste(outputOfPlots, "alignment_metrics.png", sep=""), res=111, width=860*2.02, height=860)#
qcdf_melt=melt(qcdf[,c("sampleNames", "totalReads", "alignedReads", "uniquelyAlignReads", "uniquelyAlignReadsMinusdup")])
g = ggplot(data=qcdf_melt)
g = g + geom_bar(mapping=aes(x=sampleNames, y=value, fill=variable),
                 stat="identity",
                 position="identity",
                 alpha=0.88)
g = g + theme(axis.text.x = element_text(size  = 10,
                                         angle = 90,
                                         hjust = 1,
                                         vjust = 1))
g = g + ggtitle("reads aligned") + ylab("number of reads")
g
dev.off()


### this shows duplication fractions

CairoPNG(file=paste(outputOfPlots, "percentDuplicates.png", sep=""), res=111, width=860*2.02, height=860)
qcdf_melt=melt(qcdf[,c("sampleNames", "percentDuplicates")])
g = ggplot(data=qcdf_melt)
g = g + geom_bar(mapping=aes(x=sampleNames, y=value),
                 stat="identity",
                 position="identity",
                 alpha=0.88)
g = g + theme(axis.text.x = element_text(size  = 10,
                                         angle = 45,
                                         hjust = 1,
                                         vjust = 1))
g = g + ggtitle("fraction of duplicates") + ylim(0,1.1) + ylab("percent of duplicate reads")
g
dev.off()
