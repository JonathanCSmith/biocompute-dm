######
rm(list=ls())
library("ggplot2")
library("reshape")
library("Cairo")

# This Rscript is made to run with command line arguments.

# Arg 1 is the path to where the table with post alignment metrics is.
# Arg 1 example: "/home/<user>/<project>/qctable.csv"

# Arg 2 is the path to where the plots will be writen to.
# Arg 2 example: "/home/<user>/<project>/qcFigs/"

# it assumes certain names for the table columns, as these should match the QC metrics we decided.

# the result is the production of several PNG files which are the ones we had decided,


###############################################################
####### THIS IS THE TABLE THAT CONTAINS SAMPLE RUN INFO #######
###############################################################
args = commandArgs(trailingOnly = TRUE)
pathToTable = args[1]
outputOfPlots = args[2]
qcdf = read.csv(file = pathToTable)

print("friend, you cookin?")
print(head(qcdf))

# a multiplication to get the non-duplicates


## ^^ this dataframe is what could be a standard entry point for this type of data
# something like:
# qcdf = read.table(file="/home/user/blabla/blabla.csv")
##



####################################
######## NOW do plotings ###########
####################################
# i commented out the Cairo stuff that saves the image


### this puts the readParis<alignedReads<withoutDuplicates together (notice that multiple aligments are still there)
qcdf_melt=melt(qcdf[,c("sampleNames", "totalReads", "alignedReads", "uniquelyAlignReads", "uniquelyAlignReadsMinusdup")])
CairoPNG(file=paste(outputOfPlots, "readsOverall.png", sep=""), res=111, width=860*2.02, height=860)
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


qcdf_melt=melt(qcdf[,c("sampleNames", "fractionDuplicates")])
CairoPNG(file=paste(outputOfPlots, "fractionDuplicates.png", sep=""), res=111, width=860*2.02, height=860)
g = ggplot(data=qcdf_melt)
g = g + geom_bar(mapping=aes(x=sampleNames, y=value),
stat="identity",
position="identity",
alpha=0.88)
g = g + theme(axis.text.x = element_text(size  = 10,
angle = 45,
hjust = 1,
vjust = 1))
g = g + ggtitle("fraction of duplicates") + ylim(0,1.1) + ylab("fraction of duplicate reads")
g
dev.off()
