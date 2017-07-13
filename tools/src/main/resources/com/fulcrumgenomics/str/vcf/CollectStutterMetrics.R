#################################################################################
# The MIT License                                                               #
# Copyright (c) 2017 Fulcrum Genomics LLC                                       #
# Permission is hereby granted, free of charge, to any person obtaining a copy  #
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      #
# The above copyright notice and this permission notice shall be included in    #
# all copies or substantial portions of the Software.                           #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     #
# THE SOFTWARE.                                                                 #
#################################################################################

#################################################################################
# R script to generate QC plots CollectStutterMetrics
#
# Three inputs should be given.  The first two should be the *.stutter.tab and
# *.stutter.raw.tab, in that order, produced by CollectStutterMetrics.  Next,
# should be the pdf in which to store the plots.
#################################################################################


options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)
library(scales)  # Require by ggplot
library(reshape2)

# Standard colors
fgblue  = "#2E9DD7"
fggreen = "#155936"
fgred   = "firebrick3"
fpurple = "purple"

args          = commandArgs(trailingOnly=T)
consensusFile = args[1]
rawFile       = args[2]
outputFile    = args[3]


consensusData = read.table(consensusFile, sep="\t", header=F, stringsAsFactors=FALSE)
rawData       = read.table(rawFile, sep="\t", header=F, stringsAsFactors=FALSE)

pdf(outputFile, width=11, height=8.5)

stopifnot(nrow(consensusData) == nrow(rawData))
consensusX = as.numeric(consensusData[1,2:ncol(consensusData)])
rawX       = as.numeric(rawData[1,2:ncol(rawData)])
stopifnot(all(consensusX == rawX))

allConsensusY = as.numeric(consensusData[nrow(consensusData),2:ncol(consensusData)])
allRawY       = as.numeric(consensusData[nrow(consensusData),2:ncol(consensusData)])

# Normalize to 1
allConsensusY = allConsensusY / sum(allConsensusY)
allRawY       = allRawY / sum(allRawY)

# all but the last row ("all_strs")
for (i in 2:(nrow(consensusData)-1)) {
    name = consensusData[[i,1]]
    stopifnot(rawData[i,1] == name)

    consensusY = as.numeric(consensusData[i,2:ncol(consensusData)])
    rawY       = as.numeric(rawData[i,2:ncol(rawData)])

    # Normalize to 1
    consensusY = consensusY / sum(consensusY)
    rawY       = rawY / sum(rawY)

     df = data.frame(
        stutter      = consensusX,
        strRaw       = rawY,
        allRaw       = allRawY,
        strConsensus = consensusY,
        allConsensus = allConsensusY
    )

    df = melt(df, id.vars="stutter", value.name="counts")

    p = ggplot(df, aes(x=stutter, y=counts, fill=variable)) +
        geom_bar(width = 0.9, position=position_dodge(width = 0.9), stat="identity", colour="black") +
        scale_fill_manual(
            labels=c("STR raw reads", "All raw reads", "STR consensus calls", "All consensus calls"),
            values=c(fgblue, fggreen, fgred, fpurple)
        ) +
        labs(x="Stutter distance", y="Normalized counts", title=name) +
        theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())
    print(p)
}
dev.off()