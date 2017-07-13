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
# R script to generate QC plots VcfToGenotypes
#
# The input should be the *.info.txt file produced by VcfToGenotypes and the
# output should be the pdf in which to store the plots.
#################################################################################

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)
library(scales)  # Require by ggplot

# Standard colors
fgblue  = "#2E9DD7"
fggreen = "#155936"
fgred   = "firebrick3"
fgcolors = c(fgblue, fggreen, fgred)

args       = commandArgs(trailingOnly=T)
inputFile  = args[1]
outputFile = args[2]

data = read.table(inputFile, sep="\t", header=F)

pdf(outputFile, width=11, height=8.5)

for (i in 1:nrow(data)) {
    row = data[i,]
    chromosome  = row[1]
    start       = row[2]
    end         = row[3]
    unitLength  = row[4]
    refNumUnits = row[5]
    strName     = row[[6]]
    known       = c()
    for (j in 7:ncol(row)) {
        if (grepl(":", row[[j]])) {
            stopifnot(j == ncol(row))
            rowData     = toString(row[[j]])
        }
        else {
            known = append(known, row[[j]])
        }
    }

    # Row data contains a comma-seperated list of key-value pairs (colon-delimited), each representing the repeat length
    # and observation counts
    numUnits = c()
    counts   = c()
    tuples = strsplit(rowData, ",")[[1]]
    for (j in 1:length(tuples)) {
        tuple = tuples[j]
        tupleData = strsplit(tuple, ":")[[1]];
        numUnits[j] = as.numeric(tupleData[1])
        counts[j]   = as.numeric(tupleData[2])
    }
    df = data.frame(numUnits, counts)
    p = ggplot(df) + aes(x=numUnits) +
        geom_point(aes(y=counts, color=counts)) +
        scale_x_continuous() +
        scale_fill_manual(values=fgcolors) +
        scale_color_gradient(low=fgblue, high=fggreen) +
        labs(x="# of Repeat Units", y="Raw counts", title=paste(strName, " (", chromosome, ":", start, "-", end, " motif of length ", unitLength, ")", sep="")) +
        theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())
    for (k in known) {
        p = p + geom_vline(xintercept = k, linetype="dotted")
    }
    print(p)
}

dev.off()