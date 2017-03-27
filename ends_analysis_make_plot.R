#!/usr/bin/env Rscript

# Copyright 2017 Colette L Picard

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

# version 1 (05/9/2015)
# helper script for ends_analysis.sh, makes the ends-analysis plot.
# Needs to know the numIn, numOut and width values from ends_analysis.sh.

# v.1.1 (08/17/2016)
#	- small update to add option to output plot in PDF format

# Usage: ends_analysis_make_plot.R infile outprefix numIn numOut width

args = commandArgs(trailingOnly = TRUE)
if (length(args) <= 6) {
	stop("Error: 7 arguments required. \nUsage: ends_analysis_make_plot.R infile yupper ylabel title outprefix numIn numOut width [options]")
}

option_list <- list(
	make_option("--colors", default="",
		help = "Colors for samples in plot in order"),
	make_option("--makebarchart", default=FALSE, action="store_true",
		help = "Make a barchart instead of a line plot"),
	make_option("--PDF", default=FALSE, action="store_true",
		help = "Output plot as PDF instead of PNG")
	)

arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = 9)
opt = arguments$options
infile = arguments$args[1]
yupper = arguments$args[2]
ylabel = arguments$args[3]
ftitle = arguments$args[4]
outprefix = arguments$args[5]
numIn = as.numeric(arguments$args[6])
numOut = as.numeric(arguments$args[7])
width = as.numeric(arguments$args[8])
linewidth = as.numeric(arguments$args[9])
colors = opt$colors
makebarchart = opt$makebarchart

counts = read.table(infile, header = FALSE, sep = '\t')
df_counts = data.frame(counts)
if (colors != "") {
	colors = strsplit(colors," ")[[1]]
}

# summary of calls to this script
cat("Counts file is:",infile,"\n")
names(df_counts)[names(df_counts)=="V1"] <- "binID"
names(df_counts)[names(df_counts)=="V2"] <- "counts"
names(df_counts)[names(df_counts)=="V3"] <- "sample"

upstream_out = numOut/width
middle = (numIn/width)+(numOut/width)
downstream_in = (numIn/width) + middle
downstream_out = numOut/width + downstream_in

if (numIn == 0) {			# corner case
	cat("bins 5' upstream outside of feature: 1 to",upstream_out[1],"\n")
	cat("bins 5' upstream inside of feature:",upstream_out[1],"to",middle[1],"\n")
	cat("bins 3' downstream inside of feature:",middle[1],"to",downstream_in[1],"\n")
	cat("bins 3' downstream outside of feature:",downstream_in[1]+1,"to",downstream_out[1],"\n")
} else {
	cat("bins 5' upstream outside of feature: 1 to",upstream_out[1],"\n")
	cat("bins 5' upstream inside of feature:",upstream_out[1]+1,"to",middle[1],"\n")
	cat("bins 3' downstream inside of feature:",middle[1]+1,"to",downstream_in[1],"\n")
	cat("bins 3' downstream outside of feature:",downstream_in[1]+1,"to",downstream_out[1],"\n")
}

factororder = unique(df_counts$sample)
df_counts$sample = factor(df_counts$sample, levels = factororder)

# make x values for a "gap" to represent middle of feature
if (numIn == 0) {			# corner case
	gaplen = 0
} else {
	gaplen = floor(length(unique(df_counts$binID))/8)
	if (gaplen == 0) { gaplen = 0.5 }
}
pregap = subset(df_counts, binID <= middle)
postgap = subset(df_counts, binID > middle)
postgap$binID = postgap$binID+gaplen

if(numOut < 1000) { outKB = paste(numOut,"bp",sep="") } else { outKB = paste(numOut/1000,"kb",sep="") }
if(numIn < 1000) { inKB = paste(numIn,"bp",sep="") } else { inKB = paste(numIn/1000,"kb",sep="") }

# output plot
if (opt$PDF == TRUE) {
	pdf(paste(outprefix,'_plot.pdf',sep=""), width = 12, height = 6.5, useDingbats=FALSE)
} else {
	png(paste(outprefix,'_plot.png',sep=""), width = 12, height = 6.5, units = 'in', res = 300)
}
if (makebarchart == TRUE) {
	a = ggplot() + 
		geom_bar(data = pregap, aes(x=binID, y=counts, fill=sample), stat="identity", position=position_dodge(), colour="black") +
		geom_bar(data = postgap, aes(x=binID, y=counts, fill=sample), stat="identity", position=position_dodge(), colour="black") +
		xlab("") + ylab(ylabel) +
		scale_x_continuous(breaks=c(1,upstream_out + 0.5,middle,middle+gaplen+1,downstream_in + gaplen + 0.5,downstream_out+gaplen), labels=c(paste("-",outKB,sep=""),"feature\nstart 5'",paste("+",inKB,sep=""),paste("-",inKB,sep=""),"feature\nend 3'",paste("+",outKB,sep=""))) + 
		theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle(ftitle) +
		geom_vline(xintercept=upstream_out + 0.5, linetype="dotted") +
		geom_vline(xintercept=downstream_in + gaplen + 0.5, linetype="dotted") +
		theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
		theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
	if (length(colors) != 1) {
		a = a + scale_fill_manual(name = "Sample", values = colors)
	}
	if (yupper != "NULL") {
		a = a + expand_limits(x = 1, y = c(0,as.numeric(yupper)))
	} else {
		a = a + expand_limits(x = 1, y = 0)
	}	
	print(a)
} else {
	yupper = as.numeric(yupper)
	if (length(colors) != 1) {
		if (args[2] != "NULL") {
			cat("Outputting plot with custom color scheme and y axis limits\n")
			a = ggplot(pregap, aes(x=binID, y=counts, colour=sample)) +
				geom_line(size=linewidth) +
				geom_line(data=postgap, aes(x=binID, y=counts),size=linewidth) +
				theme_bw() + theme(panel.grid.minor = element_blank()) + 
				scale_x_continuous(breaks=c(1,upstream_out + 0.5,middle,middle+gaplen+1,downstream_in + gaplen + 0.5,downstream_out+gaplen), labels=c(paste("-",outKB,sep=""),"feature\nstart 5'",paste("+",inKB,sep=""),paste("-",inKB,sep=""),"feature\nend 3'",paste("+",outKB,sep=""))) +
				xlab("") + ylab(ylabel) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
				theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
				ggtitle(ftitle) + ylim(0, yupper) +
				scale_color_manual("", values=colors)
			if (numIn != 0) {
				a = a + geom_vline(xintercept=upstream_out + 0.5, linetype="dotted") + geom_vline(xintercept=downstream_in + gaplen + 0.5, linetype="dotted")
			}
			print(a)
		} else {
			cat("Outputting plot with custom color scheme and default y axis\n")
			a = ggplot(pregap, aes(x=binID, y=counts, colour=sample)) +
				geom_line(size=linewidth) +
				geom_line(data=postgap, aes(x=binID, y=counts),size=linewidth) +
				theme_bw() + theme(panel.grid.minor = element_blank()) + 
				scale_x_continuous(breaks=c(1,upstream_out + 0.5,middle,middle+gaplen+1,downstream_in + gaplen + 0.5,downstream_out+gaplen), labels=c(paste("-",outKB,sep=""),"feature\nstart 5'",paste("+",inKB,sep=""),paste("-",inKB,sep=""),"feature\nend 3'",paste("+",outKB,sep=""))) +
				xlab("") + ylab(ylabel) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
				theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
				ggtitle(ftitle) + scale_color_manual("", values=colors)
			if (numIn != 0) {
				a = a + geom_vline(xintercept=upstream_out + 0.5, linetype="dotted") + geom_vline(xintercept=downstream_in + gaplen + 0.5, linetype="dotted")
			}
			print(a)
		}
	} else {
		if (args[2] != "NULL") {
			cat("Outputting plot with default color scheme and custom y axis\n")
			a = ggplot(pregap, aes(x=binID, y=counts, colour=sample)) +
				geom_line(size=linewidth) +
				geom_line(data=postgap, aes(x=binID, y=counts),size=linewidth) +
				theme_bw() + theme(panel.grid.minor = element_blank()) + 
				scale_x_continuous(breaks=c(1,upstream_out + 0.5,middle,middle+gaplen+1,downstream_in + gaplen + 0.5,downstream_out+gaplen), labels=c(paste("-",outKB,sep=""),"feature\nstart 5'",paste("+",inKB,sep=""),paste("-",inKB,sep=""),"feature\nend 3'",paste("+",outKB,sep=""))) +
				xlab("") + ylab(ylabel) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
				theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
				ggtitle(ftitle) + ylim(0, yupper)
			if (numIn != 0) {
				a = a + geom_vline(xintercept=upstream_out + 0.5, linetype="dotted") + geom_vline(xintercept=downstream_in + gaplen + 0.5, linetype="dotted")
			}
			print(a)
		} else {
			cat("Outputting plot with default color scheme and default y axis\n")
			a = ggplot(pregap, aes(x=binID, y=counts, colour=sample)) +
				geom_line(size=linewidth) +
				geom_line(data=postgap, aes(x=binID, y=counts),size=linewidth) +
				theme_bw() + theme(panel.grid.minor = element_blank()) + 
				scale_x_continuous(breaks=c(1,upstream_out + 0.5,middle,middle+gaplen+1,downstream_in + gaplen + 0.5,downstream_out+gaplen), labels=c(paste("-",outKB,sep=""),"feature\nstart 5'",paste("+",inKB,sep=""),paste("-",inKB,sep=""),"feature\nend 3'",paste("+",outKB,sep=""))) +
				xlab("") + ylab(ylabel) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
				theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
				ggtitle(ftitle)
			if (numIn != 0) {
				a = a + geom_vline(xintercept=upstream_out + 0.5, linetype="dotted") + geom_vline(xintercept=downstream_in + gaplen + 0.5, linetype="dotted")
			}
			print(a)
		}
	}
}
dev.off()





