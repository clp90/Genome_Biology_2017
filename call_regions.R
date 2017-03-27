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
   

library(optparse)

# version 1.2 (04/26/2016)
# -------------------------
# Version history:
# v.1.0: initial build
# v.1.1 (02/04/2015)
# 	- updates fixing corner cases involving the beginning and end of interval lists
# 	modified default parameter settings
# v.1.2 (04/26/2016)
#	- added "both" category to be able to identify heterozygous regions
# -------------------------

# This script accepts an input bed file with depth information from two samples binned in windows
# (5 columns: chromosome, start of window, end of window, count_sample1, count_sample2). The two
# samples are expected to represent two different strains (e.g. strain1 and strain2). The depth information
# will be used to infer regions that are homozygous for strain1, homozygous for strain2, heterozygous,
# or with too few total reads to determine. All windows in input file will be classified into these
# four categories. This script is designed to handle a single chromosome at a time.

# Usage: call_regions.R [options] infile.bed outprefix

# infile.bed = bed file with 5 fields: chr, start, end, depth_sample1 and depth_sample2)
# outprefix = prefix for output files

# This script will generate the following output files:
# TODO

option_list <- list(
	make_option("--smoothn", default=25,
		help = "For smoothing using adjacent windows, use this many bins to the left and right of current bin (e.g. for --windowsize 5, 5 bins to the left and 5 to the right will be used, for 11 bins total)"),
	make_option("--minFoldChange", default=1.5,
		help = "For a bin to be declared homozygous, must have this many times more reads than the other strain (a > minFoldChange*b)"),
	make_option("--mindif", default=2,
		help = "For a bin to be declared homozygous, must have this many more reads that the other strain (abs(a-b) > mindif"),
	make_option("--mindepth", default=5,
		help = "Bins covered with fewer than mindepth reads total in both strains are considered to have too few reads to categorize" ),
	make_option("--minregionsize", default=2000,
		help = "Minimum size of homozygous region" ),
	make_option("--strain1", default="strain1",
		help = "Name of strain representing column 4 of input file"),
	make_option("--strain2", default="strain2",
		help = "Name of strain representing column 5 of input file"),
	make_option("--maxmerge", default=50000,
		help = "Maximum size (in bps) for a region to be merged into its flanking regions")
	)

args <- commandArgs(trailingOnly = TRUE)	
if (length(args) <= 1) {
	cat("Usage: call_regions.R [options] infile.bed outprefix\n")
	cat("----------------------\n")
	cat("infile.bed contains depth information for two different strains in BED format (5 fields)\n")
	cat("outprefix is the prefix used for all output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--smoothn <int> : [default 25] use this many bins to the left and right of a bin when smoothing using a sliding window\n")
	cat("--minFoldChange <float> : [default 1.5] For a bin to be declared homozygous, must have this many times more reads than the other strain (a > minFoldChange*b)\n")
	cat("--mindif <float> : [default 2] For a bin to be declared homozygous, must have this many more reads that the other strain (abs(a-b) > mindif\n")
	cat("--mindepth <int> : [default 5] Bins covered with fewer than mindepth reads total in both strains are considered to have too few reads to categorize\n")
	cat("--minregionsize <int> : [default 2000] Minimum size of homozygous region\n")
	cat("--strain1 <str> : [default \"strain1\"] name of strain representing column 4 in input file\n")
	cat("--strain2 <str> : [default \"strain2\"] name of strain representing column 5 in input file\n")
	cat("--maxmerge <int> : [default 50000] small regions of less than maxmerge bp can be merged into larger regions if both flanking regions agree\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}
options(scipen=999)
arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = 2)
opt = arguments$options
infile = arguments$args[1]
outprefix = arguments$args[2]
smoothn = opt$smoothn
minFoldChange = opt$minFoldChange
maxmerge = opt$maxmerge
strain1 = opt$strain1
strain2 = opt$strain2
mindepth = opt$mindepth
mindif = opt$mindif
minregionsize = opt$minregionsize

# -------------------
# Helper functions
# -------------------
movingAverage <- function(x, n=1, centered=FALSE) {
# from R cookbook: http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
    if (centered) {
        before <- floor  ((n-1)/2)
        after  <- ceiling((n-1)/2)
    } else {
        before <- n-1
        after  <- 0
    }

    # Track the sum and count of number of non-NA items
    s     <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data 
    new <- x
    # Add to count list wherever there isn't a 
    count <- count + !is.na(new)
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new   <- c(rep(NA, i), x[1:(length(x)-i)])

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new   <- c(x[(i+1):length(x)], rep(NA, i))

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new

        i <- i+1
    }

    # return sum divided by count
    s/count
}

getRegions <- function(dataf,mark) {
# for a set of intervals marked by 0,1 in a data frame, groups together consecutive ones marked by 1
# into single interval. Use mark to add a numeric label for this particular dataset (e.g. strain1 = 1, strain2 = 2 etc.)
	cur_region = 1
	if (dataf$calls[1] == 1) {
		regions = cbind(dataf$start[1],0,mark)			# first region starts at position 0
		end = TRUE										# we are now looking for END of current region
	} else {
		regions = cbind()								# first region not yet found
		end = FALSE										# NOT looking for END of current region
	}
	
	diffs <- dataf$calls[-1L] != dataf$calls[-length(dataf$calls)]		# find all indices where 0 -> 1 or vice versa
	
	for (i in 2:length(dataf$calls)-1) {
		if (diffs[i] == TRUE & end == TRUE) {
			regions[cur_region,2] = dataf$start[i+1]
			end = FALSE
			cur_region = cur_region + 1
		} else { 
			if (diffs[i] == TRUE & end == FALSE) {
				regions = rbind(regions,cbind(dataf$start[i+1],0,mark))
				end = TRUE
			}
		}
	}
	# if after loop, end = TRUE, we are still waiting for the end of the interval, set to final pos
	if (end == TRUE) {
		regions[cur_region,2] = dataf$end[length(dataf$calls)]
	}
	
	regions
}	

mergeRegions <- function(curRegions, min_region_size) {
# given a set of regions, merges small regions flanked by larger regions of the
# same type (e.g. a small "none" flanked by strain2) into the larger region
# Note that only regions smaller than min_region_size will be reassigned
	curRegions_st = curRegions[order(curRegions[,1]),]		# sort regions
	
	# reassign smaller nested regions to the flanking region type	
	# for first interval, if it is small enough and flanking is large enough, merge
	if (curRegions_st[1,2] - curRegions_st[1,1] <= min_region_size && curRegions_st[2,2] - curRegions_st[2,1] > 100000) {
		curRegions_st[1,3] = curRegions_st[2,3]		# assign same label as 2nd interval
	}
	for (i in 2:(length(curRegions_st[,1])-1)) {
		# for each interval, confirm current left endpoint = previous right endpoint (no gaps)
		stopifnot(curRegions_st[i,1] == curRegions_st[i-1,2])	# check no gaps in intervals
		stopifnot(curRegions_st[i,3] != curRegions_st[i+1,3])	# check adjacent intervals not same label (strain1/strain2/unk)
				
		if (curRegions_st[i,2] - curRegions_st[i,1] <= min_region_size) {		# region is small enough (see min_region_size)
			if (curRegions_st[i-1,3] == curRegions_st[i+1,3]) {					
				curRegions_st[i,3] = curRegions_st[i-1,3]						# merge small region if flanking regions are same
			} else if (curRegions_st[i-1,2] - curRegions_st[i-1,1] > 100000) {	# merge to left region only if very large
				curRegions_st[i,3] = curRegions_st[i-1,3]	
			} else if (curRegions_st[i+1,2] - curRegions_st[i+1,1] > 100000) {	# merge to right region only if very large
				curRegions_st[i,3] = curRegions_st[i+1,3]	
			}				
		} 
	}
	# for last interval, if it is small enough and flanking is large enough, merge
	if (curRegions_st[length(curRegions_st[,1]),2] - curRegions_st[length(curRegions_st[,1]),1] <= min_region_size && curRegions_st[length(curRegions_st[,1]-1),2] - curRegions_st[length(curRegions_st[,1])-1,1] > 100000) {
		curRegions_st[length(curRegions_st[,1]),3] = curRegions_st[length(curRegions_st[,1])-1,3]		# assign same label as 2nd to last interval
	}
	
	curRegions_st
}

combineByType <- function(curRegions) {
# given a set of intervals in the form (start,end,type), combines all adjacent
# intervals with the same type (aka mark)
# returns a list of intervals in the same format, with no two adjacent intervals having the same type
	curRegions_st = curRegions[order(curRegions[,1]),]				# sort regions	
	mark = curRegions_st[1,3]										# label of 1st region
	i = 1															# index into cur pos in curRegions_st
	cur_region = 1													# index into cur pos in refined_regions
	refined_regions = cbind(curRegions_st[1,1],0,mark)				# set initial interval, use 0 as placemarker for end
	while (i < length(curRegions_st[,1])) {
#		message(paste("find all intervals w/ same mark as interval ",i))
		while (i <= length(curRegions_st[,1]) && curRegions_st[i,3] == mark) {
			i = i+1
#			message(paste("interval ",i-1))
		}
		# after this loop, i-1 is index of last interval w/ same mark
		refined_regions[cur_region,2] = curRegions_st[i-1,2]		# add to list of intervals
		# check if reached end of interval list; break if yes
		if (i == length(curRegions_st[,1])+1) {
			break
		}
		# otherwise, get mark of next interval, start new interval in refined_regions
		mark = curRegions_st[i,3]
		refined_regions = rbind(refined_regions,cbind(curRegions_st[i,1],0,mark))
		cur_region = cur_region+1
	}
	refined_regions[cur_region,2] = curRegions_st[length(curRegions_st[,1]),2]	# fill in final entry
	
	refined_regions = matrix(refined_regions,ncol=3)							# convert to matrix and return
	refined_regions
}

refineRegions <- function(curRegions, min_region_size) {
# given a set of regions, refines them by letting small regions flanked by larger regions of the
# same type (e.g. a small "none" flanked by strain2) be incorporated into the larger region
# Note that only regions smaller than min_region_size will be reassigned
	if (length(curRegions) == 3) {						# if only one region, return it
		return(curRegions)
	}	
	# otherwise, merge nearby small intervals into adjacent larger ones and combine
	merged = mergeRegions(curRegions, min_region_size)
	combined = combineByType(merged)
	combined
}

# -------------------
# MAIN
# -------------------
# output options for user
cat("\ncall_regions.R v.1.2 - 03/29/2016 by Colette L Picard\n")
cat("----------------------\n")
cat("Input file: ",infile,"\n")
cat("Prefix for output files: ",outprefix,"\n")
cat("Strain 1 name: ",strain1,"\n")
cat("Strain 2 name: ",strain2,"\n")
cat("----------------------\n")
cat("Additional parameters:\n")
cat("Bins to use on either side for smoothing: ",smoothn,"\n")
cat("Min depth for bin to be assigned to either strain1, strain2 or both: ",mindepth,"\n")
cat("Min fold change difference between sample1 and sample2 for bin to be considered homozygous: ",minFoldChange,"\n")
cat("Min absolute difference (in # reads) between sample1 and sample2 for bin to be considered homozygous: ",mindif,"\n")
cat("Max size (in bp) of unassigned region that can be merged into flanking assigned regions: ",maxmerge,"\n")
cat("Min region size for adjacent small unassigned regions to be merged into this bigger region: ",minregionsize,"\n")
cat("----------------------\n\n")

# -------------------
# Read in the depth file into a table and smooth using a symmetrical moving average
# -------------------
depth <- read.table(infile)
if (length(levels(depth$V1)) > 1) {
	stop("input files should contain data from a single chromosome")
}
chr = depth$V1[1]

# Perform smoothing using windowsize provided (default 21)
strain1_smooth <- movingAverage(depth$V4, smoothn*2 + 1, centered=TRUE) 
strain2_smooth <- movingAverage(depth$V5, smoothn*2 + 1, centered=TRUE) 

# output the smoothed datasets for viewing in IGV
write.table(data.frame(depth$V1,depth$V2,depth$V3,strain1_smooth), paste(outprefix,"_",strain1,"_smoothed_",chr,".bedgraph",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(data.frame(depth$V1,depth$V2,depth$V3,strain2_smooth), paste(outprefix,"_",strain2,"_smoothed_",chr,".bedgraph",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# -------------------
# Based on smoothed data, make initial call of most likely strain1 and strain2 regions
# -------------------
assn_both <- ifelse(abs(strain1_smooth - strain2_smooth) < mindif & (strain1_smooth+strain2_smooth > mindepth) ,1,0)
assn_strain1 <- ifelse(!assn_both & (strain1_smooth > minFoldChange*strain2_smooth) & (strain1_smooth+strain2_smooth > mindepth) ,1,0)
assn_strain2 <- ifelse(!assn_both & (strain2_smooth > minFoldChange*strain1_smooth) & (strain1_smooth+strain2_smooth > mindepth) ,1,0)
assn_none <- ifelse(!assn_both & !assn_strain1 & !assn_strain2,1,0)

# get blocks of regions with same assignment
dataf = data.frame(start = depth$V2, end = depth$V3, calls = assn_strain1)
strain1_regions = getRegions(dataf,1)
dataf = data.frame(start = depth$V2, end = depth$V3, calls = assn_strain2)
strain2_regions = getRegions(dataf,2)
dataf = data.frame(start = depth$V2, end = depth$V3, calls = assn_none)
none_regions = getRegions(dataf,0)
dataf = data.frame(start = depth$V2, end = depth$V3, calls = assn_both)
both_regions = getRegions(dataf,3)

#COMwrite.table(cbind(as.character(chr),strain1_regions), paste(outprefix,"_testStrain1.bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#COMwrite.table(cbind(as.character(chr),strain2_regions), paste(outprefix,"_testStrain2.bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#COMwrite.table(cbind(as.character(chr),both_regions), paste(outprefix,"_testboth.bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# -------------------
# Now refine initial calls by iteratively merging small windows into bigger flanking regions
# -------------------
i = 1
message(paste("Refining intervals:"))
all_regions = rbind(strain1_regions,strain2_regions,none_regions,both_regions)
regions = combineByType(all_regions)		# make sure no adjacent regions of same type
message(paste("Iteration",i,"..."))
refined = refineRegions(regions,minregionsize)

while (!identical(refined,regions)) {
	# repeat until convergence...
	if (length(refined) == 3) {
		break		# only one line left, no more to do
	}
	# otherwise, refine the regions
	i = i+1
	message(paste("Iteration",i,"..."))
	regions = refined
	refined = refineRegions(regions, maxmerge)
}


# -------------------
# Output regions as bed files - split by 1,2,0 (strain1,strain2,unk)
# -------------------
strain1_regions_final = refined[refined[,3] == 1,] 
strain2_regions_final = refined[refined[,3] == 2,] 
het_regions_final = refined[refined[,3] == 3,] 
unk_regions_final = refined[refined[,3] == 0,] 

if (length(strain1_regions_final) == 3 & length(strain1_regions_final) != 0) {
	out = matrix(strain1_regions_final,ncol=3)
	write.table(cbind(as.character(chr),out[,1],out[,2]), paste(outprefix,"_",strain1,"_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else if (length(strain1_regions_final) != 0) {
	write.table(cbind(as.character(chr),strain1_regions_final[,1],strain1_regions_final[,2]), paste(outprefix,"_",strain1,"_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else {
	cat("No regions for",strain1,"identified\n")
}

if (length(strain2_regions_final) == 3 & length(strain2_regions_final) != 0) {
	out = matrix(strain2_regions_final,ncol=3)
	write.table(cbind(as.character(chr),out[,1],out[,2]), paste(outprefix,"_",strain2,"_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else if (length(strain2_regions_final) != 0) {
	write.table(cbind(as.character(chr),strain2_regions_final[,1],strain2_regions_final[,2]), paste(outprefix,"_",strain2,"_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else {
	cat("No regions for",strain2,"identified\n")
}

if (length(unk_regions_final) == 3 & length(unk_regions_final) != 0) {
	out = matrix(unk_regions_final,ncol=3)
	write.table(cbind(as.character(chr),out[,1],out[,2]), paste(outprefix,"_unk_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else if (length(unk_regions_final) != 0) {
	write.table(cbind(as.character(chr),unk_regions_final[,1],unk_regions_final[,2]), paste(outprefix,"_unk_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else {
	cat("No unknown regions identified\n")
}

if (length(het_regions_final) == 3 & length(het_regions_final) != 0) {
	out = matrix(het_regions_final,ncol=3)
	write.table(cbind(as.character(chr),out[,1],out[,2]), paste(outprefix,"_het_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else if (length(het_regions_final) != 0) {
	write.table(cbind(as.character(chr),het_regions_final[,1],het_regions_final[,2]), paste(outprefix,"_het_regions_final_",chr,".bed",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
} else {
	cat("No heterozygous regions identified\n")
}




