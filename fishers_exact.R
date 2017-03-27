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

suppressPackageStartupMessages(library(optparse))

# --------------
# fishers_exact.R v.1.0		by Colette L Picard		05/10/2016

# Version history:
# --------------
# v.1.0: initial version		05/10/2016
# --------------

# This simple script will calculate Fisher's exact test for you using the R builtin fisher.test
# --------------

# Usage 1: provide a set of 4 numbers using --numbers:
# fishers_exact.R [options] --numbers xx,xy,yx,yy
# Example:
# fishers_exact.R --numbers 5,4,1,9
# Will create the following contingency table:
#			 col1	col2
# my_table = | 5	4 | row1
#			 | 1	9 | row2
# and will test the hypothesis that the rows and columns are independent (e.g. there is no significant
# difference between distribution of col1 and col2 in row1 vs row2 or vice versa). Note that
# fishers_exact.R --numbers 5,4,1,9 is the same as fishers_exact.R --numbers 5,1,4,9
#			 col1	col2
# my_table = | 5	1 | row1
#			 | 4	9 | row2

# --------------
# Usage 2: provide an input file and a set of 4 numbers representing columns to use as elements in
# the contingency table. fisher.test will be run for each row, and the pval reported in a new column.
# Once all tests run, the pvals will be corrected using the Benjamini-Hochberg method.
# fishers_exact.R --infile infile --numbers 2,3,4,5. 
#			 col1	col2
# my_table = | 2	3 | row1		<- here the numbers in the cells represent the column that the actual value would be taken from
#			 | 4	5 | row2

option_list <- list(
	make_option("--numbers", help = "Comma-separated list of 4 numbers to create contingency table (if --infile provided, these are column numbers)"),
	make_option("--infile", help = "Name of input file if not using --numbers (one Fisher's exact test per row)"),
	make_option("--header", default=FALSE, action="store_true", help = "Input file contains a header"),
	make_option("--outfile", help = "Name for output file"),
	make_option("--name", default="all", help = "Suffix for newly created variables with Fisher pvalues pval_* and BHpval_*")
	)

args <- commandArgs(trailingOnly = TRUE)	
if (length(args) <= 1) {
	cat("Usage: fishers_exact.R [options]\n")
	cat("----------------------\n")
	cat("--numbers : comma-separated list of 4 numbers to create contingency table (if --infile specified, these represent column numbers)\n")
	cat("--infile : Name of input file if not using --numbers (one Fisher's exact test per row)\n")
	cat("--header : input file contains a header\n")
	cat("--name : name for new column (will have prefix pval_*)\n")
	cat("--outfile : Name for output file\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}

arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = 0)
opt = arguments$options

if (length(opt$numbers) != 0 && length(opt$infile) == 0) {
	numbers = strsplit(opt$numbers,",")[[1]]
	if (length(numbers) != 4) {
		stop("list of --numbers must contain 4 elements, provided list has ",length(numbers))
	}
	mm = matrix(c(as.numeric(numbers[1]),as.numeric(numbers[3]),as.numeric(numbers[2]),as.numeric(numbers[4])), nrow=2)
	rownames(mm) = c("row1","row2")
	colnames(mm) = c("col1","col2")
	cat("Performing Fisher's Exact Test on the following contingency table:\n")
	print(mm)
	cat("H0: rows and columns are independent\n")
	
	fisher_res = fisher.test(mm)
	cat("p-value under H0 is",signif(fisher_res$p.value, 5),"\n")	
} else if (length(opt$infile) != 0) {
	if (length(opt$numbers) == 0) {
		stop("If providing input file, must use --numbers to indicate which columns to use for contingency table")
	}
	if (length(opt$outfile) == 0) {
		stop("If providing input file, must use --outfile to indicate where to save output")
	}
	
	cat("Performing fisher's exact test on the rows of input file:",opt$infile,"\n")
	cat("Columns to be used for contingency table:",opt$numbers,"\n")
	
	numbers = strsplit(opt$numbers,",")[[1]]
	if (length(numbers) != 4) {
		stop("list of --numbers must contain 4 elements, provided list has ",length(numbers))
	}

	alldata = read.table(opt$infile, header=opt$header, sep="\t", stringsAsFactors = TRUE)
	alldata[[paste("pval_",opt$name,sep="")]]=NA
	cat("Done reading in data. Performing tests...\n")
	for (i in 1:nrow(alldata))	{
		if (i %% 10000 == 0 || i == 1) {
			cat("Processing row #",i,"\n")
		}		
		contingency = matrix(c(alldata[i,as.numeric(numbers[1])], alldata[i,as.numeric(numbers[3])], alldata[i,as.numeric(numbers[2])], alldata[i,as.numeric(numbers[4])]),nrow=2)
		fisher_res = fisher.test(contingency)
		alldata[[paste("pval_",opt$name,sep="")]][i] = signif(fisher_res$p.value, 5)
	}
	# apply Benjamini-Hochberg correction
	alldata[[paste("BHpval_",opt$name,sep="")]] = p.adjust(alldata[[paste("pval_",opt$name,sep="")]], method="BH")
	
	write.table(alldata, opt$outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=opt$header, na="")	
}	
	
	
	
	
	
	
	
	
	
	
	
	

