#!/bin/bash

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

# ------------------------------------------------------------------------------------
# v1.7 by Colette L. Picard
# 03/27/2017
# ------------------------------------------------------------------------------------

# Usage:
# ends_analysis.sh [options] -i reads1.bed,reads2.bed,...,readsN.bed -r regions.bed -o outprefix

# -------------------------
# Version history:
# v.1.0: initial build - 04/13/2015
# v.1.1: added option to do analysis based on mean value in VALUE field of .bed files
# (fifth field, e.g. BED file format is: chr start end name value strand). This allows
# methylation-based analyses. Will now output an intermediate file for each sample
# containing values for each bin in each feature. Format of this
# intermediate file has been changed to allow clustering (e.g. via Cluster 3.0).
# User can now request kmeans clustering on the file described above for each sample
# and generates ends analysis plots for each cluster. Uses Cluster 3.0 for clustering.
# The cluster feature has not been well tested and can be considered incomplete.
# v.1.2: added "verbose" option -v for debugging purposes. If -v is set, script will
# not delete any intermediate files.
# v.1.3: 10/24/2015
#	- changed plot width in ends_analysis_make_plot.R
#	- removed horizontal and vertical bars on plot
#	- greatly improved speed by running jobs simultaneously on cluster where possible
# v.1.4: 11/02/2015
#	- added -P option to parallelize job as described in 1.3, but won't by default (avoids spamming the 
#	cluster with too many jobs simultaneously when running many instances of ends_analysis.sh)
# v.1.5: 11/13/2015
#	- added functionality that lets users plot a single methylation profile over multiple
# types of regions instead of multiple methylation profiles over a single type of region
#	- to do this, provide a comma-separated list of regions to -r and a single file to -i
# (script will throw error if both -r and -i contain more than one item)
# 	- if multiple regions are detected (and a single -i), names and colors provided by -n and -c
# are expected to correspond to the regions files, instead of the -i files. By default, if only
# one file is provided to both -i and -r, names/colors are expected to correspond to -i.
# v.1.6: 04/21/2016
#	- added option -U to allow switching between weighted and unweighted averages across features
# 	- added option to censor points in plot where fewer than X genes contributed to average value
# 	- added -B option to plot barchart instead of line plot (don't use with small binsize! try -w 1000)
# v.1.7: 03/27/2017
#	- added -l option to modify thickness of lines in plot
# -------------------------

# Plots average coverage of reads in each input file seperately across the start and
# end of the regions in regions.bed. User can specify how many basepairs before and 
# after start. Example below is for 3bp outside of feature and 3bp into feature, bin width 1.
# (Note actual plot will be line not scatter or bar).

# 			1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22
#			-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
# regions:								[			  a region				]
# 
# reads:		1	1	1	1		2	2	2	2							3	3	3	3	3
#						4	4	4	4	4		5	5	5	5	5		
#									6	6	6	6	6			7	7	7	7	7	7
#
#
# result:
#									o	o		o							
#							o		o	o	o	o							o	o	o	
#							o	o	o	o	o	o					o	o	o	o	o	o
#						----------------------------		    ------------------------------
#										|									| 
#									  start								   end

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
ends_analysis 	v.1.7		03/27/2017
Plots average coverage of the provided input reads across the start and end of the 
provided regions. Number of bases to go outside the feature and into the feature
can be specified by the user (both default to 1kb). See options -O and -I below.
Calculates coverage seperately over each provided input file. Can also choose to output
an intermediate file which contains a matrix, where every row is a feature (e.g. gene)
and every column is a bin in the feature. This intermediate file is useful for downstream
analysis, including clustering - generate it using the -C option.

New in v.1.5: instead of plotting multiple methylation datasets on a single type of region (e.g. genes),
users can now request that a single methylation dataset be plotted on the same plot for multiple region
types (e.g. genes and TEs). If using this option, provide a single file to -i and multiple files to
-r. Names and colors provided with -n and -c will be expected to correspond to the regions files in
-r instead of the -i file. User cannot provide more than one file in both -i and -r simultaneously.

Example for -O 3, -I 5 and -w 1:
                  
      feature:        5'----------------------3'
      plot range:   ---[-----            -----]---

Takes strandedness into account, so "start" on left of plot is always the 5' end
of the feature, and "end" is always the 3' end (if no stranded info provided for features, assumes +).
					
Usage:
ends_analysis.sh [options] -i medata1.bed,medata2.bed,...,medataN.bed -r regions.bed -o outprefix
or
ends_analysis.sh [options] -i medata.bed -r regions1.bed,regions2.bed,...,regionsN.bed -o outprefix

User-specified options:
Required arguments:
	-i inlist : single or comma-separated list of BED files containing read or methylation data
	-r regions : single or comma-separated list of BED files containing regions to calculate average depth over
	-o outprefix : prefix for output files (including path)
Additional arguments:
	-s : path to folder containing all required helper scripts (see list below) - default scriptdir
	-n : list of names corresponding to the -i input file list - default use filename
	-O : number of bases outside of feature to plot at each end (see diagram above) - default 1000
	-I : number of bases inside of feature to plot at each end (see diagram above) - default 1000
	-w : width of the windows that average depth will be calculated over - default 25
	-N : normalization factors for each infile, comma-seperated and in same order as -i - default RPM normalization, ignored if -V is set
	-u : set upper limit for y axis manually in output plot - default does not set the upper limit
	-c : list of colors for plotting each sample, can be any R color format (including hex)
	-t : title for output plot - default "Average depth across features"
	-y : label for y axis of plot (useful when -V specified for describing the VALUE field of the BED file) - default "Average depth (normalized counts)"
	-V : calculate average of value in this field of BED files, instead of average depth (e.g. -V 4 will average 4th field)
	-m : if bin average is calculated over fewer than this many separate regions in -r file, omit from plot [0]
	-l : width of lines in plot (default 1)
Options:
	-B : turns on barchart mode - instead of plotting ends analysis line plot, will plot barchart (suggested bin width for -I = -O = 1000 is -w 1000, 4 bins total)
	-U : calculate unweighted average over all genes with data in bin, instead of weighting average so that genes with more data impact mean more (default weighted avg)
	-M : output matrix of values (one for each sample) with rows = regions and columns = window within that region
	-P : parallelize job where possible by submitting independent tasks to LSF cluster
	-R : allow overwrite of existing output files
	-v : verbose mode - useful for debugging purposes. Script will not delete any intermediate files.
	-h : prints this version and usage information
	
Required helper scripts (if not in same dir as this script, specify location with -s):
	- ends_analysis_get_intervals.py - by Colette L Picard
	- ends_analysis_process_intersect.py - by Colette L Picard
	- ends_analysis_make_plot.R - by Colette L Picard
	- ends_analysis_make_matrix.py - by Colette L Picard (only required if using -M or -C options)

Required installed on user PATH:
	- bedtools

		
------------------------------------------------------------------------------------
EOF

scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 

# Initiate environment
# ----------------------
# Required arguments:
# ----------------------
infilelist=""									# list of input files (reads in BED format)
regionlist=""									# list of regions to calculate average depth over
outprefix=""									# prefix for output files
# Set additional options to default values:
# ----------------------
path_to_scripts="$scriptDir"				# path to required scripts (default is current script's dir)
namelist=""									# list of names to use for the input files in plots, etc.
numOut=1000									# number of bases outside features to calculate over
numIn=1000									# number of bases inside features
width=25									# width of windows
yupper="NULL"								# upper limit of y axis in plot (by default, will use R default)
ylabel="Average depth (normalized counts)"	# label for y axis
norm=""										# list of factors for normalizing each library
colors=""									# list of colors for plot of each sample
title="Average depth across features"		# title for output plot
useValue=0									# calculate average over this field (if useValue = 0, uses depth instead)
minN=0										# don't plot point for this bin if fewer than this many genes/features were used to calc. it
linewidth=1
Verbose=false
writemat=false
overwrite=false
parallel=false
unweighted=false
barchart=false

# ----------------------
while getopts "i:r:o:s:n:O:I:w:y:t:u:N:c:V:m:l:vRBMPUh" opt; do
	case $opt in
		i)	# input file list
			infilelist="$OPTARG"
			;;
		r)	# regions file
			regionlist="$OPTARG"
			;;
		o)	# output directory
			outprefix="$OPTARG"
			;;
		s)	# path to scripts
			path_to_scripts="$OPTARG"
			;;
		n)	# name list
			namelist="$OPTARG"
			;;
		O)	# num bases out of feature
			numOut="$OPTARG"
			;;
		I)	# num bases in feature
			numIn="$OPTARG"
			;;
		w)	# width of windows
			width="$OPTARG"
			;;
		t)	# plot title
			title="$OPTARG"
			;;
		y)	# y axis title
			ylabel="$OPTARG"
			;;
		N)	# norm factors
			norm="$OPTARG"
			;;
		u)	# yupper
			yupper="$OPTARG"
			;;
		c)	# sample colors
			colors="$OPTARG"
			;;
		m)	# min bin N
			minN="$OPTARG"
			;;
		R)	# allow overwrite of output files
			overwrite=true
			;;
		M)	# output matrices for each sample
			writemat=true
			;;
		V)	# use VALUE field of BED instead of depth
			useValue="$OPTARG"
			;;
		l)	# line width in plot
			linewidth="$OPTARG"
			;;
		v)	# keep all intermediate files
			Verbose=true
			;;
		U)	# weight all genes equally across bin
			unweighted=true
			;;
		P)	# parallelize
			parallel=true
			;;
		B)	# barchart
			barchart=true
			;;
		h)	# print usage and version information to stdout and exit
			echo "$usage"
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# ----------------------
# Check that all files, scripts etc. required by this pipeline can be located
# ----------------------

# Check all required inputs are provided, if none then print help info
# ----------------------
if [[ $# -eq 0 ]]; then printf "%s\n" "$usage"; exit 0
elif [ -z "$infilelist" ]; then echo "Error: must provide a BED file or comma-seperated list of BED files with -i"; exit 1
elif [ -z "$regionlist" ]; then echo "Error: must provide a BED file with a list of regions with -r" "$log"; exit 1
elif [ -z "$outprefix" ]; then echo "Error: must provide a prefix for output files using -o"; exit 1; fi

# Check that all required helper scripts can be located
# ----------------------
if [ ! -f "$path_to_scripts"/ends_analysis_get_intervals.py ]; then
	echo "Error: required helper script $path_to_scripts/ends_analysis_get_intervals.py could not be located, see option -s"; exit 1
elif [ ! -f "$path_to_scripts"/ends_analysis_process_intersect.py ]; then
	echo "Error: required helper script $path_to_scripts/ends_analysis_process_intersect.py could not be located, see option -s"; exit 1
elif [ ! -f "$path_to_scripts"/ends_analysis_make_plot.R ]; then
	echo "Error: required helper script $path_to_scripts/ends_analysis_make_plot.R could not be located, see option -s"; exit 1
fi
if "$writemat"; then
	if [ ! -f "$path_to_scripts"/ends_analysis_make_matrix.py ]; then
		echo "Error: required helper script $path_to_scripts/ends_analysis_make_matrix.py could not be located, see option -s"; exit 1; fi; fi

# Check that all programs required on PATH are installed
# ----------------------
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required but was not found on PATH"; exit 1; }

# Check that window size is ok
# ----------------------
if [ ! $(( $numOut % $width )) -eq 0 ]; then echo "Error: outside feature length $numOut must be multiple of width $width"; exit 1; fi
if [ ! $(( $numIn % $width )) -eq 0 ]; then echo "Error: insize feature length $numIn must be multiple of width $width"; exit 1; fi

# If overwrite is not permitted, check no output files exist
# ----------------------
if [ "$overwrite" = "false" ]; then
	if [ -f "${outprefix}_plot.png" ]; then echo "Error: output file ${outprefix}_plot.png already exists, allow overwrite with -R"; exit 1; fi
	if [ -f "${outprefix}_average_counts.txt" ]; then echo "Error: output file ${outprefix}_average_counts.txt already exists, allow overwrite with -R"; exit 1; fi
else
	echo "Overwriting existing output files..."
	if ! "$Verbose"; then rm -f "${outprefix}_plot.png"; fi
	if ! "$Verbose"; then rm -f "${outprefix}_average_counts.txt"; fi
fi

# Check that all provided input files exist
# ----------------------
infiles=$(echo $infilelist | tr "," " ")
regions=$(echo $regionlist | tr "," " ")
numinfiles=0; numregions=0
for infile in $infiles; do
	if [ "${infile##*.}" != "bed" ]; then 
		if [ "${infile##*.}" != "bedgraph" ]; then
			echo "Error: input file $infile was expected to be in .bed or .bedgraph format but has extension .${infile##*.}"; exit 1; fi; fi
	if [ ! -f "$infile" ]; then echo "Error: input file $infile could not be opened"; exit 1; fi
	numinfiles=$(( $numinfiles + 1 ))
done
for region in $regions; do
	if [ "${regions##*.}" != "bed" ]; then 
		if [ "${regions##*.}" != "bedgraph" ]; then
			echo "Error: input file $region was expected to be in .bed or .bedgraph format but has extension .${region##*.}"; exit 1; fi; fi
	if [ ! -f "$region" ]; then echo "Error: region file $region could not be opened"; exit 1; fi
	numregions=$(( $numregions + 1 ))
done

# numregions and numinfiles are not allowed to both be >1
if [[ $numinfiles -gt 1 && $numregions -gt 1 ]]; then
	echo "Error: cannot provide more than one file in both -i and -r simultaneously"; exit 1
elif [[ $numregions -gt 1 ]]; then
	mode="regions"
else
	mode="infiles"
fi

# Check if names were provided, if not use input file stub names
# ----------------------
[ "$mode" == "regions" ] && filestouse=$regions || filestouse=$infiles
if [ -z "$namelist" ]; then
	# no name list provided, go through each file, add stub to namelist
	for f in $filestouse; do
		stub=${f##*/}; stub=${stub%.*}
		[ "$names" == "" ] && names="$stub" || names="$names $stub"
	done
else
	# make  sure length of list here == number of input files
	names=$(echo $namelist | tr "," " ")
	numnames=0
	for name in $names; do
		numnames=$(( $numnames + 1 ))
	done
	if [ "$mode" == "regions" ]; then
		[ ! $numnames -eq $numregions ] && { echo "Error: number of provided names (${numnames}) and number of provided region files (${numregions}) not the same"; exit 1; }
	else
		[ ! $numnames -eq $numinfiles ] && { echo "Error: number of provided names (${numnames}) and number of provided input files (${numinfiles}) not the same"; exit 1; }
	fi	
fi

# Get normalization factor for each sample, either provided in $norm or RPM, if using depth
# Skip if using -V
# ----------------------
if [[ "$mode" == "infiles" && "$useValue" == "0" ]]; then
	tot_reads=""
	for infile in $infiles; do
		if [ -z "$tot_reads" ]; then
			tot_reads=$( wc -l $infile | awk {'print $1'} )
		else
			tot_reads="$tot_reads $( wc -l $infile | awk {'print $1'} )"
		fi
	done
	
	if [ -z "$norm" ]; then
		for treads in $tot_reads; do
			if [ -z "$norm" ]; then
				norm=$(echo "scale=3; 1000000 / $treads" | bc)
			else
				norm="$norm $(echo "scale=10; 1000000 / $treads " | bc)"
			fi
		done
	else
		# make sure length of norm list == number of input files
		norm=$(echo $norm | tr "," " ")
		numnorms=0
		for n in $norm; do
			numnorms=$(( $numnorms + 1 ))
		done
		if [ ! $numnorms -eq $numinfiles ]; then
			echo "Error: number of provided normalization factors (${numnorms}) and number of provided input files (${numinfiles}) not the same"; exit 1; fi
	fi
fi


# Regions file must have this format
# ----------------------
for region in $regions; do
	cols_in_regions=$( awk -F$'\t' '{print NF}' "$region" | sort -nu | tail -n 1 )
	[[ "$cols_in_regions" != "6" && "$cols_in_regions" != "4" ]] && { echo "Error: regions file must be in format chr, start, end, name, with optional score, strand - 6 fields total if strand is important, otherwise 4"; exit 1; }
done

# If using -V, make sure it's not one of the first 3 fields, since those are chr/start/end
# ----------------------
if [[ "$useValue" -le 3 && "$useValue" != 0 ]]; then
	echo "Error: first three fields in .bed files are chr, start, and end - cannot use one of these as value field for -V option"; exit 1; fi

# Get colors for each sample
# ----------------------
if [ ! -z "$colors" ]; then
	colors=$(echo $colors | tr "," " ")
	numcol=0
	for n in $colors; do
		numcol=$(( $numcol + 1 ))
	done
	if [ "$mode" == "regions" ]; then
		[ ! $numcol -eq $numregions ] && { echo "Error: number of provided colors (${numcol}) and number of provided regions files (${numregions}) not the same"; exit 1; }
	else
		[ ! $numcol -eq $numinfiles ] && { echo "Error: number of provided colors (${numcol}) and number of provided regions files (${numinfiles}) not the same"; exit 1; }
	fi	
fi	

# ----------------------
# Helper functions for this script
# ----------------------
# prints an error message both to stdout and to the log file, then exits
# usage: err_msg msg logfile
err_msg ()
{
	printf "Error: $1 \n"
	exit 1	
}

displaytime () {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  [[ $D > 0 ]] && printf '%d days ' $D
  [[ $H > 0 ]] && printf '%d hours ' $H
  [[ $M > 0 ]] && printf '%d minutes ' $M
  [[ $D > 0 || $H > 0 || $M > 0 ]] && printf 'and '
  printf '%d seconds\n' $S
}

# Make arrays of all inputs
# ----------------------
infilearray=( $infiles ) 
regionarray=( $regions )
namearray=( $names ) 
colarray=( $colors ) 
tot_readsarray=( $tot_reads ) 
[ "$useValue" == "0" ] && { normarray=( $norm ); tot_readsarray=( $tot_reads ) ; }
mkdir "${outprefix}_LSF_logs"

# Output user-derived options to stdout and to log file
# ----------------------
time_start=$(date)	# time run was started
time_ss=$(date +%s)	# time run was started (in seconds)
echo "Running ends_analysis.sh v1.0 (05/08/2015):"
echo "Run start on: $time_start"
echo "-------------------------"
if [[ "$mode" == "infiles" && "$useValue" == "0" ]]; then
	echo "Calculating average depth over the following files(s):"
	echo "Name	File	Num Reads	Norm Factor	Color"
	for ((i=0;i<${#infilearray[@]};++i)); do
		printf "%s\t%s\t%s\t%s\t%s\n" "${namearray[i]}" "${infilearray[i]}" "${tot_readsarray[i]}" "${normarray[i]}" "${colarray[i]}"
	done
elif [[ "$mode" == "infiles" && "$useValue" != "0" ]]; then
	echo "Calculating average of field $useValue for the following file(s):"
	echo "Name	File	Color"
	for ((i=0;i<${#infilearray[@]};++i)); do
		printf "%s\t%s\t%s\t%s\t%s\n" "${namearray[i]}" "${infilearray[i]}" "${colarray[i]}"
	done
elif [[ "$mode" != "infiles" && "$useValue" == "0" ]]; then
	echo "Calculating average depth for the following file: ${infilearray[0]}"
else
	echo "Calculating average of field $useValue for the following file: ${infilearray[0]}"
fi
echo "-------------------------"
[ "$mode" == "regions" ] && echo "Plotting ends analysis over these regions:" || echo "Plotting ends analysis over regions in: ${regionarray[0]}"
for ((i=0;i<${#regionarray[@]};++i)); do
	[ "$mode" == "regions" ] && echo "${namearray[i]}	${regionarray[i]}"
done
echo "Output file prefix: $outprefix"
echo "-------------------------"
echo "Number of bases outside gene: $numOut"
echo "Number of bases inside gene: $numIn"
echo "Width of intervals: $width"
if [ "$useValue" == "0" ]; then
	[ "$unweighted" = "true" ] && echo "Calculating unweighted average across all genes for each bin" || echo "Calculating weighted average across all genes for each bin"
fi
[ "$barchart" = "true" ] && echo "Plotting results as barchart" || echo "Plotting results as line plot"
echo "-------------------------"
echo ""


# ----------------------
# Step 1: from regions file, make BED file with all intervals
# ----------------------
echo "Step 1: Converting provided regions BED file into intervals of width ${width}..."
if [ "$mode" == "regions" ]; then
	for ((i=0;i<${#regionarray[@]};++i)); do
		echo "Processing ${namearray[i]}..."
		"$path_to_scripts"/ends_analysis_get_intervals.py "${regionarray[i]}" "${outprefix}_${namearray[i]}.bed" "$width" "$numIn" "$numOut" > "${outprefix}_LSF_logs/ends_analysis_get_intervals_${namearray[i]}.txt"
		[ $? != 0 ] && { echo "Error: ends_analysis_get_intervals failed for ${namearray[i]}, see ${outprefix}_LSF_logs/ends_analysis_get_intervals_${namearray[i]}.txt"; exit 1; }
	done
else
	"$path_to_scripts"/ends_analysis_get_intervals.py "${regionarray[0]}" "${outprefix}_regions.bed" "$width" "$numIn" "$numOut"
	[ $? != 0 ] && { echo "Error: ends_analysis_get_intervals failed"; exit 1; }
fi
echo ""

# ----------------------
# Step 2: intersect regions with input files
# ----------------------
echo "Step 2: Intersecting intervals from step 1 with provided input BED files..."

if [ "$mode" == "regions" ]; then
	cols_in_input=$( awk '{print NF}' "${infilearray[0]}" | sort -nu | tail -n 1 )
else
	cols_in_regions=$( awk '{print NF}' "${outprefix}_regions.bed" | sort -nu | tail -n 1 )
fi

pid=()
for ((i=0;i<${#namearray[@]};++i)); do
	[ "$mode" == "regions" ] && cols_in_regions=$( awk '{print NF}' "${outprefix}_${namearray[i]}.bed" | sort -nu | tail -n 1 ) || cols_in_input=$( awk '{print NF}' "${infilearray[i]}" | sort -nu | tail -n 1 )
	if "$parallel"; then
		if ! [ "$useValue" == "0" ]; then
			if [ "$mode" == "regions" ]; then
				cmd="bedtools intersect -a ${infilearray[0]} -b ${outprefix}_${namearray[i]}.bed -wa -wb | cut -f ${useValue},$(( ${cols_in_input}+1 ))-$(( ${cols_in_input} + ${cols_in_regions} )) | awk '{ printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$1 }' > ${outprefix}_${namearray[i]}_intersect.bed"
			else
				cmd="bedtools intersect -a ${infilearray[i]} -b ${outprefix}_regions.bed -wa -wb | cut -f ${useValue},$(( ${cols_in_input}+1 ))-$(( ${cols_in_input} + ${cols_in_regions} )) | awk '{ printf \"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n\", \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$1 }' > ${outprefix}_${namearray[i]}_intersect.bed"
			fi
		else
			if [ "$mode" == "regions" ]; then
				cmd="coverageBed -a ${infilearray[0]} -b ${outprefix}_${namearray[i]}.bed | cut -f 1,2,3,4,5,6,7,8 > ${outprefix}_${namearray[i]}_intersect.bed"
			else
				cmd="coverageBed -a ${infilearray[i]} -b ${outprefix}_regions.bed | cut -f 1,2,3,4,5,6,7,8 > ${outprefix}_${namearray[i]}_intersect.bed"
			fi
		fi
		bsub -o "${outprefix}_LSF_logs/intersect_${namearray[i]}.txt" -K "$cmd" & pid[i]=$!
	else
		echo "Processing ${namearray[i]}..."
		if ! [ "$useValue" == "0" ]; then
			if [ "$mode" == "regions" ]; then
				bedtools intersect -a "${infilearray[0]}" -b "${outprefix}_${namearray[i]}.bed" -wa -wb | cut -f ${useValue},$(( ${cols_in_input}+1 ))-$(( ${cols_in_input} + ${cols_in_regions} )) | awk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $5, $6, $7, $8, $1 }' > "${outprefix}_${namearray[i]}_intersect.bed"
				[ $? != 0 ] && { echo "Error: bedtools intersect failed for ${namearray[i]}"; exit 1; }
			else
				bedtools intersect -a "${infilearray[i]}" -b "${outprefix}_regions.bed" -wa -wb | cut -f ${useValue},$(( ${cols_in_input}+1 ))-$(( ${cols_in_input} + ${cols_in_regions} )) | awk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $5, $6, $7, $8, $1 }' > "${outprefix}_${namearray[i]}_intersect.bed"
				[ $? != 0 ] && { echo "Error: bedtools intersect failed for ${namearray[i]}"; exit 1; }
			fi
		else
			if [ "$mode" == "regions" ]; then
				coverageBed -a "${infilearray[0]}" -b "${outprefix}_${namearray[i]}.bed" | cut -f 1,2,3,4,5,6,7,8 > "${outprefix}_${namearray[i]}_intersect.bed"
				[ $? != 0 ] && { echo "Error: coverageBed failed for ${namearray[i]}"; exit 1; }
			else
				coverageBed -a "${infilearray[i]}" -b "${outprefix}_regions.bed" | cut -f 1,2,3,4,5,6,7,8 > "${outprefix}_${namearray[i]}_intersect.bed"
				[ $? != 0 ] && { echo "Error: coverageBed failed for ${namearray[i]}"; exit 1; }
			fi			
		fi
	fi
done

for ((i=0;i<${#pid[@]};++i)); do
	wait "${pid[i]}" || { echo "Error: step 2 failed for ${namearray[i]}, see ${outprefix}_LSF_logs/intersect_${namearray[i]}.txt"; exit 1; }
done
		
if ! "$Verbose"; then 
	if [ "$mode" == "regions" ]; then
		for ((i=0;i<${#namearray[@]};++i)); do
			rm "${outprefix}_${namearray[i]}.bed"
		done
	else
		rm "${outprefix}_regions.bed"
	fi
fi

echo ""

# ----------------------
# Step 3: process output from step 2 to get counts per bin
# ----------------------
echo "Step 3: Getting average counts per bin..."
for ((i=0;i<${#namearray[@]};++i)); do
#	awk -F$'\t' 'BEGIN{ OFS = FS } { tt[$7]+=1; ss[$7]+=$8 } END { for (id in tt) {print id,ss[id]/tt[id]} } ' "${outprefix}_${namearray[i]}_intersect.bed" | sort -k1n,1 > "${outprefix}_${namearray[i]}_counts.txt"
	if [ "$useValue" != "0" ]; then
		[ "$unweighted" = "true" ] && weight=" --unweighted" || weight=""
	else
		weight=" --unweighted"
	fi
	$path_to_scripts/ends_analysis_process_intersect.py "${outprefix}_${namearray[i]}_intersect.bed" "${outprefix}_${namearray[i]}_counts.txt"${weight}	
done
echo ""

# ----------------------
# Step 3.5: if requested by user, output matrix of values in each bin by feature (one per sample)
# ----------------------
if "$writemat"; then
	echo "Step 3.5: Converting results of step 2 into matrix format..."
	if "$parallel"; then
		pid=()
		for ((i=0;i<${#namearray[@]};++i)); do
	#		echo "Converting results of step 2 into matrix format for ${namearray[i]}..."
			cmd="$path_to_scripts/ends_analysis_make_matrix.py ${outprefix}_${namearray[i]}_intersect.bed ${outprefix}_${namearray[i]}_mat.txt $numIn $numOut $width"
			bsub -o "${outprefix}_LSF_logs/make_matrix_${namearray[i]}.txt" -K "$cmd" & pid[i]=$!
		done
		for ((i=0;i<${#namearray[@]};++i)); do
			wait "${pid[i]}" || { echo "Error: ends_analysis_make_matrix failed for ${namearray[i]}, see ${outprefix}_LSF_logs/make_matrix_${namearray[i]}.txt"; exit 1; }
		done
	else
		for ((i=0;i<${#namearray[@]};++i)); do
			echo "Converting results of step 2 into matrix format for ${namearray[i]}..."
			$path_to_scripts/ends_analysis_make_matrix.py "${outprefix}_${namearray[i]}_intersect.bed" "${outprefix}_${namearray[i]}_mat.txt" $numIn $numOut $width
		done
	fi
	echo ""
fi

# ----------------------
# Step 4: if using counts, or norm factors provided, normalize the average counts according to the provided factors
# ----------------------
if [[ "$useValue" == "0" && mode == "infiles" ]]; then
	echo "Step 4: Multiplying average counts by normalization factors..."
	for ((i=0;i<${#namearray[@]};++i)); do
		normfactor="${normarray[i]}"
		echo "Normalizing average counts for ${namearray[i]}, norm factor = ${normarray[i]}..."
		awk -v n=$normfactor '{printf "%s\t%s\n",$1,$2*n}' "${outprefix}_${namearray[i]}_counts.txt" > "${outprefix}_${namearray[i]}_norm.txt"
		if [ $? != 0 ] ; then echo "Error: normalizing failed for ${outprefix}_${namearray[i]}_counts.txt"; exit 1; fi
		[ "$Verbose" = "false" ] && rm "${outprefix}_${namearray[i]}_counts.txt"
	done
	echo ""
else
	for ((i=0;i<${#namearray[@]};++i)); do
		mv "${outprefix}_${namearray[i]}_counts.txt" "${outprefix}_${namearray[i]}_norm.txt"
	done
	echo "Step 4: Normalization is not performed in -V mode, skipping normalization"
	echo ""
fi

# ----------------------
# Step 5: append all count files and make plots using R
# ----------------------
echo "Step 5: Making plot and outputting results..."
echo "Merging all counts files together..."

sed "s/$/\t${namearray[0]}/" "${outprefix}_${namearray[0]}_norm.txt" > "${outprefix}_average_counts.txt"

[ "$Verbose" = "false" ] && rm "${outprefix}_${namearray[0]}_norm.txt"
for ((i=1;i<${#namearray[@]};++i)); do
	sed "s/$/\t${namearray[i]}/" "${outprefix}_${namearray[i]}_norm.txt" >> "${outprefix}_average_counts.txt"
	[ "$Verbose" = "false" ] && rm "${outprefix}_${namearray[i]}_norm.txt"
done

awk -F$'\t' -v minN="$minN" '{OFS=FS} {if ($3 >= minN) {print $1,$2,$4} else {print $1,"NA",$4}}' "${outprefix}_average_counts.txt" > "${outprefix}_average_counts_filt.txt"

[ "$barchart" = "true" ] && bb=" --makebarchart" || bb=""
[ -z "$colors" ] && cc=" --colors $colors" || cc=""

if [ -z "$colors" ]; then
	$path_to_scripts/ends_analysis_make_plot.R "${outprefix}_average_counts_filt.txt" "$yupper" "$ylabel" "$title" "$outprefix" "$numIn" "$numOut" "$width" "$linewidth"${bb} > /dev/null
else
	$path_to_scripts/ends_analysis_make_plot.R "${outprefix}_average_counts_filt.txt" "$yupper" "$ylabel" "$title" "$outprefix" "$numIn" "$numOut" "$width" "$linewidth" --colors "$colors"${bb} > /dev/null
fi
[ $? != 0 ] && { echo "Error: ends_analysis_make_plot failed"; exit 1; }

echo ""

# ----------------------
# Clean up, remove intermediate files
# ----------------------
if ! "$Verbose"; then 
	for ((i=0;i<${#namearray[@]};++i)); do
		rm "${outprefix}_${namearray[i]}_intersect.bed"
	done
fi
if ! "$Verbose"; then rm -rf "${outprefix}_LSF_logs"; fi

time_end=$(date)	# time run was started
time_es=$(date +%s)	# time run was started
echo "Run ended $time_end"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )"















