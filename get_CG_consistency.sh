#!/usr/bin/env bash

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
# v1.0 by Colette Picard
# 10/12/2016
# ------------------------------------------------------------------------------------
	
# Usage:
# get_CG_consistency.sh [options] -i infile.bed -g genome.fa -o outprefix

# -------------------------
# Version history:
# v.1.0: initial build											10/12/2016
# -------------------------

# See below for description, or run this script without any arguments.

# -------------------------
# Parameters and required programs/scripts/files:
# -------------------------
read -d '' usage <<"EOF"
v1.0 by Colette Picard, 10/12/2016
This script takes a BED file containing CG methylation data in the format
chr, 0-based start, 0-based end, # unmethylated reads, # methylated, % or fraction methylated
(no header); example:

Chr1	108	109	0	25	100
Chr1	109	110	3	23	88.4615
Chr1	114	115	2	22	91.6667
Chr1	115	116	0	27	100

and for CG sites where there is data on both strands, looks for cases where the two strands
are significantly different (this should be very rare; if it is not it could be an interesting
biological phenomenon, or could just indicate that there are issues with the data). Sites
are flagged according to if they had data on only the forward, only the reverse, on both
(and they were similar) and on both (and they disagreed). Total occurrences of each of the
four categories are counted. Finally, for sites with data on both strands that were not significantly
different, data are consolidated into a single record containing the sum of the number of
unmethylated and methylated reads, and the weighted mean methylation. The start/end coordinates of all
records are adjusted to include the Cs on both strands. Significance of difference between the two strands at each site is evaluated using Fisher's 
exact test. Example output based on input above:

Chr1	108	110	3	48	94.1
Chr1	114	116	2	49	96.1


Usage:
get_CG_consistency.sh [options] -i infile.bed -g genome.fa -o outprefix

User-specified options (defaults indicated in brackets):
Required arguments:
	-i infile : methylation data from CG sites in BED format
	-g genome : genome sequence in FASTA format
	-o outprefix : prefix for all output files
Additional options:
	-s path_to_scripts : path to folder containing all required scripts (note - $scriptDir in default is the location of this script) ["$scriptDir"]
	-d mindif : minimum difference in methylation for sites on opposite strands to be considered significantly different [0.4]
	-p pval : p-value cutoff for sites on opposite strands to be considered significantly different [0.05]
Flag options:
	-F : if script detects a site that is not a C in the reference, ignore it instead of exiting [force=false]
	
Must be installed on your PATH:
	- bedtools (tested on v2.23.0)
Must be in path_to_scripts - set path to this folder using -s option:
	- merge_by_column.R (by Colette Picard)
	- fishers_exact.R (by Colette Picard)
	
------------------------------------------------------------------------------------
EOF
	
[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Required arguments:
# ----------------------
infile=""							# methylation data from CG sites in BED format
genome=""							# genome sequence in FASTA format
outprefix=""							# prefix for all output files

# Additional options:
# ----------------------
path_to_scripts="$scriptDir"							# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
mindif=0.4								# minimum difference in methylation for sites on opposite strands to be considered significantly different
pval=0.05								# p-value cutoff for sites on opposite strands to be considered significantly different

# Flag options:
# ----------------------
force=false							# if script detects a site that is not a C in the reference, ignore it instead of exiting


# ----------------------
while getopts "i:g:o:s:d:p:Fh" opt; do
	case $opt in
		i)	# methylation data from CG sites in BED format
			infile="$OPTARG"
			;;
		g)	# genome sequence in FASTA format
			genome="$OPTARG"
			;;
		o)	# prefix for all output files
			outprefix="$OPTARG"
			;;
		s)	# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
			path_to_scripts="$OPTARG"
			;;
		d)	# minimum difference in methylation for sites on opposite strands to be considered significantly different
			mindif="$OPTARG"
			;;
		p)	# p-value cutoff for sites on opposite strands to be considered significantly different
			pval="$OPTARG"
			;;
		F)	# if script detects a site that is not a C in the reference, ignore it instead of exiting
			force=true
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

# Check that all programs required on PATH are installed
# ----------------------
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required on PATH but was not found"; exit 1; }

# Check that required files are in path_to_scripts (set path to this location with option -s)
# ----------------------
[ ! -f "$path_to_scripts/merge_by_column.R" ] && { echo "Error: could not find required file merge_by_column.R in provided folder (${path_to_scripts})"; exit 1; }
[ ! -f "$path_to_scripts/fishers_exact.R" ] && { echo "Error: could not find required file fishers_exact.R in provided folder (${path_to_scripts})"; exit 1; }

# Check all required inputs are provided
# ----------------------
[ -z "$infile" ] && { echo "Error: -i infile is a required argument (methylation data from CG sites in BED format)"; exit 1; }
[ -z "$genome" ] && { echo "Error: -g genome is a required argument (genome sequence in FASTA format)"; exit 1; }
[ -z "$outprefix" ] && { echo "Error: -o outprefix is a required argument (prefix for all output files)"; exit 1; }

# Check input files exist
# ----------------------
[ -s "$infile" ] || { echo "Error: provided input file $infile is empty"; exit 1; }
[ -s "$genome" ] || { echo "Error: provided genome file $genome is empty"; exit 1; }

# Functions
# ----------------------
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

# Output user-derived options
# ----------------------
time_start=$(date)	# time run was started
time_ss=$(date +%s)	# time run was started (in seconds)
echo "Running get_CG_consistency v1.0 by Colette Picard (10/12/2016)"
echo "Run start on: $time_start"
echo "-------------------------"
echo "Input methylation file: $infile"
echo "Genome FASTA file: $genome"
echo "Output file prefix: $outprefix"
echo "-------------------------"
echo "Min difference in methylation for significance: $mindif"
echo "Max pval for significance: $pval"
echo "-------------------------"
echo ""


# For each cytosine in data, determine if a C (forward strand) or G (reverse strand) in genome
printf "Getting strand of each cytosine in data... "
awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$1":"$2"-"$3}' "$infile" > "${outprefix}_tmp.bed"
bedtools getfasta -fi "$genome" -bed "${outprefix}_tmp.bed" -name -fo "${outprefix}_tmp.fa"
sed '$!N;s/\n/\t/' "${outprefix}_tmp.fa" | sed 's/>//' | sed 's/:/\t/' | sed 's/-/\t/' | sed 's/\tC/\t1/' | sed 's/\tG/\t-1/' | sed 's/\tA/\tERROR/' | sed 's/\tT/\tERROR/' > "${outprefix}_tmp.bed"
	
err=$( grep -F "ERROR" "${outprefix}_tmp.bed" | head -2 )
numErr=$( grep -F "ERROR" "${outprefix}_tmp.bed" | wc -l )
	
if [ "$numErr" -gt 0 ]; then
	if [ "$force" = "false" ]; then
		echo ""
		echo "Error: $numErr records in the input file contain positions that do not correspond to either a C or G in the provided genome"
		echo "Example:"
		echo "$err"
		echo "Are you sure the correct input files were provided?"
		echo "This can also be caused by SNPs between your data and the real genome. To ignore this error, use the -F option."
		exit 1
	else
		echo ""
		echo ""
		echo "WARNING --------------"
		echo "$numErr records in the input file containing positions that do not correspond to either a C or G in provided genome were saved to ${outprefix}_notC.bed and omitted from rest of analysis"
		echo "WARNING --------------"
		echo ""
		awk '!(/ERROR/)' "${outprefix}_tmp.bed" > "${outprefix}_tmp_fx.bed"
		awk '(/ERROR/)' "${outprefix}_tmp.bed" > "${outprefix}_notC.bed"
		rm "${outprefix}_tmp.bed"; mv "${outprefix}_tmp_fx.bed" "${outprefix}_tmp.bed"
	fi
fi

awk '$4 == "1" {print $0}' "${outprefix}_tmp.bed" > "${outprefix}_tmp_f.bed"
awk '$4 == "-1" {print $0}' "${outprefix}_tmp.bed" > "${outprefix}_tmp_r.bed"
rm "${outprefix}_tmp.bed" "${outprefix}_tmp.fa"
echo "DONE"


# Merge back to original data
printf "Merging back strand information to the original data..."
echo "chr	start	end	strand	unme	me	revstart	revend" > "${outprefix}_f.txt"
echo "chr	revstart	revend	revstrand	revunme	revme" > "${outprefix}_r.txt"
bedtools intersect -a "${outprefix}_tmp_f.bed" -b "$infile" -wa -wb -f 1 -r | cut -f1-4,8,9 | awk -F$'\t' '{OFS=FS} {print $0,$2+1,$2+2}' >> "${outprefix}_f.txt"
bedtools intersect -a "${outprefix}_tmp_r.bed" -b "$infile" -wa -wb -f 1 -r | cut -f1-4,8,9 >> "${outprefix}_r.txt"
rm "${outprefix}_tmp_f.bed" "${outprefix}_tmp_r.bed"
echo " DONE"
	
	
# Combine CGs where both the C and G (C on reverse) are in the data
printf "Combining Cs on the forward and reverse strand..."
$path_to_scripts/merge_by_column.R "${outprefix}_f.txt" "${outprefix}_r.txt" chr,revstart,revend "${outprefix}_merged.txt" > "${outprefix}_mergelog.txt"
[ $? != 0 ] && { echo "error when running merge_by_column, see ${outprefix}_mergelog.txt"; exit 1; }
rm "${outprefix}_mergelog.txt"

awk -F$'\t' '{OFS=FS} $9 == "" {print $1,$4,$5+1,$6,$7,$8}' "${outprefix}_merged.txt" > "${outprefix}_fonly.txt"
awk -F$'\t' '{OFS=FS} $4 == "" {print $1,$2-1,$3,$9,$10,$11}' "${outprefix}_merged.txt" > "${outprefix}_ronly.txt"		
awk -F$'\t' '{OFS=FS} $4 != "" && $9 != "" {print $1,$4,$3,"both",$7,$8,$10,$11}' "${outprefix}_merged.txt" > "${outprefix}_both.txt"

rm "${outprefix}_merged.txt" "${outprefix}_f.txt" "${outprefix}_r.txt"
echo " DONE"


# Run Fisher's exact test to identify sites where the two strands are significantly different
printf "Running Fisher's exact test to identify sites where the two strands disagree..."
$path_to_scripts/fishers_exact.R --numbers 5,6,7,8 --infile "${outprefix}_both.txt" --outfile "${outprefix}_both_fishers.txt" --header --name CG_consistency > "${outprefix}_fisherlog.txt"
[ $? != 0 ] && { echo "error when running fishers_exact, see ${outprefix}_fisherlog.txt"; exit 1; }
rm "${outprefix}_fisherlog.txt"
echo " DONE"

# Identify sites that fail the test and flag; add back the single-strand sites to get one output file
printf "Identifying sites that disagree between strands and summarizing..."
awk -F$'\t' -v p="$pval" -v d="$mindif" '{OFS=FS} {
	if (NR==1) {
		print "chr","start","end","unme","me","p_me","CG_status"
	} else if($10 < p && (($6/($5+$6)) - ($8/($7+$8)) > d || ($8/($7+$8)) - ($6/($5+$6)) > d)) {
		print $1,$2,$2+1,$5,$6,($6/($6+$5))*100,"-2"
		print $1,$3-1,$3,$7,$8,($8/($8+$7))*100,"-2"
	} else {
		print $1,$2,$3,$5+$7,$6+$8,(($6+$8)/($5+$7+$6+$8))*100,"2"
	}}' "${outprefix}_both_fishers.txt" | tail -n +2 > "${outprefix}.bed"

consistent=$( awk -F$'\t' 'BEGIN{num=0} $7=="2" {num+=1} END {print num}' "${outprefix}.bed" )
inconsistent=$(( $( awk -F$'\t' 'BEGIN{num=0} $7=="-2" {num+=1} END {print num}' "${outprefix}.bed" ) / 2 ))		# records output separately per strand so 2 records per CG

# filter out inconsistent sites and save to separate file
if [ "$inconsistent" -gt 0 ]; then
	awk -F$'\t' '{OFS=FS} $7==-2 {print $0}' "${outprefix}.bed" > "${outprefix}_inconsistent.bed"
	awk -F$'\t' '{OFS=FS} $7!=-2 {print $0}' "${outprefix}.bed" > "${outprefix}_pass.bed"
else
	cat "${outprefix}.bed" > "${outprefix}_pass.bed"
fi

rm "${outprefix}.bed"

awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$5,$6,($6/($5+$6))*100,$4}' "${outprefix}_fonly.txt" >> "${outprefix}_pass.bed"
awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$5,$6,($6/($5+$6))*100,$4}' "${outprefix}_ronly.txt" >> "${outprefix}_pass.bed"

fonly=$( wc -l "${outprefix}_fonly.txt" | awk '{print $1}' )
ronly=$( wc -l "${outprefix}_ronly.txt" | awk '{print $1}' )

rm "${outprefix}_fonly.txt" "${outprefix}_ronly.txt" "${outprefix}_both_fishers.txt" "${outprefix}_both.txt"
echo " DONE"

echo ""
echo "Summary of results:"
echo "Sites with data only on forward strand: $fonly"
echo "Sites with data only on reverse strand: $ronly"
echo "Sites with data on both strands that are -not- significantly different: $consistent"
echo "Sites with data on both strands that -are- significantly different: $inconsistent"
echo ""

pass=$(( $fonly + $ronly + $consistent ))

echo "All $pass sites -excluding- CGs significantly different on both strands saved to ${outprefix}_pass.bed"
[ "$inconsistent" -gt 0 ] && echo "$inconsistent sites where CGs significantly different on both strands saved to ${outprefix}_inconsistent.bed"

echo "7th field indicates result (1 = forward only, -1 = reverse only, 2 = both+agree, -2 = both+disagree)"
	
time_end=$(date)	# time run was started
time_es=$(date +%s)	# time run was started
echo ""
echo "Run ended $time_end"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )"
	
	
