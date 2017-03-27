# Genome_Biology_2017

Brief description of all scripts used in Picard and Gehring, 2017 Genome Biology.
All scripts written by Colette L Picard (cpicard AT mit DOT edu) and licensed under the Apache License, version 2.0:

------------------------
Copyright 2017 Colette L Picard

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
------------------------

Any questions or issues can be directed to CLP. Some scripts require additional
tools to be installed, which will be indicated where possible.

Note that all scripts noted here can be called without arguments for more detail on options and usage.
Scripts included with Additional file 2 but not described here are only provided as helper
scripts and are not separately described.

------------------------
assign_to_allele.py
------------------------
v.1.3, python script, requires Python 2, tested on 2.7.6
	- required packages sys, os, argparse, re
	
	The purpose of this script is to sort mapped reads into reads from strain1 and reads from strain2
	based on SNPs. Usage:
	
	assign_to_allele.py SNPs.bed reads.sam prefix [options]
	
	To get more detail on options, run this script with no arguments.
	
	This script takes as input a set of mapped reads in SAM format (paired or unpaired),
	a set of SNPs between strain1 and strain2 in BED-like format, and sorts reads into four categories:
		- strain1 = overlapped at least 1 SNP with strain1 allele
		- strain2 = overlapped at least 1 SNP with strain2 allele
		- none = read did not overlap any SNPs
		- conflicted = read overlapped 2+ SNPs and had strain1 allele at at least 1 SNP and strain2 allele at at least 1 SNP
		
	This script has "bisulfite" mode which causes it to ignore C->T and G->A SNPs where appropriate.
	Additionally, paired-end reads that are properly paired (default behavior, see --chgflag option)
	according to the SAM file's SAM flag field will be treated as a single read, so that any SNPs overlapping
	only one member of a pair will apply to the other member. Note that reads in SAM file must be sorted
	by position (this can be done using samtools: samtools sort -n -o outfilename.sam infilename.sam).
	
	Example input files:
	SNP file:
	Chr1	345	346	C>T
	Chr1	501	502	T>C
	Chr1	507	508	T>C
	...	
	where C>T indicates that strain1 allele is C and strain2 allele is T.
	
	SAM file:
	WIGTC-HISEQ2:6:1212:18039:21169#CTTGTA/1;1	99	scaffold_1	1412235	255	89M	=	1412390	247	TGTTTTTTTTTTTATGTATTTTGTGATTATTAAATTTTTTAGGTAAAGGGTTTTTTGATTTTGAGTTGTTTTTGGTTGAGATTTGGTTT	fffffffffffffO__fa_eff_c_bffedfdffffffffOYP^ef]bbfaa]O\NNN\defbNefcefffffNbef[ObWedfdffff	NM:i:24	MD:Z:2C1C1C3C5C3C0C1C6C3C3C4C6C1C0A1C5C3C5C0C3C4C0C5C0	XM:Z:..h.h.h...h.....h...xz.z......h...h...h....h......h.h..z.....z...x.....xz...z....hx.....z	XR:Z:CT	XG:Z:CT	XD:i:1	XA:i:1
	WIGTC-HISEQ2:6:2313:7290:19274#CTTGTA/1;1	99	scaffold_1	1490728	255	91M	=	1490752	116	TTGAGAGAGATGTTTTGTTTTTTAGTATTTTTTTGTGTTTTTTTTTTTAAATTTGTTTAAAGTTTTTAATTTTAGGTTTAATTGGGTTGTT	ffffffffffffffffffffffffffffffffffffffffffffffffdffdff_effObb[effffOO]effOOOeffOOcOZWOefMcf	NM:i:21	MD:Z:1C10C0C1C1C1C2C2C2C13C3C5C0C3C5C2C2C8C3A4C2C0	XM:Z:.z..........hh.z.h.h..x..h..h.............h...h.....xz...h.....h..h..h.................z..z	XR:Z:CT	XG:Z:CT	XD:i:0	XA:i:1	XL:Z:scaffold_1.1490806,scaffold_1.1490810
	WIGTC-HISEQ2:6:2313:7290:19274#CTTGTA/1;1	147	scaffold_1	1490752	255	92M	=	1490728	-116	GTATTTTTTTGTGTTTTTTTTTTTAAATTTGTTTAAAGTTTTTAATTTTAGGTTTAATAGGGTTGTTGTTTATATATTTTTTTAGTGTTATT	cc[fcfffffffbffffffffc\\PPPffbP]a\P\Pefffd]PPffd]]PPeb[PYPYNNObYOYPPdcdP^PPPffffe^PPP\PfePfe	NM:i:20	MD:Z:1C2C13C3C5C0C3C5C2C2C8C8C2C3C1C1G3C0C1C0G9	XM:Z:.h..h.............h...h.....xz...h.....h..h..h.................z..z...h.x.....hh.z..........	XR:Z:GA	XG:Z:CT	XD:i:1	XA:i:0	XL:Z:scaffold_1.1490806,scaffold_1.1490810
	WIGTC-HISEQ2:6:1101:7791:95134#CTTGTA/1;1	99	scaffold_1	1578078	255	92M	=	1578359	365	GATGGAGATAGATGAAGAATAATTTGAATTGATTAGATAGAATTAAATGATTAATAATATTTTGGGGATTTTTTTTTTTTTTATGAAAGGTT	fffffffffffdfffffffffffffffffffffffffffffffffffffffffeffdffffffbfffcfffffffffffM\WbeOObebaOe	NM:i:20	MD:Z:2C5C3C6C3C9C8C0C3C0A1C0C2C2C1C2C8C1C5G1C10	XM:Z:..z.....x...z......h...x.........x........hh...h..hh..h..h.h..z........h.h.......h..........	XR:Z:CT	XG:Z:CT	XD:i:2	XA:i:1
	

------------------------
call_regions.R
------------------------
v.1.2, R script, requires R, tested on 3.3.2
	- requires package optparse
	
	This script was used to determine parent-of-origin regions for each RIL based on
	reads overlapping Col SNPs and Cvi SNPs (see assign_to_allele.py above). Usage:
	
	call_regions.R [options] infile.bed outprefix
	
	To get more detail on options, run this script with no arguments.
	
	This script takes as input a BED-like file:
	Chr1	0	200	0	0
	Chr1	200	400	0	0
	Chr1	400	600	1	14
	Chr1	600	800	1	24
	Chr1	800	1000	0	8

	This file should contain per-window (here window size is 200 bp) counts of strain1 
	reads overlapping the window (field 4) and strain2 reads overlapping the window (field5).
	Note that this script is currently designed to process only a single chromosome at a
	time, and will produce an error if more than one chromosome is included in the input file.
	The script outputs four different sets of regions:
		- regions assigned to strain1 (*_strain1_regions_final_*)
		- regions assigned to strain2 (*_strain2_regions_final_*)		
		- heterozygous regions (*_het_regions_final_*)		
		- unknown regions (*_unk_regions_final_*)		
		
	Parameters can be tweaked to improve performance depending on sequencing experiment
	depth, etc. Parameters used in our analyses are indicated in the methods section.
	
	
------------------------
classify_variation_across_samples.py
------------------------
v.1.0, python script, requires Python 2, tested on 2.7.6
	- required packages sys, os, argparse, re, matplotlib, numpy, scipy
	
	This script was used to classify CGs according to their CG methylation variability
	across 927 ecotypes from Kawakatsu et al. 2016. Usage:
	
	classify_variation_across_samples.py infile outprefix [options]
	
	This script takes as input a file containing frac. methylation data across any number
	of strains (in the example below, four strains), and sorts CGs into different categories
	based on variability across those strains. Example input file:
	
	locus	strain1	strain2	strain3	strain4
	Chr1:188-190	0	0	0	0
	Chr1:214-216	1.0	1.0	1.0	.9
	Chr1:239-241	0	0.9	0.9	0
	
	First column must consist of a single unique location identifier (locus) indicating
	the location of a CG, and then a unique name for each of the strains (e.g. strain1, strain2...).
	Most of the parameters used in this script to define different categories are hardcoded,
	and correspond to Additional file: Figure S9. The exception is the maximum number of strains
	missing data for a CG that are allowed before the CG is censored from the analysis, which can
	be modified with the --maxmissing option (run this script without arguments to see options). 
	Default is no CGs will be censored. A number of parameters used to delineate the different categories
	are listed in the script beginning on line 444, but are not presently accessible through
	command-line options, though they can be modified in the script itself.
	
	This script outputs each CG with its assigned category and subcategory, as well as up to
	5 plots per combination of category and subcategory showing the distributions that have
	been assigned to each category.
	
	
------------------------
ends_analysis.sh
------------------------
v.1.7, bash script, requires Python 2, tested on 2.7.6, and R, tested on 3.3.2
	- required helper scripts (must be in same directory as this script):
		- ends_analysis_get_intervals.py - by Colette L Picard
		- ends_analysis_process_intersect.py - by Colette L Picard
		- ends_analysis_make_plot.R - by Colette L Picard
		- ends_analysis_make_matrix.py - by Colette L Picard (only required if using -M or -C options)
	- required installed on user PATH:
		- bedtools (tested on v2.23.0)
	
	This script was used to generate the average methylation profiles around CGs (plots shown in Fig. 4B).
	Usage:
	
	ends_analysis.sh [options] -i medata1.bed,medata2.bed,...,medataN.bed -r regions.bed -o outprefix
	or
	ends_analysis.sh [options] -i medata.bed -r regions1.bed,regions2.bed,...,regionsN.bed -o outprefix

	Call this script without arguments for more details on options. The two different usage examples are
	for (case 1) plotting multiple methylation profiles over a single set of regions (e.g. genes) and
	(case 2) plotting a single methylation profile over multiple sets of regions (e.g. RIL gain sites and RIL gain background).
	This script was originally designed to plot average coverage around intended intervals, but
	it can also plot average methylation profiles or any other values around intervals of interest.
	
	Example of inputs and usage for methylation data:
	
	Input methylation files (BED format):
	Chr5	6291458	6291459	3	0	0.0
	Chr5	6291459	6291460	3	0	0.0
	Chr5	14680069	14680070	0	3	100.0
	
	For input files, the only crucial thing is that the first 3 fields correspond to location of a
	site with methylation data, in BED-like format (chr,start,end, using 0-based indexing). Fields beyond
	the first 3 can be anything, and you can specify which one you want to use to create the plot. In
	this example, since the % methylation at each position is in the 6th field, you can use the -V 6 option.
	
	Input region files:
	Chr1	5394004	5394006	CG	0	+
	Chr1	29862465	29862467	CG	0	+
	Chr2	5834435	5834437	CG	0	+
	Chr2	7393233	7393235	CG	0	+
	
	These should be "true" BED files with the first 6 required fields. The fourth and fifth field
	are not used and can be anything (but not missing), but the 6th field must have strand information.
	If strand doesn't matter for your analysis, just give all intervals a '+' in that field.
	
	This script allows you to specify how far outside each interval in the region file, and how far
	into each interval you want to make the plot. These are set by -O for # bp outside feature, -I for # bp inside,
	and -w for width of windows to calculate averages over. Example for -O 3, -I 5 and -w 1:
                  
      feature:        5'----------------------3'
      plot range:   ---[-----            -----]---

	In fig. 4B, we only wanted methylation profiles for 100bp upstream and downstream of each feature,
	so we used -O 100 -I 0 -w 20.
	
	Example usage for fig. 4B:
	$path_to_scripts/ends_analysis.sh -r ril_gain_regions.bed,ril_loss_regions.bed,ril_gain_bkd.bed,ril_loss_bkd.bed -i methylation_data.bed -o outprefix -n RILgain,RILloss,hime,lome -O 100 -I 0 -w 20 -V 6 -c "#BE4F51","#56ACC5",gray40,gray88 -t "CGs inherited from $parent" -y "Average % $context methylation" -l 2
	
	
------------------------
get_CG_consistency.sh
------------------------
v.1.0, bash script
	- required helper scripts (must be in same directory as this script):
		- merge_by_column.R (by Colette Picard)
		- fishers_exact.R (by Colette Picard)
	- required installed on user PATH:
		- bedtools (tested on v2.23.0)
	
	Takes as input a set of per-position methylation data in BED-like format:
		Chr1	108	109	0	25	100
		Chr1	109	110	3	23	88.4615
		Chr1	114	115	2	22	91.6667
		Chr1	115	116	0	27	100
	corresponding to fields chr, start, end, # unmethylated reads, # methylated reads and % m
	
	and converts into one record per CG, combining the data from both strands if present:
		Chr1	108	110	3	48	94.1
		Chr1	114	116	2	49	96.1
	
	Requires users to input the genome sequence (FASTA format) that reads were mapped to, so that the
	script can look up each cytosine in the input file and determine if it is on the forward
	strand or the reverse strand.
	
	Usage: get_CG_consistency.sh [options] -i infile.bed -g genome.fa -o outprefix
	
	For all CGs with data on both strands, this script will test whether or not the two
	symmetric sites significantly differ in methylation status, and will output "discordant" sites
	in their own file, without combining them.
	
	Call this script without arguments for a more detailed explanation and list of options.
	

------------------------
predict_logit_train_test.do
------------------------
v.1.0, Stata script, requires Stata version 14

This script runs the logistic regression analysis discussed in the methods section "Logistic regression model fitting",
and requires the input file Additional file 3: full_dataset.txt to be in the same directory
as this script (or, change line 44 to the correct location of this file). Results will be saved
to your current directory, and named "regression_results.dta" and "regression_results.txt"
(a Stata-formated data file and a human-readable tab-del file, respectively). WARNING! this 
script will overwrite both output files if they already exist. After finishing, the script
will appear to exit with an error because of the "exit 1" at the last line - this is just so that
the current data (e.g. the results) will be left in memory instead of deleted.

		
------------------------
predict_logit_train_test_mini.do, predict_logit_train_test_maize.do and predict_logit_train_test_distachyon.do
------------------------
v.1.0, Stata script, requires Stata version 14

These three scripts all run very similar analyses - the process resembles predict_logit_train_test.do
except only seven models, representing all combinations of the three methylation predictors, are tested, and the
dependent var is different (variability across strains instead of dynamic CGs). *_mini.do runs this
analysis on our thaliana data and requires Additional file 3: full_dataset.txt. *_maize.do runs this
analysis on data from maize (Li et al. 2015) and requires Additional file 3: full_maize_B73_data.txt.
*_distachyon.do runs this analysis on B. distachyon data (Eichten et al. 2016) and requires 
Additional file 3: full_distachyon_Bd1_1_data.txt.
	
