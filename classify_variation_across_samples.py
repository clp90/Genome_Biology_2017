#!/usr/bin/env python

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
   
''' 
-------------------------
v.1.0	03/27/2017
by Colette L Picard

Usage: classify_variation_across_samples.py infile outprefix [options]

See below for details.

Version history:
v.1.0 - initial build (03/27/2017)

-------------------------
'''
 
import sys, os, argparse, re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm
from scipy.stats import chisquare

if len(sys.argv) == 1:
	print "-------------------------"
	print "classify_variation_across_samples v1.0		by Colette L. Picard, 10/18/2016"
	print "-------------------------"
	print """This script attempts to take methylation data across many samples (can be for individual
sites or larger regions) and to classify each data point (row) according to the type
of variability in methylation exhibited across the samples. Originally designed to look
at variation across all the strains of the 1001 methylomes project.

Example input:
locus	strain1	strain2	strain3	strain4
Chr1:188-190	0	0	0	0
Chr1:214-216	1.0	1.0	1.0	.9
Chr1:239-241	0	0.9	0.9	0

Input file must have the first column as some sort of unique identifier of the rows (all
other rows are assumed to represent samples). Input file must also have header with
sample names.

Each row will be classified into one of these categories based on variation across the samples:
- strong unimodal : all samples' values are tightly clustered around a single peak:
			
            x
            x
            x
            xx
           xxxx
        <------------------------>
        
    characterized by: single maximum, low std devation
    modes: low, intermediate, or high methylation; lorange or hirange (hirange = only couple strains differ, lorange = all same)

- weak unimodal : all samples' values are loosely clustered, but there is still one peak

             x
            xxx
            xxxx
           xxxxxxx
         xxxxxxxxxxxxx
        <------------------------>

    characterized by: single maximum, intermediate std devation
    modes: low, intermediate, or high methylation
        
- weak bimodal : all samples' values are loosely clustered, and there are two peaks

            
             xx
            xxxx              x
           xxxxxxx      xxxxxxxxx
         xxxxxxxxxxxxxxxxxxxxxxxxx
        <------------------------>

    characterized by: two maxima, intermediate std devation around both
    modes: more highly methylated samples, more low, or equal
        
- strong bimodal : all samples' values are tightly clustered around one of two peaks

            x
            x                x
            x                x
            xx               xx
           xxxx             xxxxx
        <------------------------>

    characterized by: two maxima, intermediate std devation
    modes: mostly highly methylated (with a few unmethylated), mostly unmethylated, equal

- uniform/random : no well defined, obvious peaks; most samples are different from the others

    
            xxxx x  xx   xx xx
           xxxxxxxxxxxxxxxxxxxxxx
         xxxxxxxxxxxxxxxxxxxxxxxxx
        <------------------------>

	characterized by: no strong maxima
	modes: none
		
- L-shaped : strong peak at one extreme, with all remaining samples uniformly distributed

         x
         xx
         xx                
         xx               
         xxxxx xxxxxxx xxxxx xxxx
         xxxxxxxxxxxxxxxxxxxxxxxxx
        <------------------------>

	characterized by: one strong maxima; once removed, remaining distribution is roughly uniform
	modes: peak at 0% (unmethylated) or peak at 100% (methylated)
		

"""
	print "\nUsage: classify_variation_across_samples.py infile outprefix [options]"
	print "Options:"
	print "--maxmissing : "
	print "-------------------------"
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', help = 'input file, first column is row ID, all other columns are data for analysis')
parser.add_argument('outprefix', help = 'Prefix for output files')
parser.add_argument('--maxmissing', default = None, help = 'max number of strains with missing values for a CG allowed before CG is censored')

args = parser.parse_args()

infile = args.infile
outprefix = args.outprefix
maxmissing = args.maxmissing

print "Running classify_variation_across_samples v1.0		by Colette L. Picard, 10/18/2016"
print "-------------------------"
print "Input file:",infile
print "Prefix for output files:",outprefix
print "-------------------------"

#-------------------------------------------------------------
# Helper functions
def plot_histogram(bins,data,width,ymax,color,title,outfilename):
	plt.bar(bins,data,width=width,color=color)
	plt.ylim(0, ymax)
	plt.title(title)
	plt.xlabel('% methylation')
	plt.ylabel('fraction of samples')
	plt.savefig(outfilename)
	plt.clf()	

def fxdblpeaks(peaks):
	# if peaks contains more than one peak in a row (e.g. because of ties, identical counts in adjacent bins),
	# re-assigns peak to single bin (leftmost bin if pair includes bin 2, rightmost if pair includes bin 4)	
	#	P	P	.	.	.	-> 	P	.	.	.	.
	# 	.	P	P	.	.	-> 	.	P	.	.	.
	#	.	.	P	P	.	-> 	.	.	.	P	.
	#	.	.	.	P	P	-> 	.	.	.	.	P
	
	# For runs of 3 peaks, keep middle value:
	#	P	P	P	.	.	-> 	.	P	.	.	.
	# 	.	P	P	P	.	-> 	.	.	P	.	.
	#	.	.	P	P	P	-> 	.	.	.	P	.
	
	for i in range(0,3):
#		print "Looking at bins",i,"and",i+1,"and",i+2
#		print "Values are:",peaks[i],"and",peaks[i+1],"and",peaks[i+2]	
		if peaks[i] == peaks[i+1] and peaks[i] == peaks[i+2] and peaks[i] == True:
			if i == 0:
				peaks[0] = False; peaks[2] = False
			elif i == 1:
				peaks[1] = False; peaks[3] = False
			elif i == 2:
				peaks[2] = False; peaks[4] = False
				
	for i in range(0,4):
#		print "Looking at bins",i,"and",i+1
#		print "Values are:",peaks[i],"and",peaks[i+1]	
		if peaks[i] == peaks[i+1] and peaks[i] == True:
			if i == 0:
				peaks[1] = False
			elif i == 1:
				peaks[2] = False
			elif i == 2:
				peaks[2] = False
			elif i == 3:
				peaks[3] = False
				
	return peaks
				
def get_bimodal_subcat(p1,p1val,p2,p2val):
	if p1val > foldchg_for_mostly * p2val:
		subcat = 'mostly_'+binsstr[p1]+'_rem_'+binsstr[p2]
	elif p1val > foldchg_for_biased * p2val:
		subcat = 'biased_'+binsstr[p1]+'_rem_'+binsstr[p2]
	else:
		v1 = sorted([p1,p2])[0]; v2 = sorted([p1,p2])[1]
		subcat = 'similar_'+binsstr[v1]+'_'+binsstr[v2]	
	return subcat
	
def get_Lshaped_subcat(p1,p1val):
	if p1val >= min_peak_Lshaped_mostly:
		subcat = 'mostly_'+binsstr[p1]
	elif p1val >= min_peak_Lshaped_biased:
		subcat = 'biased_'+binsstr[p1]
	else:
		subcat = 'similar_'+binsstr[p1]
	
	return subcat
	
	
def classify_unimodal(binneddata,peaks):
	peakloc = peaks.index(True)
	peakval = binneddata[peakloc]
	aroundpeakval = sum(binneddata[max(0,peakloc-1):min(peakloc+2,len(binneddata))])		
	dif = aroundpeakval - peakval
	if peakloc == 1 or peakloc == 2 or peakloc == 3:
		aroundpeakval = peakval + dif/2			# add only half of outer bin value, since when peak is in
												# middle this value sums two peaks, but if peak is on edge only one												
#	print peakloc, peakval, aroundpeakval
	if peakval >= min_binfrac:
		classify = 'unimodal_sharp'
		subcat = binsstr[peakloc]
	elif aroundpeakval >= min_binfrac_inter:
		classify = 'unimodal_inter'
		subcat = binsstr[peakloc]
	else:
		# a peak could be so broad it's essentially uniform
		SSerrall = sum([(x-0.2)**2 for x in binneddata])
		binneddatacopy = binneddata[:].tolist()
		del binneddatacopy[peakloc]; ss = sum(binneddatacopy)
		deltaYall = peakval - min(binneddata)
		deltaYnopeak = max(binneddatacopy) - min(binneddatacopy)
		binneddatacopyadj = [x/ss for x in binneddatacopy]
		SSerrnopeak = sum([(x-0.25)**2 for x in binneddatacopyadj])
		
#		print "SSerr =",SSerrall," deltaY =",deltaYall
#		print "SSerrnopeak =",SSerrnopeak," deltaYnopeak =",deltaYnopeak

		if SSerrall <= max_SSerr_for_uniform and deltaYall <= max_yrange_for_uniform:
			classify = 'uniform'
			subcat = 'none'
		elif SSerrnopeak <= max_SSerr_for_uniform and (peakloc == 0 or peakloc == 4) and deltaYnopeak <= max_yrange_for_uniform and SSerrall/SSerrnopeak >= min_FC_SSerr_for_Lshaped:
			classify = 'L-shaped'
			subcat = get_Lshaped_subcat(peakloc,peakval)
		else:
			classify = 'unimodal_broad'
			subcat = binsstr[peakloc]

	return classify,subcat
	
	
def classify_bimodal(binneddata,peaks):
	# classify bimodal distributions
	# -----------------------
	peaklocs = [i for i, x in enumerate(peaks) if x == True]
	if binneddata[peaklocs[0]] > binneddata[peaklocs[1]]:
		peak1=peaklocs[0]; peak2=peaklocs[1]
	else:	
		peak1=peaklocs[1]; peak2=peaklocs[0]
		
	peak1val = binneddata[peak1]; peak2val = binneddata[peak2]
	if peak2 == 0:
		peak2adjacent = binneddata[1]
	elif peak2 == 4:
		peak2adjacent = binneddata[3]
	else:
		peak2adjacent = max(binneddata[peak2-1],binneddata[peak2+1])
		
	SSerrall = sum([(x-0.2)**2 for x in binneddata])
	binneddatacopy = binneddata[:].tolist()
	del binneddatacopy[peak1]; ss = sum(binneddatacopy)
	deltaYall = peak1val - min(binneddata)
	deltaYnopeak = max(binneddatacopy) - min(binneddatacopy)
	binneddatacopy = [x/ss for x in binneddatacopy]
	SSerrnopeak = sum([(x-0.25)**2 for x in binneddatacopy])
	
#	print "Peak 1:",peak1,peak1val
#	print "Peak 2:",peak2,peak2val
#	print "Peak 1+2:",peak1val+peak2val
#	print "Peak2 - adj:",peak2val-peak2adjacent
#	print "Peak1/2 ratio:",peak1val/peak2val
#	print "SSerr =",SSerrall," deltaY =",deltaYall
#	print "SSerrnopeak =",SSerrnopeak," deltaYnopeak =",deltaYnopeak
#	print "SSerr/SSerrnopeak =",SSerrall/SSerrnopeak

	# if peak2 is small, identify cases where this secondary peak shouldn't really count as a peak
	if peak2val < 0.1 and peak2val-peak2adjacent < censor2_dif_to_adj and peak1val/peak2val >= censor2_min_FC:
#		print "Treating",lineNo,"as unimodal"
		classify,subcat = classify_unimodal(binneddata, peaks)
	elif abs(peak1 - peak2) < 2:
		print "Error: two adjacent bins were called as peaks!"
		sys.exit(1)
	elif abs(peak1 - peak2) == 2:			# peaks are just one bin apart
		# check if peak is 'real' if it's very similar to adjacent peak
		if peak2val - peak2adjacent <= censor2_dif_to_adj:
			binneddatacopy = binneddata[:].tolist()
			del binneddatacopy[peak2]
			binneddatacopy = np.array(binneddatacopy)
			pp2 = np.r_[True, binneddatacopy[1:] >= binneddatacopy[:-1]] & np.r_[binneddatacopy[:-1] >= binneddatacopy[1:], True]
			peaksfilt = [True if p == True and v >= min_peak_height else False for p,v in zip(pp2,binneddatacopy)]

			# deal with rare case where after removing peak 2, there are still two bins with same value (which might be wrongly called as peaks)
			peaklocs = [i for i, x in enumerate(peaksfilt) if x == True]
			if peaksfilt.count(True) >= 2:
				if binneddata[peaklocs[0]] > binneddata[peaklocs[1]]:
					p1=peaklocs[0]; p2=peaklocs[1]
				else:	
					p1=peaklocs[1]; p2=peaklocs[0]
				p2val = binneddata[p2]
				if len([x for x in binneddata if x == p2val]) >= 2:
					peaksfilt[p2] = False

		if peak2val - peak2adjacent <= censor2_dif_to_adj and peaksfilt.count(True) == 1:
			peaks[peak2] = False
			classify,subcat = classify_unimodal(binneddata, peaks)
		elif SSerrall <= max_SSerr_for_uniform and deltaYall <= max_yrange_for_uniform:
			classify = 'uniform'
			subcat = 'none'
		elif SSerrnopeak <= max_SSerr_for_uniform and (peak1 == 0 or peak1 == 4) and deltaYnopeak <= max_yrange_for_uniform and SSerrall/SSerrnopeak >= min_FC_SSerr_for_Lshaped:
			classify = 'L-shaped'
			subcat = get_Lshaped_subcat(peak1,peak1val)
		else:
			classify = 'bimodal_short'
			subcat = get_bimodal_subcat(peak1,peak1val,peak2,peak2val)
	elif peak1val + peak2val >= min_binfrac_bimodal:
		classify = 'bimodal_sharp'
		subcat = get_bimodal_subcat(peak1,peak1val,peak2,peak2val)		
	elif SSerrall <= max_SSerr_for_uniform and deltaYall <= max_yrange_for_uniform:
		classify = 'uniform'
		subcat = 'none'
	elif SSerrnopeak <= max_SSerr_for_uniform and (peak1 == 0 or peak1 == 4) and deltaYnopeak <= max_yrange_for_uniform and SSerrall/SSerrnopeak >= min_FC_SSerr_for_Lshaped:
		classify = 'L-shaped'
		subcat = get_Lshaped_subcat(peak1,peak1val)
	elif peak1val + peak2val >= min_binfrac_bimodal_inter:
		classify = 'bimodal_inter'
		subcat = get_bimodal_subcat(peak1,peak1val,peak2,peak2val)		
	else:
		classify = 'bimodal_broad'
		subcat = get_bimodal_subcat(peak1,peak1val,peak2,peak2val)		

	return classify,subcat
	

def classify_trimodal(binneddata,peaks):
	peaklocs = [i for i, x in enumerate(peaks) if x == True]
	if binneddata[peaklocs[0]] > binneddata[peaklocs[1]] and binneddata[peaklocs[0]] > binneddata[peaklocs[2]]:
		peak1=peaklocs[0]
		if binneddata[peaklocs[1]] > binneddata[peaklocs[2]]:
			peak2=peaklocs[1]; peak3 = peaklocs[2]
		else:
			peak2=peaklocs[2]; peak3 = peaklocs[1]
	elif binneddata[peaklocs[1]] > binneddata[peaklocs[0]] and binneddata[peaklocs[1]] > binneddata[peaklocs[2]]:
		peak1=peaklocs[1]
		if binneddata[peaklocs[0]] > binneddata[peaklocs[2]]:
			peak2=peaklocs[0]; peak3 = peaklocs[2]
		else:
			peak2=peaklocs[2]; peak3 = peaklocs[0]					
	else:	
		peak1=peaklocs[2]
		if binneddata[peaklocs[0]] > binneddata[peaklocs[1]]:
			peak2=peaklocs[0]; peak3 = peaklocs[1]
		else:
			peak2=peaklocs[1]; peak3 = peaklocs[0]					
		
	peak1val = binneddata[peak1]; peak2val = binneddata[peak2]; peak3val = binneddata[peak3]
	if peak3 == 0:
		peak3adjacent = binneddata[1]
	elif peak3 == 4:
		peak3adjacent = binneddata[3]
	else:
		peak3adjacent = max(binneddata[peak3-1],binneddata[peak3+1])

	SSerrall = sum([(x-0.2)**2 for x in binneddata])
	binneddatacopy = binneddata[:].tolist()
	del binneddatacopy[peak1]; ss = sum(binneddatacopy)
	deltaYall = peak1val - min(binneddata)
	deltaYnopeak = max(binneddatacopy) - min(binneddatacopy)
	binneddatacopy = [x/ss for x in binneddatacopy]
	SSerrnopeak = sum([(x-0.25)**2 for x in binneddatacopy])
		
#	print "Peak 1:",peak1,peak1val
#	print "Peak 2:",peak2,peak2val
#	print "Peak 3:",peak3,peak3val
#	print "Peak 1+2:",peak1val+peak2val
#	print "Peak3 - adj:",peak3val-peak3adjacent
#	print "Peak1/2 ratio:",peak1val/peak2val
#	print "Peak1/3 ratio:",peak1val/peak3val
#	print "SSerr =",SSerrall," deltaY =",deltaYall
#	print "SSerrnopeak =",SSerrnopeak," deltaYnopeak =",deltaYnopeak

	if peak3val-peak3adjacent < censor2_dif_to_adj_peak3:
		# censor peak 3 and treat as bimodal
		peaksfilt[peak3] = False
		classify,subcat = classify_bimodal(binneddata, peaksfilt)
	elif SSerrall <= max_SSerr_for_uniform and deltaYall <= max_yrange_for_uniform:
		classify = 'uniform'
		subcat = 'none'
	elif SSerrnopeak <= max_SSerr_for_uniform and (peak1 == 0 or peak1 == 4) and deltaYnopeak <= max_yrange_for_uniform and SSerrall/SSerrnopeak >= min_FC_SSerr_for_Lshaped:
		classify = 'L-shaped'
		subcat = get_Lshaped_subcat(peak1,peak1val)
	else:
		classify = 'trimodal'
		subcat = 'biased_'+binsstr[peak1]

	return classify,subcat


#-------------------------------------------------------------
try:
	f = open(infile, 'r') 
except IOError, e:
	print e
	print 'Could not open input file',infile,'...'
	sys.exit(2)

try:
	outf = open(outprefix+"_assignments.txt", 'w') 
except IOError, e:
	print e
	print 'Could not create output file',outprefix+'_assignments.txt...'
	sys.exit(2)

header = f.readline()			# this is the header
numsamples = len(header.strip().split('\t')) - 1

print "Header contains",numsamples,"samples"
if maxmissing is None:
	minobs = 0
	print "All rows will be processed regardless of missing values"
else:
	maxmissing = int(maxmissing)
	minobs = numsamples - maxmissing
	if minobs < 0:
		print "Error: max number of allowable missings ("+str(maxmissing)+") is greater than number of samples ("+str(numsamples)+")"
		sys.exit(1)
	print "Only rows with fewer than",maxmissing," missing samples ("+str(minobs)+" non-missing) will be processed"	
print "-------------------------"

# params to add:
global min_peak_height; global min_binfrac; global min_binfrac_inter; global min_dif_to_adj; global max_yrange_for_uniform; global max_SSerr_for_uniform
global min_FC_SSerr_for_Lshaped
min_peak_height = 0.05					# peaks below this height (in fraction) are censored
min_binfrac = 0.90						# min peak height
min_binfrac_inter = 0.8					# min peak height for intermediate unimodal
max_yrange_for_uniform = 0.15			# max difference between highest and lowest peak for 'uniform' distribution
max_SSerr_for_uniform = 0.05			# max Std squared error from uniform dist. to consider approximately uniform
min_FC_SSerr_for_Lshaped = 3			# if L-shaped, SSerrall / uncorrected (peak1 removed) SSerr must be at least this high

# for bimodal distributions, if secondary peak doesn't meet these criteria and has less than 5% of obs, it is censored
# and distribution is treated as unimodal
global censor2_dif_to_adj_peak3; global censor2_min_FC; global censor2_dif_to_adj_peak3
censor2_dif_to_adj = 0.01				# secondary peak must be less than this much greater than lowest bin
censor2_min_FC = 5						# primary peak must be this many times greater than secondary peak
censor2_dif_to_adj_peak3 = 0.05			# for trimodal dists, tertiary peak must be at least this much greater than lowest bin to count

global min_binfrac_bimodal; global min_binfrac_bimodal_inter
min_binfrac_bimodal = 0.85
min_binfrac_bimodal_inter = 0.7

global foldchg_for_mostly; global foldchg_for_biased
foldchg_for_mostly = 4
foldchg_for_biased = 1.5

global min_peak2_Lshaped; global min_peak_Lshaped_mostly; global min_peak_Lshaped_biased
min_peak2_Lshaped = 0.1
min_peak_Lshaped_mostly = 0.8
min_peak_Lshaped_biased = 0.5

maxplots = 5
	
counts = {}						# count occurrences of each classification
processed = 0

binsstr = {0: 'lo', 1: 'medlo', 2: 'med', 3: 'medhi', 4: 'hi'}		# what the 5 bins correspond to (in terms of methylation)

line = f.readline()
processed = 0; skipped = 0; lineNo = 0
while line:
	lineNo+=1
	
	if lineNo % 10000 == 1:
		print "Processing line",lineNo,"..."
	
	ll = line.rstrip('\n').split('\t')
	if len(ll) != numsamples + 1:
		print "Error: line",lineNo,"has only",len(ll),"fields, but should have",(numsamples+1),"according to header"
		sys.exit(1)
	
	data=ll[1:]
	mis = data.count('')

	if (numsamples - mis) >= minobs:
		processed+=1
		
#	if lineNo == 30721:
		
#		print "ID="+ll[0]
#		print "Processing sample #"+str(processed)

		dd = [float(x) for x in data if x != '']
			
		# bin data into histogram with 5 bins
		hist, edges = np.histogram(dd, bins=5, range=(0,100))
		tot = float(sum(hist)); binfrac = [x/tot for x in hist]
		binneddata = np.array(binfrac)			
	
		# get local maxima
		peaks = np.r_[True, binneddata[1:] >= binneddata[:-1]] & np.r_[binneddata[:-1] >= binneddata[1:], True]
		peaksfilt = [True if p == True and v >= min_peak_height else False for p,v in zip(peaks,binneddata)]
		
		# deal with peaks made up of more than one bin due to identical values
		# all of these categories should be really rare (runs of 5 or 4 peaks)
		if peaksfilt == [True, True, True, True, True]:
#			print "Wow, all 5 bins are identical!  Classifying as uniform"
			classify = 'uniform'
			subcat = 'none'
		elif peaksfilt == [True, True, True, True, False]:
			print "4 bins are peaks! Rightskewed"
			if binneddata[4] < 0.5*np.mean(binneddata[0:3]):		# if nonpeak is less than half the avg value of the others
				classify = 'unimodal_broad'
				subcat = 'lo'
			else:													# nonpeak not too dif. from others, class as uniform
				classify = 'uniform'
				subcat = 'none'
		elif peaksfilt == [False, True, True, True, True]:
			print "4 bins are peaks! Leftskewed"
			if binneddata[0] < 0.5*np.mean(binneddata[1:4]):		# if nonpeak is less than half the avg value of the others
				classify = 'unimodal_broad'
				subcat = 'hi'
			else:													# nonpeak not too dif. from others, class as uniform
				classify = 'uniform'
				subcat = 'none'
		else:
			# fix runs of peaks due to identical values by picking one (see fxdblpeaks)
			peaksfilt = fxdblpeaks(peaksfilt)
			
			# classify according to number of peaks identified
			if peaksfilt.count(True) == 1:
				classify,subcat = classify_unimodal(binneddata, peaksfilt)
								
			elif peaksfilt.count(True) == 2:
				classify,subcat = classify_bimodal(binneddata, peaksfilt)
				
			else:
				classify,subcat = classify_trimodal(binneddata, peaksfilt)

		# save result
		print "Classified as:",classify,"	subcat:",subcat			
		outf.write(ll[0]+'\t'+classify+'\t'+subcat+'\n')
		
		# add this result to counts; output plot if fewer than maxplots generated
		if classify not in counts:
			counts[classify] = {}
		if subcat not in counts[classify]:
			counts[classify][subcat] = 1
		else:
			counts[classify][subcat]+=1
			
		if counts[classify][subcat] <= maxplots:
			# bin data into histogram with 50 bins and plot to show original distribution
			hist, newedges = np.histogram(dd, bins=50, range=(0,100))
			newtot = float(sum(hist)); newbinfrac = [x/newtot for x in hist]
			plot_histogram(newedges[:-1], newbinfrac, 2, max(newbinfrac)+0.01, '0.75', ll[0]+', line '+str(lineNo)+', class:'+classify+' '+subcat, outprefix+"_"+classify+"_"+subcat+"_"+str(counts[classify][subcat])+".png")

			# now make 5-bin plot, highlight peaks
			ccs = ['red' if x == True else '0.75' for x in peaksfilt]			
			plot_histogram(edges[:-1], binfrac, 20, max(binfrac)+0.01, ccs, ll[0]+', line '+str(lineNo)+', class:'+classify+' '+subcat, outprefix+"_"+classify+"_"+subcat+"_"+str(counts[classify][subcat])+"_peaks.png")
			plot_histogram(edges[:-1], binfrac, 20, max(binfrac)+0.01, '0.75', ll[0]+', line '+str(lineNo)+', class:'+classify+' '+subcat, outprefix+"_"+classify+"_"+subcat+"_"+str(counts[classify][subcat])+"_peaks_grey.png")
			# make version with fixed ymax/min
			plot_histogram(edges[:-1], binfrac, 20, 1, '0.75', ll[0]+', line '+str(lineNo)+', class:'+classify+' '+subcat, outprefix+"_"+classify+"_"+subcat+"_"+str(counts[classify][subcat])+"_peaks_grey_ymax.png")

	else:
		skipped+=1
		if lineNo > 30721:
			break
	
	line = f.readline()

# output final counts
print ""
print "Done. Summary of results:"
print "-------------------"
print "class	subtype	count"
for cc in counts:
	for ss in counts[cc]:
		print cc+'\t'+ss+'\t'+str(counts[cc][ss])
	






	
	
	
	
