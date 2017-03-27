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

# v.1.0 - 06/01/2015
# v.1.1 - 08/12/2015
# Changes the way averages within bins are computed. Previously, averages were
# simply taken over all entries corresponding to that bin, which effectively
# leant more weight to genes/etc. (identified by "name" field) with more observations
# within that bin. Now, calculates average over all sites/reads within a bin per
# gene first, then takes average over all genes (so all genes weighted equally).
# v.1.2 - 04/15/2016
# added third column in output that will contain # of genes/regions that the bin
# average was calculated over

# Given a BED file with a set of intersected regions from step 2 of ends_analysis.sh, in the format:
# chr	start	end		name			str	bin	count
# Chr1	32750	32775	Chr1_32001	0	+	23	2
# Chr1	32750	32775	Chr1_48001	0	+	23	10
# Chr1	48150	48175	Chr1_48001	0	+	39	2
# Chr1	49150	49175	Chr1_49001	0	+	15	2

# Returns average of the last field (either COUNT or the value of VALUE field in BED file),
# over bin, where each gene (identified by different name) has equal weight.

import sys, os, argparse

if len(sys.argv)<=1:
	print "Usage: ends_analysis_process_intersect.py [options] infile.bed outfile.txt"
	sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument('infile', help = 'input file')
parser.add_argument('outfile', help = 'output file')
parser.add_argument('--unweighted', default = False, action="store_true",
	help = 'Calculate average of bin over all genes as unweighted average, so all genes weighted equally (default weighted by gene coverage in bin)')
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
unweighted = args.unweighted

counts = {}		# store counts per bin in dict
bins = []		# store a list of all bins in file

f = open(infile, 'r') 
line = f.readline()
while line:
	r = line.strip().split('\t')
	currentbin = int(r[6]); geneID = r[3]; value = float(r[7])
	if currentbin not in counts:
		counts[currentbin] = {}
		counts[currentbin][geneID] = [1,value]		# in counts, store 2-list with left value = # sites in bin, right value = sum of value
		bins.append(currentbin)
	else:
		if geneID not in counts[currentbin]:
			counts[currentbin][geneID] = [1,value]
		else:
			counts[currentbin][geneID][0]+=1
			counts[currentbin][geneID][1]+=value
	line = f.readline()
f.close()

# print results
f = open(outfile, 'w')
bins.sort()

for bin in bins:
	# calculate average over all genes according to requested method (weighted or unweighted)
	totvalue = 0; N = 0; totGenes = 0
#	print counts[bin]
	for gene in counts[bin]:
		totGenes += 1
		if unweighted == True:
			totvalue += counts[bin][gene][1]/float(counts[bin][gene][0])
			N += 1
		else:
			totvalue += counts[bin][gene][1]
			N += float(counts[bin][gene][0])

#	print "Printing:",str(bin)+'\t'+str(totvalue/N)+'\t'+str(totGenes)
	f.write(str(bin)+'\t'+str(totvalue/N)+'\t'+str(totGenes)+'\n')

f.close()


'''
# Alternative method, puts more weight on genes with more observations instead of calculating per-gene
# average and showing that.

f = open(infile, 'r') 
line = f.readline()
while line:
	r = line.strip().split('\t')
	if int(r[6]) not in counts:
		counts[int(r[6])] = [0,0]
	counts[int(r[6])][0]+=1
	counts[int(r[6])][1]+=float(r[7])
	line = f.readline()
	
f.close()

# print results
f = open(outfile, 'w')

for bin in counts.keys():
	# calculate average over all genes
#	print "For bin",bin,"total counts are:",counts[bin][0],"and sum is",counts[bin][1],"and avg frac is",float(counts[bin][1])/counts[bin][0]
	f.write(str(bin)+'\t'+str(float(counts[bin][1])/counts[bin][0])+'\n')

f.close()

'''