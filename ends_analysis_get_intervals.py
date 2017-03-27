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

# --------------------------

# Given an input .bed file containing intervals, returns another .bed file with
# all windows of width w going O bp outside the feature and I bp inside the feature.
# Entries keep their current ID and add a 7th column with the number of the interval.
# Interval numbering is always 5' to 3', uses +/- strand info if present.

# Intervals of length w are marked below along with the start and end pos of the feature.
# Interval numbering is shown below for feature on + strand vs. - strand.
# Assume O = 3 intervals and I = 2.
#								S						E
# 		|-w-|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|
# Feature is +:	5'	  1   2   3   4   5			  6   7   8   9   10	3'
# Feature is -:	3'	  10  9   8   7   6			  5   4   3   2   1		5'

import sys

if len(sys.argv)<=4:
	print "Usage: ends_analysis_get_intervals.py infile.bed outfile.bed width numIn numOut"
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
width = int(sys.argv[3])
numIn = int(sys.argv[4])
numOut = int(sys.argv[5])

print "Number of bins being output: (S = feature start, E = feature end)"
print "5'  "+str(numOut/width)+"---->S<----"+str(numIn/width)+"  Feature  "+str(numIn/width)+"---->E<----"+str(numOut/width)+"  3'"
tot_bins = (2*numOut/width)+(2*numIn/width)
print "Total bins per feature =",str(tot_bins)

f = open(infile, 'r') 
o = open(outfile, 'w')
line = f.readline()
while line:
	r = line.strip().split('\t')
	if len(r) < 4:
		print "Error in ends_analysis_get_intervals: BED file has less than 4 fields. Sample line:"
		print line.strip()
		sys.exit(1)
	elif len(r) == 4:
		r.append("0")
		r.append("+")		# default is + strand
				
	# start with smallest coordinate (feature start - numOut)
	start = int(r[1]) - numOut
	for i in range(1,1+(numIn + numOut)/width):
		leftC = start+(width*(i-1))
		rightC = start+(width*i)
		if r[5] == "+":
			intervalNo = i
		else:
			intervalNo = tot_bins - i + 1
#		print "interval is [",leftC,",",rightC,"], intervalNo =", intervalNo
		# interval must not go below zero, or go more than halfway through feature
		if leftC < int(r[1]) + ((int(r[2])-int(r[1]))/2) and not (leftC < 0 or rightC < 0):
			# print to output file
			o.write(r[0]+'\t'+str(leftC)+'\t'+str(rightC)+'\t'+r[3]+'\t'+r[4]+'\t'+r[5]+'\t'+str(intervalNo)+'\n')

	# now go to end of feature
	start = int(r[2]) - numIn
	for i in range(1,1+(numIn + numOut)/width):
		leftC = start+(width*(i-1))
		rightC = start+(width*i)
		if r[5] == "+":
			intervalNo = (numIn + numOut)/width + i
		else:
			intervalNo = (numIn + numOut)/width - i + 1
#		print "interval is [",leftC,",",rightC,"], intervalNo =", intervalNo
		if leftC >= int(r[1]) + ((int(r[2])-int(r[1]))/2) and not (leftC < 0 or rightC < 0):
			# print to output file
			o.write(r[0]+'\t'+str(leftC)+'\t'+str(rightC)+'\t'+r[3]+'\t'+r[4]+'\t'+r[5]+'\t'+str(intervalNo)+'\n')
			
	line = f.readline()
