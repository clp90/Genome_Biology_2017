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

# Given a BED file with a set of intersected regions from step 2 of ends_analysis.sh, in the format:
# chr	start	end		name			str	bin	count
#Chr1	849239	849289	AT1G03420	0	+	151	0
#Chr1	849289	849339	AT1G03420	0	+	152	0
#Chr1	2069086	2069136	AT1G06740	0	+	8	66.7
#Chr1	2069136	2069186	AT1G06740	0	+	9	60
#Chr1	2069136	2069186	AT1G06740	0	+	9	57.1

# Returns a matrix with each row = a feature (specified by name field above) and each column
# = a bin (specified by "bin" field above). For each unique name:bin, the value for the matrix
# is the average of the final field ("count" in example above). Also adds an empty column between
# the bins flanking the start and end of the feature, as well as 5 empty columns in the middle.

import sys

if len(sys.argv)<=1:
	print "Usage: ends_analysis_make_matrix.py infile.txt outfile.txt numIn numOut width"
	sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]
numIn = int(sys.argv[3])
numOut = int(sys.argv[4])
width = int(sys.argv[5])

# determine where to add empty columns corresponding to start, end and middle
bin_5p = numOut/width
bin_mid = (numIn/width)+(numOut/width)
bin_3p = (numIn/width) + bin_mid
bin_last = numOut/width + bin_3p

# get values for matrix
counts = {}		# store counts per bin in nested dict, in form form {gene: {bin: (N,totcounts)}}

f = open(infile, 'r') 
line = f.readline()
while line:
	r = line.strip().split('\t')
	# r[3] = gene ID, r[6] = bin, r[7] = value
	feature = r[3]
	bin = int(r[6])
	value = float(r[7])
	
	# add value to the dict
	if feature not in counts:
		counts[feature] = {}
		counts[feature][bin] = [1,value]
	else:
		if bin not in counts[feature]:
			counts[feature][bin] = [1,value]
		else:
			counts[feature][bin][0]+=1
			counts[feature][bin][1]+=float(r[7])
	line = f.readline()
	
f.close()

# print results
f = open(outfile, 'w')
f.write("ID")

# write positions of bins along first line
for bin in range(1,bin_5p+1):
	f.write('\t'+"5p:-"+str(numOut-(width*(bin-1)))+"-"+str(numOut-(width*bin)))
f.write('\tfeature_start')
for bin in range(bin_5p+1,bin_mid+1):
	f.write('\t'+"5p:"+str(width*(bin-bin_5p-1))+"-"+str(width*(bin-bin_5p)))
f.write('\tfeature_mid\tfeature_mid\tfeature_mid\tfeature_mid')
for bin in range(bin_mid+1,bin_3p+1):
	f.write('\t'+"3p:-"+str(numIn-(width*(bin-bin_mid-1)))+"-"+str(numIn-(width*(bin-bin_mid))))
f.write('\tfeature_end')
for bin in range(bin_3p+1,bin_last+1):
	f.write('\t'+"3p:"+str(width*(bin-bin_3p-1))+"-"+str(width*(bin-bin_3p)))
f.write('\n')

for gene in counts:
	f.write(gene)
	for bin in range(1,bin_5p+1):
		if not bin in counts[gene]:
			f.write('\t')		# add missing value
		else:
			f.write('\t'+str(counts[gene][bin][1]/float(counts[gene][bin][0])))
	# write in one empty field to mark start position of feature
	f.write('\t')
	for bin in range(bin_5p+1,bin_mid+1):
		if not bin in counts[gene]:
			f.write('\t')		# add missing value
		else:
			f.write('\t'+str(counts[gene][bin][1]/float(counts[gene][bin][0])))
	# write in four empty fields to mark middle of feature
	f.write('\t\t\t\t')
	for bin in range(bin_mid+1,bin_3p+1):
		if not bin in counts[gene]:
			f.write('\t')		# add missing value
		else:
			f.write('\t'+str(counts[gene][bin][1]/float(counts[gene][bin][0])))
	# write in one empty field to mark end position of feature
	f.write('\t')
	for bin in range(bin_3p+1,bin_last+1):
		if not bin in counts[gene]:
			f.write('\t')		# add missing value
		else:
			f.write('\t'+str(counts[gene][bin][1]/float(counts[gene][bin][0])))
	f.write('\n')
	

f.close()















