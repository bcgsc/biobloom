import csv
import sys

# taxid : [count, [seqid, ..., ...]]
count = {}

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	for row in reader:
		taxid = int(row[2]);
		if taxid in count:
			count[taxid][0] += 1
		else:
			count[taxid] = [1, []]

keylist = count.keys()
keylist.sort()

for taxid in keylist:
	with open('seqid2taxid.map', 'r') as d:
		reader = csv.reader(d, delimiter='\t')
		for row in reader:
			if taxid == int(row[1]):
				count[taxid][1].append(row[0])

# print out: taxid	count	[seqid, ..., ...]
with open(outfile, 'w') as g:
	writer = csv.writer(g, delimiter='\t')
	for key in keylist:
		writer.writerow([key, count[key][0], count[key][1]])






