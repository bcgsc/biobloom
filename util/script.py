import csv

# taxid : [[sequid, .., ...], count]
count = {}

with open('SRR172903.kraken', 'rb') as f:
	reader = csv.reader(f, delimiter='\t')
	for row in reader:
		taxid = int(row[2]);
		seqid = row[1];
		if taxid in count:
			count[taxid][0].append(seqid)
			count[taxid][1] += 1
		else:
			count[taxid] = [[seqid], 1]

keylist = count.keys()
keylist.sort()

with open('counts.txt', 'wb') as g:
	writer = csv.writer(g, delimiter='\t')
	for key in keylist:
		writer.writerow([key, count[key][1]])

with open('sequences.txt', 'wb') as g:
	writer = csv.writer(g, delimiter='\t')
	for key in keylist:
		writer.writerow([key, count[key][0]])






