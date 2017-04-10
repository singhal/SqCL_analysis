import os
import pandas as pd
import re

dir = '/Volumes/sosi/brazil/'
covdir = os.path.join(dir, 'coverage')

# get samples
sfile = os.path.join(dir, 'samples.csv')
d = pd.read_csv(sfile)
inds = sorted(d['sample'].tolist())

# get locus names
lfile = os.path.join(dir, 'squamate_AHE_UCE_genes_loci.fasta')
locs = {}
f = open(lfile, 'r')
for l in f:
	if re.search('>', l):
		name = re.search('>([^_]+)', l).group(1)
		locs[name] = 1
locs = sorted(locs.keys())
f.close()

# initialize the hash
cov = {}
for ind in inds:
	cov[ind] = {}
	for loc in locs:
		cov[ind][loc] = {'count': 0, 'bases': 0}

# start parsing 
for ind in inds:
	cfile = os.path.join(covdir, ind)
	f = open(cfile, 'r')
	head = f.next()

	for l in f:
		d = re.split(',', l.rstrip())
		loc = re.search('(\S+):', d[0]).group(1)
		c = int(d[1])

		cov[ind][loc]['count'] += 1
		cov[ind][loc]['bases'] += c
	f.close()
	print(ind + " done!")

# print results
out = os.path.join(covdir, 'coverage_by_locus_sample.csv')
o = open(out, 'w')
o.write('locus,' + ','.join(inds) + '\n')
for loc in locs:
	d = []
	d.append(loc)
	for ind in inds:
		if cov[ind][loc]['count'] > 0:
			c = cov[ind][loc]['bases'] / float(cov[ind][loc]['count'])
			c = round(c, 2)
		else:
			c = 0
		d.append(c)
	o.write(','.join([str(x) for x in d]) + '\n')
o.close()
