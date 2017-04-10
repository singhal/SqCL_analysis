import re
import pandas as pd
import os
import numpy as np

dir = '/Volumes/sosi/brazil/'
sfile = os.path.join(dir, 'samples.csv')

probes = os.path.join(dir, 'squamate_AHE_UCE_genes_loci.fasta')
p = {}
f = open(probes, 'r')
for l in f:
	if re.search('>', l):
		id = re.search('>([^_]+)_', l).group(1)
		p[id] = ''
	else:
		p[id] += l.rstrip()
f.close()
for x,y in p.items():
	p[x] = len(y)

d = pd.read_csv(sfile)
for sp in d['lineage'].unique():
	prg = os.path.join(dir, 'PRG', '%s.fasta' % sp)
	
	seq = {}
	f = open(prg, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	types = {'AHE': [], 'uce': [], 'gene': []}
	for x,y in seq.items():
		type = re.search('^([^-]+)', x).group(1)
		incr = len(y) / float(p[x])
		types[type].append(incr)

	for type in types:
		print('%s,%s,%s' % (sp,type, np.mean(types[type])))
			
