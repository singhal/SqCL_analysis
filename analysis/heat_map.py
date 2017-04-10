import re
import pandas as pd

a = '/Volumes/sosi/brazil/target_loci/squamate_AHE_UCE_genes_loci.fasta'
d = pd.read_csv('/Volumes/sosi/brazil/samples.csv')
inds = d['lineage'].unique().tolist()
out = '/Volumes/sosi/brazil/metadata/heat_map.csv'

loc = {}
a = open(a, 'r')
for l in a:
	if re.search('>', l):
		id = re.search('>([^_]+)', l).group(1)
		loc[id] = 1
a.close()

seq = {}
for ind in inds:
	seq[ind] = {}
	for l in loc:
		seq[ind][l] = 0

for ind in inds:
	a = '/Volumes/sosi/brazil/PRG/%s.fasta' % ind
	a = open(a, 'r')
	for l in a:
        	if re.search('>', l):
                	id = re.search('>(\S+)', l).group(1)
                	seq[ind][id] = 1
	a.close()

o = open(out, 'w')
locs = sorted(loc.keys())
o.write('lineage,%s\n' % ','.join(locs))
for ind in inds:
	val = [ind] + [seq[ind][loc] for loc in locs]
	o.write('%s\n' % ','.join([str(x) for x in val]))
o.close()
		
	
