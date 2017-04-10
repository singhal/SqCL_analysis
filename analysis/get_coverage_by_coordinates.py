import re
import pandas as pd
import gzip

sp = 'Bothrops_moojeni'
# sp = 'Colobosaura_modesta'
samp = pd.read_csv('/Volumes/sosi/brazil/samples.csv')
samp = samp[samp['sample'] != 'CHUNB45362']
c = '/Volumes/sosi/brazil/metadata/%s_to_probe_coords.csv' % sp
out = '/Volumes/sosi/brazil/coverage/%s_probes.csv' % sp

coords = {}
c = pd.read_csv(c)
for c, s, e in zip(c.locus, c.start, c.stop):
	coords[c] = int((int(s) + int(e)) / 2.0)

inds = samp.ix[samp.lineage == sp, 'sample'].tolist()
inds = sorted(inds)
all = {}

for ind in inds:
	print(ind)
	c = '/Volumes/sosi/brazil/coverage/%s' % ind
	cov = {}
	f = open(c, 'r')
	head = f.next()
	for l in f:
		d = re.split(',', l.rstrip())
		gene, pos = re.split(':', d[0])
		type = re.search('^([^-]+)', gene).group(1)
		if type not in cov:
			cov[type] = {}
		relloc = int(pos) - coords[gene]
		if relloc not in cov[type]:
			cov[type][relloc] = {'val': 0, 'num': 0}
		cov[type][relloc]['num'] += 1
		cov[type][relloc]['val'] += int(d[1])
	f.close()

	for type in cov:
		for loc in cov[type]:
			avg = cov[type][loc]['val'] / float(cov[type][loc]['num'])
			avg = round(avg, 4)
			
			if type not in all:
				all[type] = {}
			if loc not in all[type]:
				all[type][loc] = {}
			if ind not in all[type][loc]:
				all[type][loc][ind] = {}
			all[type][loc][ind] = avg

o = open(out, 'w')
o.write('type,loc,%s\n' % (','.join(inds)))
for type in all:
	for loc in all[type]:
		vals = [all[type][loc][ind] if ind in all[type][loc] else 'NA' for ind in inds]
		vals = [str(x) for x in vals]
		vals = [type, str(loc)] + vals
		o.write('%s\n' % ','.join(vals))
o.close()
