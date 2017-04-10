import re
import os
import glob
import subprocess
import numpy as np
from collections import Counter

alndir = '/Volumes/sosi-1/brazil/phylogeny/alignments/'
coords = '/Volumes/sosi-1/brazil/metadata/alignment_to_probe_coords.csv'
out1 = '/Volumes/sosi-1/brazil/metadata/amount_missing.csv'
miss_cap = 0.33

coor = {}
f = open(coords, 'r')
# head = f.next()
head = f.__next__()
for l in f:
	d = re.split(',', l.rstrip())
	mid = int((int(d[1]) + int(d[2])) / 2.0) 
	coor[d[0]] = mid
f.close()

fout = open(out1, 'w')
fout.write('locus,abs_pos,rel_pos,per_miss\n')

def parse_col(col):
	bases = Counter(col)
	bases = [x for x in bases if x != '-' and bases[x] > 1]
	return bases

def get_relpos(a):
	relpos = {}
	cur = 0
	for ix, col in enumerate(a):
		miss = (col == '-').sum() / float(len(col))
		if miss <= miss_cap:
			relpos[ix + 1] = cur + 1
			cur += 1
		else:
			relpos[ix + 1] = np.nan
	return relpos


def get_mid(relpos, coor, name):
	relmid = [relpos[x] for x in range(coor[name] - 5, coor[name] + 6) if x in relpos]
	relmid = [x for x in relmid if np.isfinite(x)]
	if len(relmid) > 0:
		relmid = int(np.mean(relmid))
	else:
		relmid = None

	return relmid


def write_aln(a, relpos, relmid, coor, name, out):
	for ix, col in enumerate(a):
		per_miss = (col == '-').sum() / float(len(col))

		out.write('%s,%s,%s,%s\n' % (name, (ix + 1) - coor[name], relpos[ix + 1] - relmid, per_miss))


aln = glob.glob(alndir + '*fasta.aln')
for f in aln:
	name = re.search('([^/]+)\.fasta', f).group(1)

	# seq with gaps
	gseq = {}
	id = ''
	f = open(f, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			gseq[id] = ''
		else:
			gseq[id] += l.rstrip()
	f.close()

	if 'Gallus_gallus' in gseq:
		del gseq['Gallus_gallus']

	seq = {}
	for c, s in gseq.items():
		seq[c] = re.sub('-', '', s)
	
	# all aln 
	a = [list(gseq[c]) for c in gseq]
	a = np.transpose(np.array(a))

	if a.shape[0] > 0:
		relpos1 = get_relpos(a)
		relmid1 = get_mid(relpos1, coor, name)
		if relmid1:
			write_aln(a, relpos1, relmid1, coor, name, fout)

	print(name)

fout.close()
