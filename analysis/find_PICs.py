import re
import os
import glob
import subprocess
import numpy as np
from collections import Counter

alndir = '/Volumes/sosi/brazil/phylogeny/alignments/'
coords = '/Volumes/sosi/brazil/metadata/alignment_to_probe_coords.csv'
out1 = '/Volumes/sosi/brazil/metadata/number_of_pics.csv'
out2 = '/Volumes/sosi/brazil/metadata/locus_lengths_pics.csv'
miss_cap = 0.33

snakes = ['Trilepida_brasiliensis', 'Typhlops_brongersmianus', 'Liotyphlops_ternetzii', 
		'Bothrops_lutzi', 'Bothrops_moojeni', 'Bothrops_pauloensis', 
		'Micrurus_brasiliensis', 'Chironius_exoletus', 'Tantilla_melanocephala',
		'Sibynomorphus_mikanii', 'Leptodeira_annulata', 'Imantodes_cenchoa', 
		'Oxyrhopus_petolarius', 'Oxyrhopus_trigeminus', 'Phimophis_guerini', 
		'Pseudoboa_neuwiedii', 'Pseudoboa_nigra', 'Taeniophallus_occipitalis', 
		'Apostolepis_cearensis', 'Apostolepis_polylepis', 'Philodryas_nattereri', 
		'Philodryas_olfersii', 'Thamnodynastes_hypoconia', 'Psomophis_joberti', 
		'Lygophis_paucidens', 'Xenodon_merremi', 'Erythrolamprus_almadensis', 
		'Erythrolamprus_poecilogyrus', 'Erythrolamprus_reginae']

coor = {}
f = open(coords, 'r')
head = f.next()
#head = f.__next__()
for l in f:
	d = re.split(',', l.rstrip())
	mid = int((int(d[1]) + int(d[2])) / 2.0) 
	coor[d[0]] = mid
f.close()

fout = open(out1, 'w')
fout.write('locus,type,abs_gap_pos,abs_ungap_pos,rel_gap_pos,rel_ungap_pos,num_bases\n')

o2 = open(out2, 'w')
o2.write('locus,type,min,max\n')

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
			relpos[ix + 1] = None
	return relpos


def get_mid(relpos, coor, name, o, type):
	relmid = [relpos[x] for x in range(coor[name] - 5, coor[name] + 6) if x in relpos]
	relmid = [x for x in relmid if x is not None]
	if len(relmid) > 0:
		relmid = int(np.mean(relmid))

		o.write('%s,%s,%s,%s\n' % (name, type, 
									min([int(x) for x in relpos.values() if x is not None]) - relmid,
									max([int(x) for x in relpos.values() if x is not None]) - relmid))
	else:
		relmid = None

	return relmid


def write_pics(a, relpos, relmid, coor, name, out, type):
	for ix, col in enumerate(a):
		bases = parse_col(col)

		if len(bases) > 1 and relpos[ix + 1]:
			out.write('%s,%s,%s,%s,%s,%s,%s,%s\n' %
					(name, type, ix + 1, ix + 1, (ix + 1) - coor[name],
					relpos[ix + 1], relpos[ix + 1] - relmid,
					(col != '-').sum()))


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

	# snakes only
	s = [list(gseq[c]) for c in gseq if c in snakes]
	s = np.transpose(np.array(s))

	if s.shape[0] > 0:
		relpos2 = get_relpos(s)
		relmid2 = get_mid(relpos2, coor, name, o2, 'SNAKE')
		if relmid2:
			write_pics(s, relpos2, relmid2, coor, name, fout, 'SNAKE')
	if a.shape[0] > 0:
		relpos1 = get_relpos(a)
		relmid1 = get_mid(relpos1, coor, name, o2, 'ALL')
		if relmid1:
			write_pics(a, relpos1, relmid1, coor, name, fout, 'ALL')

	print name

fout.close()
