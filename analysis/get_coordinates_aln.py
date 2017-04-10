import re
import os
import glob
import subprocess
import numpy as np

alndir = '/Volumes/sosi/brazil/phylogeny/alignments/'
db = '/Volumes/sosi/brazil/squamate_AHE_UCE_genes_loci.fasta'
cout = '/Volumes/sosi/brazil/metadata/alignment_to_probe_coords.csv'

p = {}
id = ''
d = open(db, 'r')
for l in d:
	if re.search('>', l):
		id = re.search('>(\S+)', l).group(1)
		p[id] = ''
	else:
		p[id] += l.rstrip()
d.close()

for c, s in p.items():
	p[c] = len(s)


fout = open(cout, 'w')
fout.write('locus,start,stop\n')

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

	# print out seq without gaps
	out = os.path.join(alndir, '%s.tmp.fa' % (name))
	o = open(out, 'w')
	seq = {}
	for c, s in gseq.items():
		seq[c] = re.sub('-', '', s)
		o.write('>%s\n%s\n' % (c, seq[c]))
	o.close()

	# get coordinates to probe - ungapped seq
	x = subprocess.Popen("blastn -db %s -query %s -max_target_seqs 1 -evalue 1e-20 -outfmt 6" %
				(db, out), stdout=subprocess.PIPE, shell=True)
	x = [l.decode('ascii').rstrip() for l in x.stdout]
	coords = {}
	probe = ''
	for l in x:
		d = re.split('\t', l)
		coords[d[0]] = {'pos': [int(d[6]), int(d[7])], 'frac': int(d[3]) / p[d[1]] }
		probe = d[1]	

	starts = []
	ends = []
	# convert seq to gapped coordinates
	for ind in coords:
		pos = coords[ind]['pos']
		real = 0
		aln = 0
		mappos = {}
		for i in gseq[ind]:
			aln += 1
			if i != '-':
				real += 1
				mappos[real] = aln
		starts.append(mappos[pos[0]])
		ends.append(mappos[pos[1]])

	start = int(np.median(starts))
	end = int(np.median(ends))
	frac = np.median([coords[x]['frac'] for x in coords])

	fout.write('%s,%s,%s,%.3f,%s\n' % (name, start, end, frac, p[probe]))

	os.remove(out)
fout.close()
