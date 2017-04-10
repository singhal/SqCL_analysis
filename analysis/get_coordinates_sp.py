import re
import os
import glob
import subprocess
import numpy as np

sp = 'Colobosaura_modesta'
# sp = 'Bothrops_moojeni'
prg = '/Volumes/sosi/brazil/PRG/%s.fasta' % sp
db = '/Volumes/sosi/brazil/squamate_AHE_UCE_genes_loci.fasta'
cout1 = '/Volumes/sosi/brazil/metadata/%s_to_probe_coords.csv' % sp
cout2 = '/Volumes/sosi/brazil/metadata/%s_locus_lengths_variants.csv' % sp

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

fout1 = open(cout1, 'w')
fout1.write('locus,start,stop,frac,match\n')
fout2 = open(cout2, 'w')
fout2.write('locus,min,max\n')

seq = {}
f = open(prg, 'r')
for l in f:
	if re.search('>\S+', l):
		id = re.search('>(\S+)', l.rstrip()).group(1)
		seq[id] = ''
	else:
		seq[id] += l.rstrip()
f.close()

out = '%s_probes.blast.out' % sp
# get coordinates to probe - ungapped seq
x = subprocess.call("blastn -db %s -query %s -max_target_seqs 1 -evalue 1e-20 -outfmt 6 > %s" %
			(db, prg, out), stdout=subprocess.PIPE, shell=True)

coords = {}
f = open(out, 'r')
for l in f:
	d = re.split('\t', l.rstrip())
	probe = d[1]

	coords[d[0]] = {'pos' : [int(d[6]), int(d[7])], 'frac': int(d[3]) / p[probe], 'probe': d[1]}
f.close()

for id in seq:
	if id not in coords:
		mid = int(len(seq[id]) / 2.0)
		start = mid - 100
		end = mid + 100
		fout1.write('%s,%s,%s,NA,NA\n' % (id, start, end))
		fout2.write('%s,%s,%s\n' % (id, -1 * mid, len(seq[id]) - mid))
	else:
		fout1.write('%s,%s,%s,%.3f,%s\n' % (id, coords[id]['pos'][0], coords[id]['pos'][1], 
						coords[id]['frac'], p[coords[id]['probe']]))
		mid = int((coords[id]['pos'][0] + coords[id]['pos'][1]) / 2.0)
		fout2.write('%s,%s,%s\n' % (id, -1 * mid, len(seq[id]) - mid))
fout1.close()	
os.remove(out)
fout2.close()
