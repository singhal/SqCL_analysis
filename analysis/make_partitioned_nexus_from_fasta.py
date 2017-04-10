import re
import glob
import pandas as pd
import random
import os
import argparse

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--tag', help='Description for foo argument', required=True)
parser.add_argument('--num', help='Description for bar argument', required=True)
args = parser.parse_args()

tag = args.tag
filenum = args.num
num_loci = 100

d = pd.read_csv("/Users/Sonal/Desktop/div_dating/locus_data.csv")
d = d[d.number > 42]
d = d[d.type.str.contains(tag)]
if d.shape[0] < num_loci:
	locinames = d['type'].tolist()
else:
	locinames = random.sample(d['type'].tolist(), num_loci)

indir = '/Users/Sonal/Desktop/div_dating/alignments/'
outgroup = 'Gallus_gallus'

def get_seq(seqfile, ids):
	seq = {}
	id = ''
	f = open(seqfile, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			id = re.sub('_R_', '', id)
			seq[id] = ''
			if id not in ids:
				ids.append(id)
		else:
			seq[id] += l.rstrip()

	f.close()

	for id, s in seq.items():
		seq[id] = re.sub('\s+', '', s)

	return seq, ids

ids = []
seqs = {}
tot_length = 0
parts = []
loci = {}

for uce in locinames:
	uce = os.path.join(indir, '%s.fasta.aln-gb' % uce)
	seq, ids = get_seq(uce, ids)
	ucename = re.sub('^.*\/', '', uce)
	ucename = re.sub('\..*$', '', ucename)
	seqs[ucename] = seq

	loci_len = len(seq[list(seq.keys())[0]])
	tot_length += loci_len
	loci[ucename] = loci_len

ids = [id for id in ids if id != outgroup]
num_inds = len(ids)

all_seq = {}
for ind in ids:
	all_seq[ind] = ''
	for locus in sorted(list(loci.keys())):
		if ind not in seqs[locus]:
			seq = '-' * loci[locus]
		else:
			seq = seqs[locus][ind]
		all_seq[ind] += seq

o = open("/Users/sonal/Desktop/%s%s.nex" % (tag, filenum), 'w')
o.write('#nexus\n')
o.write('begin data;\n')
o.write('\tdimensions ntax=%s nchar=%s;\n' % (num_inds, tot_length))
o.write('\tformat datatype=dna missing=? gap=-;\n')
o.write('\tmatrix\n')
for ind, seq in all_seq.items():
	indname = ind
	if len(indname) > 29:
		indname = indname[0:29]
	indname = indname + ' ' * (30 - len(indname))
	o.write('\t\t%s%s\n' % (indname, seq))
o.write('\t;\n')
o.write('end;\n')

o.write('begin sets;\n')
cur_pos = 1
for locus in sorted(list(loci.keys())):
	last_pos = loci[locus] + cur_pos - 1
	locname = re.sub('-', '_', locus)
	o.write('\tCHARSET %s = %s-%s;\n' % (locname, cur_pos, last_pos))
	cur_pos = loci[locus] + cur_pos
o.write('END;\n')
o.close()
