import argparse
import re
import glob

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--tag', help='Loci to target', required=True)
parser.add_argument('--dir', help='Directory with alignments', required=True)
args = parser.parse_args()

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

def get_div(seq1, seq2):
	denom = 0
	diff = 0
	bases = ['A', 'T', 'C', 'G']
	for x, y in zip(seq1.upper(), seq2.upper()):
		if x in bases and y in bases:
			denom += 1
			if x != y:
				diff += 1
	div = diff / denom
	return div

alns = glob.glob(args.dir + '/' + '*%s*gb' % args.tag)

seqs = {}
ids = []
loci = {}
for aln in alns:
	loci_name = re.search('.*\/(\S+)\.fasta', aln).group(1)
	seq, ids = get_seq(aln, ids)

	loci_len = len(seq[list(seq.keys())[0]])
	loci[loci_name] = loci_len
	seqs[loci_name] = seq

ids = [id for id in ids if id != outgroup]

all_seq = {}
for ind in ids:
	all_seq[ind] = ''
	for locus in sorted(list(loci.keys())):
		if ind not in seqs[locus]:
			seq = '-' * loci[locus]
		else:
			seq = seqs[locus][ind]
		all_seq[ind] += seq

o = open("%s_divergence.csv" % args.tag, 'w')
o.write('ind1,ind2,divergence\n')
inds = sorted(ids)
for ix, ind1 in enumerate(inds):
	for ind2 in inds[(ix+1):]:
		div = get_div(all_seq[ind1], all_seq[ind2])
		o.write('%s,%s,%.6f\n' % (ind1, ind2, div))
o.close()