import re
import os
import pandas as pd
import subprocess
import glob

dir = '/Volumes/heloderma3/brazil/'
ref = os.path.join(dir, 'ref', 'anolis_mt_proteins.fa')
ch = os.path.join(dir, 'ref', 'chicken_mt_cds.fa')

s = os.path.join(dir, 'samples.csv')
s = pd.read_csv(s)
outdir = os.path.join(dir, 'mito_aln')

def rev_comp(s):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G', 
		    'a':'t','t':'a','g':'c','c':'g'}
	return "".join([seq_dict[base] if base in seq_dict else base for base in reversed(s)])


for ind in s['sample']:
	# read in the sequence
	seqfile = os.path.join(dir, 'mitogenomes', '%s_mitogenome.fa' % ind)
	if os.path.isfile(seqfile):
		f = open(seqfile, 'r')
		head = f.next()
		seq = ''
		for l in f:
			seq += l.rstrip()
		f.close()

		# do the exonerate
		out = os.path.join(outdir, '%s_match' % ind)
		if not os.path.isfile(out):
			subprocess.call("exonerate --geneticcode 2 --showalignment no "
        	        	        "--showvulgar no --showcigar yes --model protein2genome "
        	        	        "%s %s > %s" % (ref, seqfile, out), shell=True)

		coords = {}
		o = open(out, 'r')
		for l in o:
			if re.search('^cigar', l):
				d = re.split('\s+', l.rstrip())
				loc = [int(d[6]), int(d[7])]
				loc = sorted(loc)

				if d[1] in coords:
					if coords[d[1]]['score'] < int(d[9]):
						coords[d[1]] = {'loc': loc, 'dir': d[8], 'score': int(d[9])}
				else:
					coords[d[1]] = {'loc': loc, 'dir': d[8], 'score': int(d[9])}
		o.close()

		# write out the seq
		for locus in coords:
			out = os.path.join(outdir, '%s.fa' % locus)
			o = open(out, 'a')
			subseq = seq[coords[locus]['loc'][0]:(coords[locus]['loc'][-1] + 1)]
			if coords[locus]['dir'] == '-':
				subseq = rev_comp(subseq)
			
			o.write('>%s\n%s\n' % (ind, subseq))
			o.close()


def read_seq(seqfile):
	s = open(seqfile, 'r')
	seq = {}
	id = ''
	for l in s:
	        if re.search('>', l.rstrip()):
			id = re.search('>(\S+)', l.rstrip()).group(1)
                	seq[id] = ''
		else:
			seq[id] += l.rstrip()
	s.close()
	return seq


# write the outgroup
og = read_seq(ch)

aln = []
for id, seq in og.items():
	out = os.path.join(outdir, '%s.fa' % id)
	o = open(out, 'a')
	o.write('>chicken\n%s\n' % (seq))
	o.close()
	aln.append(out)

# align the seq
outaln = []
for f in aln:
	out = f + '.aln'
	subprocess.call("mafft --auto %s > %s" % (f, out), shell=True)
	outaln.append(out)	

# concatenate the seq
all = {}
loci = {}
for f in outaln:
	name = re.search('([^/]+).fa', f).group(1)
	tmp = read_seq(f)
	if len(tmp) > 0:
		loci[name] = len(tmp[tmp.keys()[0]])
		for ind, s in tmp.items():
			if ind not in all:
				all[ind] = {}
			all[ind][name] = s

out = os.path.join(outdir, 'brazil_aligned_mitogenes.fasta')
o = open(out, 'w')
for ind in all:
	seq = ''
	for locus in sorted(loci.keys()):
		if locus in all[ind]:
			seq += all[ind][locus]
		else:
			seq += '-' * loci[locus]
	o.write('>%s\n%s\n' % (ind, seq))
o.close()
