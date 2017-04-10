import re
import subprocess
import glob
import os
import pandas as pd
import numpy as np

dir = '/scratch/drabosky_flux/sosi/brazil/'
pr = os.path.join(dir, 'squamate_AHE_UCE_genes_loci.fasta')
org = '/nfs/turbo/sosi/brazil/baits.fasta'
s = os.path.join(dir, 'samples.csv')
s = pd.read_csv(s)

def get_seq(seqfile):
	id = ''
	seq = {}
	f = open(seqfile, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()
	
	for id, s in seq.items():
		seq[id] = s.upper()

	return seq

print("Getting probe data ...")
# get probe data
pr = get_seq(pr)
org = get_seq(org)
stats = {}
l2p = {}
for x in pr:
	name = re.search('^([^_]+)', x).group(1)
	loctype = re.search('(\S+)-', name).group(1)
	stats[name] = {'type': loctype}
	if name not in l2p:
		l2p[name] = []
	l2p[name].append(x) 


print("Getting seq data ...")
# get ind data
all_seq = {}
for ind in s['sample']:
	seqfile = os.path.join(dir, 'trinity_assembly', '%s.fasta' % ind)
	seq = get_seq(seqfile)
	
	m = os.path.join(dir, 'matches', '%s_matches.csv' % ind)
	m = pd.read_csv(m)
	ids = {}
	for ix, row in m.iterrows():
		if row['status'] in ['easy_recip_match', 'complicated_recip_match']:
			ids[row['contig']] = row['match']

	all_seq[ind] = {}
	for old, new in ids.items():
		all_seq[ind][new] = seq[old]


print("Getting missingness ...")
# calculate missingness
for loc in stats:
	stats[loc]['complete'] = sum([1 for ind in all_seq if loc in all_seq[ind]])

print("Getting probe stats ...")
# probe density
probes = {}
for loc in pr:
	probes[loc] = {'length': len(pr[loc]), 'num': 0, 'dens': 0}
for id in org:
	id = re.sub('_\d+$', '', id)
	probes[id]['num'] += 1
for loc in pr:
	probes[loc]['dens'] = probes[loc]['num'] * 120 / float(probes[loc]['length'])
for loc in stats:
	num = np.mean([probes[probe]['num'] for probe in l2p[loc]])
	dens = np.mean([probes[probe]['dens'] for probe in l2p[loc]])
	length = np.mean([probes[probe]['length'] for probe in l2p[loc]])
	stats[loc]['num_probes'] = num
	stats[loc]['probe_density'] = dens
	stats[loc]['probe_length'] = length
print stats

print("Getting GC data ...")
# calculate gc
for loc in stats:
	gc = []
	for ind in all_seq:
		if loc in all_seq[ind]:
			count = all_seq[ind][loc].count('G') + all_seq[ind][loc].count('C')
			gc.append( count / float(len(all_seq[ind][loc])) )
	if len(gc) > 0:
		stats[loc]['seq_GC'] = np.mean(gc)
	else:
		stats[loc]['seq_GC'] = np.nan

	gc = []
	for probe in l2p[loc]:
		count = pr[probe].count('G') + pr[probe].count('C')
		gc.append( count / float(len(pr[probe])))
	stats[loc]['probe_GC'] = np.mean(gc)

def get_div(g, s):
	diff = 0
	comp = 0
	for a, b in zip(g, s):
		if a in ['A', 'T', 'G', 'C'] and b in ['A', 'T', 'G', 'C']:
			comp += 1
			if a != b:
				diff += 1
	if comp > 0:
		diff = diff / float(comp)
	else:
		diff = np.nan	

	return diff

print("Getting divergence stats ...")
# calculate divergence
gecko = ['Gymnodactylus_amarali', 'Hemidactylus_mabouia']
snake = ['Bothrops_moojeni', 'Philodryas_olfersii', 'Liotyphlops_ternetzii']
for loc in stats:
	f = os.path.join(dir, 'phylogeny/alignments', '%s.fasta.aln' % loc)
	div = []
	if os.path.isfile(f):
		seq = get_seq(f)
		for g in gecko:
			for s in snake:
				if g in seq and s in seq:
					div.append(get_div(seq[g], seq[s]))

	if len(div) > 0:
		stats[loc]['seq_div'] = np.mean(div)
	else:
		stats[loc]['seq_div'] = np.nan

	out = '%s_tmp.fa' % loc
	o = open(out, 'w')
	for probe in l2p[loc]:
		o.write('>%s\n%s\n' % (probe, pr[probe])) 
	o.close()
	aln = '%s_tmp.fa.aln' % loc
	subprocess.call("mafft %s > %s" % (out, aln), shell=True)
	alnseq = get_seq(aln)
	
	os.remove(out)
	os.remove(aln)

	if len(alnseq) > 1:
		stats[loc]['probe_div'] = get_div(alnseq.values()[0], alnseq.values()[1])
	else:
		stats[loc]['probe_div'] = np.nan

out = os.path.join(dir, 'model.csv')
o = open(out, 'w')

vals = sorted(stats[stats.keys()[0]].keys())
o.write('locus,%s\n' % (','.join(vals)))
for loc in stats:
	tvals = [loc] + [stats[loc][val] for val in vals]
	tvals = [str(x) for x in tvals]
	o.write(','.join(tvals))
o.close()

