import glob
import math
import numpy as np
import operator
import os
import pandas as pd
import random
import re
import subprocess
import xlrd
import sys

# my working dir
dir = '/Users/sonal/macroevolution/brazil/probe_design/'
# the directory with all the genomes
g_dir = dir + 'squamate_genomes/'
# output dir
out_dir = dir + 'develop/'
# directory with probes
# this is chicken sequence
UCE_file = dir + 'UCE_probes/5472-probes-targeting-5060-amniote-tetrapod-uces.fasta'
# directory with AHE coordinates
# this is with reference to Anolis Carolinesis
AHE_coord_file = dir + 'AHE_probes/AHE_Snake_Ruane/AHE_Snake_ProbeRegionCoordinates.xlsx'
# repeatmasker database
repbase = '/Users/Sonal/macroevolution/brazil/probe_design/squamate_genomes/RepBase21.01.vertebrate.fasta'

# gene models
g1_file1 = dir + 'genes/sq44locus_RAXML'
g1_file2 = dir + 'genes/sq_models_RAXML'
# gene sequence
g2_file1 = dir + 'genes/Squamata_10_28.phy'
g2_file2 = dir + 'genes/Squamata_Models_modified.txt'
# gene file outgroups
g1_out = {'Gallus_gallus', 'Mus_musculus', 'Homo_sapiens', 'Alligator_mississippiensis',
				'Chelydra_serpentina', 'Crocodylus_porosus', 'Dromaius_novaehollandiae',
				'Podocnemis_expansa', 'Sphenodon_punctatus', 'Tachyglossus_aculeatus'}
g2_out = {'Sphenodon punctatus'}

# snake species
snakes = [	'Boa_constrictor', 'Crotalus_mitchelli', 'Ophiophagus_hannah',  
			'Pantherophis_guttatus', 'Python_molurus', 'Thamnophis_sirtalis',
			'Vipera_berus' ]
# probe length
probe_len = 120
max_div = 0.05
tile = 2


def write_genes(dir, genes):
	outfile =  dir + 'genes/genes.fa'
	o = open(outfile, 'w')

	for g, seq in genes.items():
		o.write('>%s\n%s\n' % (g, seq))

	o.close()

	return outfile


def get_genes(aln, models, outgroup, genes):
	seq = {}
	f = open(aln, 'r')
	head = next(f)
	for l in f:
		d = re.split('\s+', l.rstrip())
		if d[0] not in outgroup:
			seq[d[0]] = d[1]
	f.close()

	f = open(models, 'r')
	for l in f:
		if re.search('_1', l):
			gene = re.search('DNA, (\S+)_1', l).group(1)
			start = int(re.search('=\s*(\d+)', l).group(1))
			end = int(re.search('\d+\s*-\s*(\d+)', l).group(1))
			
			tmpseq = {}
			for id, s in seq.items():
				tmpseq[id] = re.sub('[-|\?]', '', s[start - 1:end]) 
			maxlength = max([len(x) for x in tmpseq.values()])

			maxids = ([x for x, y in tmpseq.items() if len(y) == maxlength])
			if 'Anolis_carolinensis' in maxids:
				winner = 'Anolis_carolinensis'
			else:
				winner = random.choice(maxids)
			
			genes['%s %s' % (gene, winner)] = tmpseq[winner]
	
	return genes


def get_genomes(g_dir):
	files = glob.glob(g_dir + '*fa')
	g = {}
	for file in files:
		sp = re.search('([A-Z|a-z]+_[A-Z|a-z]+)\.', file).group(1)
		g[ sp ] = file
		if not os.path.isfile(file + ".fai"):
			subprocess.call("samtools faidx %s" % file, shell=True)
	return g


def make_outdir(out_dir):
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)
	
	subdirs = [out_dir + 'loci/', out_dir + 'align/', out_dir + 'blat_results/']

	for subdir in subdirs:
		if not os.path.isdir(subdir):
			os.mkdir(subdir)


def do_blat(loci, g, out_dir, type):
	out_dir = out_dir + 'blat_results/'

	for sp, genome in g.items():
		outfile = out_dir + sp + '_' + type + '.out'
		if not os.path.isfile(outfile):
			subprocess.call("blat %s %s %s -out=blast8" % (genome, loci, outfile), shell=True)


def get_loci(loci):
	f = open(loci, 'r')
	
	seq = {}
	id = ''

	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()

	f.close()
	return seq


def parse_blat(g, out_dir, type, max_length):
	out_dir = out_dir + 'blat_results/'

	res = {}

	for sp in g.keys():
		outfile = out_dir + sp + '_' + type + '.out'
	
		f = open(outfile, 'r')
		for l in f:
			d = re.split('\t+', l.rstrip())
			if float(d[10]) <= 1e-10:
				
				if d[0] not in res:
					res[d[0]] = {}

				if sp not in res[d[0]]:
					direction = 'F'
					if int(d[8]) > int(d[9]):
						direction = 'R'
						
					res[d[0]][sp] = {'chr': d[1], 'dir': direction, 'loc': [[], []]}

					res[d[0]][sp]['loc'][0].append(int(d[8]))
					res[d[0]][sp]['loc'][1].append(int(d[9]))
				else:
					if d[1] == res[d[0]][sp]['chr']:
						res[d[0]][sp]['loc'][0].append(int(d[8]))
						res[d[0]][sp]['loc'][1].append(int(d[9]))
		f.close()

	for c in res:
		for sp in res[c]:
			direction = res[c][sp]['dir']
			if direction == 'F':
				start = min(res[c][sp]['loc'][0])
				end = max(res[c][sp]['loc'][1])
			else:
				start = min(res[c][sp]['loc'][1])
				end = max(res[c][sp]['loc'][0])

			# to avoid matches that span a huge distance
			if (end - start) > max_length:
				if direction == 'F':
					start = res[c][sp]['loc'][0][0]
					end = res[c][sp]['loc'][1][0]
				else:
					start = res[c][sp]['loc'][0][0]
					end = res[c][sp]['loc'][1][0]

			res[c][sp]['loc'] = [start, end]

	return(res)


def rev_comp(seq):
	# modified from http://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

	bases = list(seq)
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)

	return bases


def create_files(res, g, loci, out_dir, type, type_sp, flank):
	out_dir = out_dir + 'loci/'

	for c in loci:
		out1 = '%s%s.%s.fa' % (out_dir, c, type)
		out2 = '%s%s.%s.aln.fa' % (out_dir, c, type)

		if not (os.path.isfile(out1) or os.path.isfile(out2)):
			o = open(out1, 'w')

			o.write('>%s\n%s\n' % (type_sp, loci[c]))

			if c in res:
				for sp in res[c]:
					start = int(res[c][sp]['loc'][0]) - flank
					end = int(res[c][sp]['loc'][1]) + flank

					out = subprocess.Popen('samtools faidx %s \'%s\':%s-%s' % (g[sp], res[c][sp]['chr'], start, end), 
								shell=True, stdout=subprocess.PIPE)
					seq = ''
					for l in out.stdout:
						l = l.decode('utf-8')
						if not re.match('>', l):
							seq += l.rstrip().upper()

					if res[c][sp]['dir'] == 'R':
						seq = rev_comp(seq)

					o.write('>%s\n%s\n' %  (sp, seq))

			o.close()


def get_AHE_loci(AHE_coord_file, g):
	d = pd.read_excel(AHE_coord_file, sheet=0)
	d = d[d.anoCar2 != 'MISSING']

	out = AHE_coord_file + '.fa'
	if not os.path.isfile(out):
		o = open(out, 'w')
		
		for loc_num, chr, beg, end in zip(d.Locus, d.anoCar2, d.beg, d.end):
			# modify the chromosome names so that they match
			if re.search('chr', chr):
				chr = re.search('chr(\S+)', chr).group(1)
			if re.search('Un', chr):
				chr = re.search('Un_(\S+)', chr).group(1)
				if re.search('GL', chr):
					chr = chr + '.1'

			out = subprocess.Popen('samtools faidx %s \'%s\':%s-%s' % (g['Anolis_carolinensis'], chr, beg, end), 
										shell=True, stdout=subprocess.PIPE)
			seq = ''
			for l in out.stdout:
				l = l.decode('utf-8')
				if not re.match('>', l):
					seq += l.rstrip().upper()

			o.write('>L%s\n%s\n' % (loc_num, seq))

		o.close()

	return out


def align_seq(outdir):
	outdir = outdir + 'loci/'

	files = glob.glob(outdir + '*.fa')
	files = [file for file in files if not re.search('aln.fa', file)]

	for file in files:
		if os.path.getsize(file) > 0:
			out = file.replace('.fa', '.aln.fa')
			if not os.path.isfile(out):
				subprocess.call('mafft --maxiterate 1000 --localpair %s > %s' % (file, out), shell=True)
				os.remove(file)


def get_div(s1, s2):
	length = 0
	diff = 0

	bases = ['A', 'T', 'C', 'G'] 

	for i, j in zip(s1.upper(), s2.upper()):
		if i in bases and j in bases:
			length += 1
			if i != j:
				diff += 1

	div = np.nan
	if length > 20:
		div = diff / float(length)
	
	return div


def calc_divergence(out_dir, g, snakes):
	
	out_dir1 = out_dir + 'loci/'
	files = sorted(glob.glob(out_dir1 + '*.fa'))
	files = [file for file in files if re.search('.aln.fa', file)]

	out1 = out_dir + 'locus_divergence.csv'
	o1 = open(out1, 'w')
	o1.write('locus,type,num_inds,lizard_mean,lizard_min,lizard_max,snake_mean,snake_min,snake_max\n')

	out2 = out_dir + 'speciesbyspecies_divergence.csv'
	o2 = open(out2, 'w')
	o2.write('locus,type,spcomp,div\n')

	for file in files:
		seq = get_loci(file)

		name = re.split('\.', file.replace(out_dir1, ''))
		locus = name[0]
		type = name[1]

		snake_div = []
		lizard_div = []

		inds = sorted(list(seq.keys()))

		for ix, id1 in enumerate(inds):
			for id2 in inds[(ix + 1):]:
				if id1 in g and id2 in g:
					seq_div = get_div(seq[id1], seq[id2])

					o2.write('%s,%s,%s,%.3f\n' % (locus, type, '-'.join(sorted([id1, id2])), seq_div))

					if id1 in snakes and id2 in snakes:
						snake_div.append(seq_div)
					elif id1 not in snakes and id2 not in snakes:
						lizard_div.append(seq_div)

		if len(lizard_div) == 0:
			lizard_div.append(np.nan)
		if len(snake_div) == 0:
			snake_div.append(np.nan)

		o1.write('%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n' % (	locus, type, len(inds) - 1, 
																np.mean(lizard_div), np.min(lizard_div),
																np.max(lizard_div), np.mean(snake_div), 
																np.min(snake_div), np.max(snake_div) ) )
	o1.close()
	o2.close()


def bp_div(bp):

	comp = 0
	diff = 0
	bases = ['A', 'C', 'T', 'G']

	for ix, bp1 in enumerate(bp):
		for bp2 in bp[(ix + 1):]:
			if bp1 in bases and bp2 in bases:
				comp += 1
				if bp1 != bp2:
					diff += 1

	div = np.nan
	if comp > 0:
		div = (comp - diff) / float(comp)

	return div


def print_UCE(out_dir, g, snakes, locus, sp, tipdiv):
	
	out_dir1 = out_dir + 'loci/'
	files = sorted(glob.glob(out_dir1 + '*.fa'))
	files = [file for file in files if re.search('%s.aln.fa' % locus, file)]

	print('locus,chicken_div1,chicken_div2,full_div,ref_div,length')
	out = '/Users/Sonal/macroevolution/brazil/probe_design_v2/UCE_loci.fa'
	o = open(out, 'w')

	for file in files:
		seq = get_loci(file)
		name = re.search('(uce-\d+_p\d+)', file).group(1)

		ref = seq[sp]
		# get rid of reference
		seq1 = dict(seq)
		del seq1[sp]

		seq2 = [list(s.upper()) for s in seq1.values()]
		seq2 = list(map(list, zip(*seq2)))

		div = [bp_div(bp) for bp in seq2]
		div = pd.rolling_mean(pd.Series(div).interpolate(), 5).tolist()

		start_trim = 0
		end_trim = len(div)
		for ix, val in enumerate(reversed(div[:len(div) // 2])):
			if val <= tipdiv:
				start_trim = len(div) // 2 - ix
				break
		for ix, val in enumerate(div[len(div) // 2:]):
			if val <= tipdiv:
				end_trim = ix + len(div) // 2
				break

		ref_start = 0
		ref_end = len(ref)
		for ix, val in enumerate(ref[:len(ref) // 2]):
			if val != '-':
				ref_start = ix
				break
		for ix, val in enumerate(reversed(ref[len(ref) // 2:])):
			if val != '-':
				ref_end = len(ref) - ix
				break
		
		if ref_start < start_trim:
			start_trim = ref_start
		if ref_end > end_trim:
			end_trim = ref_end

		ref = ''.join(ref)[start_trim:end_trim]
		div = {}
		for id, s in seq.items():
			sdiv = round(get_div(seq[id][start_trim:end_trim], ref), 3)
			div[id] = sdiv

		trimseq = {}
		for id, s in seq.items():
			s = s[start_trim:end_trim].upper()
			s = re.sub('[?|-]', '', s)
			s = re.sub('^N+', '', s)
			s = re.sub('N+$', '', s)

			if len(s) > 110:
				trimseq[id] = s

		ids1 = ['Anolis_carolinensis', 'Ophisaurus_gracilis']
		ids2 = ['Python_molurus', 'Pantherophis_guttatus', 'Thamnophis_sirtalis', 'Boa_constrictor',
		        'Ophiophagus_hannah', 'Crotalus_mitchellii', 'Vipera_berus']
		ids1 = [id for id in ids1 if id in trimseq]
		ids2 = [id for id in ids2 if id in trimseq]

		if len(ids1) > 0 and len(ids2) > 0:
			id1 = ids1[0]
			id2 = ids2[0]
		else:
			if 'chicken' in trimseq:
				id1 = id2 = 'chicken'
			else:
				if len(trimseq) > 1:
					(id1, id2) = random.sample(list(trimseq.keys()), 2) 
				else:
					id1 = id2 = list(trimseq.keys())[0]

		tmpdiv1 = round(get_div(seq[id1][start_trim:end_trim], seq[id2][start_trim:end_trim]), 3)
		tmpdiv2 = round(get_div(seq[id1][ref_start:ref_end], seq[id2][ref_start:ref_end]), 3)
		
		print('%s,%s,%s,%s,%s,%s,%s' % (name, len(ref[ref_start:ref_end]), div[id1], div[id2], tmpdiv1, tmpdiv2, (len(trimseq[id1]) * 0.5 + len(trimseq[id2]) * 0.5)))
		o.write('>%s_%s_UCE\n%s\n' % (id1, name, trimseq[id1]))
		if id1 == id2:
			o.write('>%s2_%s_UCE\n%s\n' % (id2, name, trimseq[id2]))
		else:
			o.write('>%s_%s_UCE\n%s\n' % (id2, name, trimseq[id2]))


# http://stackoverflow.com/questions/1285434/efficient-algorithm-for-string-concatenation-with-overlap
def concat(*args):
    result = ''
    for arg in args:
        result = _concat(result, arg)
    return result


def _concat(a, b):
    la = len(a)
    lb = len(b)
    for i in range(la):
        j = i
        k = 0
        while j < la and k < lb and a[j] == b[k]:
            j += 1
            k += 1
        if j == la:
            n = k
            break
    else:
        n = 0
    return a + b[n:]

def get_id(id):
	newid = ''
	if re.search('_genes', id):
		newid = re.search('([^_]+)_genes', id).group(1)
	elif re.search('L\d+', id):
		newid = re.search('(L\d+)', id).group(1)
	elif re.search('uce-\d+', id):
		newid = re.search('(uce-\d+)', id).group(1)
	if not newid:
		sys.exit("No id for %s." % id)
	return newid

def run_blat(out):
	blatout = '%s_selfBlat' % out
	if not os.path.isfile(blatout):
		subprocess.call('blat %s %s %s -out=blast8' % (out, out, blatout), shell=True)

	comp = {}

	f = open(blatout, 'r')
	for l in f:
		d = re.split('\t', l.rstrip())

		id1 = get_id(d[0])
		id2 = get_id(d[1])
		
		if id1 != id2 and float(d[10]) < 1e-20:
			sort = sorted([id1, id2])
			if sort[0] not in comp:
				comp[sort[0]] = {}
			if sort[1] not in comp[sort[0]]:
				comp[sort[0]][sort[1]] = []
			comp[sort[0]][sort[1]].append(d)
	f.close()

	for i in comp:
		for j in comp[i]:
			print(i)
			print(j)
			for x in comp[i][j]:
				print('\t'.join(x))
			print('***')


def combine_files(dir, files):
	out = 'all_loci_raw.fa'

	seq = {}

	for file in files:
		f = open(file, 'r')
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l).group(1)
				if id not in seq:
					seq[id] = next(f).rstrip()
				else:
					seq[id + '_dup'] = next(f).rstrip()

	o = open(out, 'w')
	for id, s in seq.items():
		o.write('>%s\n%s\n' % (id, s))
	o.close()

	run_blat(out)

	return out


def assemble_UCEs(file):
	seqs = {}

	f = open(file, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)_p\d+_UCE', l).group(1)
			if id not in seqs:
				seqs[id] = {}
			fullid = re.search('>(\S+)', l.rstrip()).group(1)
			seqs[id][fullid] = next(f).rstrip()
	f.close()

	newseq = {}
	for id in seqs:
		if len(seqs[id]) > 1:
			out = id + '.fa'
			f = open(out, 'w')
			for id1, s in seqs[id].items():
				f.write('>%s\n%s\n' % (id1, s))
			f.close()

			subprocess.call('cap3 %s' % out, shell=True)

			sing = get_loci(out + '.cap.singlets')
			contig = get_loci(out + '.cap.contigs')

			for id1, s in sing.items():
				newseq[id1] = s
			track = 0
			for id1, s in contig.items():
				newseq[id + '-c%s_UCE' % track] = s
				track = track + 1

			subprocess.call('rm %s*' % out, shell=True)
		else:
			for id1 in seqs[id]:
				newseq[id1] = seqs[id][id1]

	newfile = file.replace('.fa', '.assembled.fa') 
	o = open(newfile, 'w')
	for id, s in newseq.items():
		if re.search('chicken', id):
			id = re.sub('chicken', "Gallus_gallus", id)
		elif re.search('_p\d', id):
			id = re.sub('_p', '-p', id)
		o.write('>%s\n%s\n' % (id, s))
	o.close()


def drop_check(dir, out, dropfile):
	f = open(dropfile, 'r')
	head = next(f)
	matches = {}
	for l in f:
		d = re.split('\s+', l.rstrip())
		matches[d[1]] = 1
	f.close()

	seq = get_loci(out)

	out = out.replace('_raw', '_filtered')
	o = open(out, 'w')
	for id, s in seq.items():
		drop = False
		for match in matches:
			if re.search('[_|-]' + match, id):
				drop = True
		if not drop:
			o.write('>%s\n%s\n' % (id, s))
	o.close()

	run_blat(out)


def gc_length(out):
	seq = get_loci(out)

	tot_length = 0
	gc = 0

	out = out.replace('.fa', '.summary.csv')
	o = open(out, 'w')
	o.write('locus,length,GC\n')
	for id, s in seq.items():
		tot_length += len(s)

		seq_gc = s.count('C') + s.count('G')
		gc += seq_gc
		seq_gc = seq_gc / float(len(s))

		o.write('%s,%s,%.3f\n' % (id, len(s), seq_gc))
	o.close()

	print('number loci: %s' % len(seq))
	print('tot length: %s' % tot_length)
	print('GC: %.3f' % (gc / float(tot_length)))


def rep_masked(file, repbase, stem, new):
	out = file + '_repeatBlat'

	if not os.path.isfile(out):
		subprocess.call('blat %s %s %s -out=blast8' % (repbase, file, out), shell=True)

	to_drop = {}
	o = open(out, 'r')
	for l in o:
		d = re.split('\s+', l.rstrip())
		if float(d[10]) < 1e-20:
			to_drop[d[0]] = 1
	o.close()

	seq = get_loci(file)
	out = file.replace(stem, new)
	o = open(out, 'w')
	for id, s in seq.items():
		if id not in to_drop:
			o.write('>%s\n%s\n' % (id, s))
	o.close()

	return out


def make_probes(out, probe_len, tile):
	seq = get_loci(out)

	out = out.replace('.fa', '.probes.fa')
	o = open(out, 'w')
	for id, s in seq.items():
		if len(s) <= 120:
			o.write('>%s_probe0\n%s\n' % (id, s))
		elif len(s) < 160:
			starts = [0, len(s) - 120]
			probes = [s[start: probe_len + start] for start in starts]

			for ix, probe in enumerate(probes):
				o.write('>%s_probe%s\n%s\n' % (id, ix, probe))
		else:
			num_probes = round(len(s) / (probe_len / float(tile))) - 1
			end = len(s) - probe_len
			starts = list(range(0, end + 1, int(end / (num_probes - 1))))

			probes = [s[start: probe_len + start] for start in starts]

			for ix, probe in enumerate(probes):
				o.write('>%s_probe%s\n%s\n' % (id, ix, probe))

	return out

def map_to_genomes(file, g):
	out = file + '_genomeBlat'

	if not os.path.isfile(out):
		subprocess.call('blat %s %s %s -out=blast8' % (g['Anolis_carolinensis'], file, out), shell=True)

	matches = {}
	o = open(out, 'r')
	for l in o:
		d = re.split('\s+', l.rstrip())
		if float(d[10]) < 1e-20:
			if d[0] not in matches:
				matches[d[0]] = []
			matches[d[0]].append(d)
	o.close()

	seq = get_loci(file)
	out = file.replace('final', 'final_filtered')
	o = open(out, 'w')
	for id, s in seq.items():
		if id not in matches:
			print(id)
			o.write('>%s\n%s\n' % (id, s))
		else:
			if len(matches[id]) < 2:
				o.write('>%s\n%s\n' % (id, s))
	o.close()

	return out

### prep work
# make_outdir(out_dir)
g = get_genomes(g_dir)

### genes
# USELESS TOO DIVERGENT
# genes = {}
# genes = get_genes(g1_file1, g1_file2, g1_out, genes)
# genes = get_genes(g2_file1, g2_file2, g2_out, genes)
# g_file = write_genes(dir, genes)
# do_blat(g_file, g, out_dir, 'genes')
# res = parse_blat(g, out_dir, 'genes', 5000)
# loci = get_loci(g_file)
# create_files(res, g, loci, out_dir, 'genes', 'lizard', 0)

### AHEs
# USELESS TOO DIVERGENT
# AHE_file = get_AHE_loci(AHE_coord_file, g)
# do_blat(AHE_file, g, out_dir, 'AHE')
# res = parse_blat(g, out_dir, 'AHE', 2000)
# # get the original loci
# loci = get_loci(AHE_file)
# # create files that will then be used for alignments
# create_files(res, g, loci, out_dir, 'AHE', 'lizard', 0)

### UCEs
# do_blat(UCE_file, g, out_dir, 'UCE')
# res = parse_blat(g, out_dir, 'UCE', 300)
# # get the original loci
# loci = get_loci(UCE_file)
# # create files that will then be used for alignments
# create_files(res, g, loci, out_dir, 'UCE', 'chicken', 60)

### align the sequences
# align_seq(out_dir)

### evaluate the loci
# get divergence among loci
# calc_divergence(out_dir, g, snakes)
# figure out conservation and extend UCEs if appropriate
# print_UCE(out_dir, g, snakes, 'UCE', 'chicken', 0.75)
# assemble across UCEs; there are multiple sequences per locus
# assemble_UCEs('/Users/Sonal/macroevolution/brazil/probe_design_v2/UCE_loci.fa')
# combine across loci types and self blat
out = combine_files(dir, ['AHE_loci.fa', 'UCE_loci.assembled.fa', 'gene_loci.fa'])
# drop repeated loci (identified by hand) and check again
out = 'all_loci_raw.fa'
drop_check(dir, out, 'to_drop.txt')
# check repeat masked
# out = rep_masked(out, repbase, 'filtered', 'final')
# check GC, length
# gc_length(out)

# out = '/Users/sonal/macroevolution/brazil/probe_design/develop/all_loci_final.fa'
# out = make_probes(out, probe_len, tile)
# out = map_to_genomes(out, g)
# out = rep_masked(out, repbase, 'final_filtered', 'final_filtered_rm')
# gc_length(out)