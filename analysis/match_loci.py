import re
import os
import subprocess

outdir = '/Volumes/sosi/brazil/anolis_genome/'
prg = '/Volumes/sosi/brazil/PRG/Anolis_brasiliensis.fasta'
gen = '/Volumes/sosi/brazil/anolis_genome/anoCar2.fa'
blastn = '/Volumes/heloderma4/sonal/bin/blastn'

# do the blast
sp1 = re.search('.*/([^/]+).fa', prg).group(1)
sp2 = re.search('.*/([^/]+).fa', gen).group(1)
out = os.path.join(outdir, '%s_%s.blast.out' % (sp1, sp2))

if not os.path.isfile(out):
	subprocess.call("%s -db %s -query %s -out %s -outfmt 6 -max_target_seqs 1 -evalue 1e-20" % 
			(blastn, gen, prg, out), shell=True)

# get seq lengths
seq = {}
id = ''
f = open(prg, 'r')
for l in f:
	if re.search('>', l):
		id = re.search('>(\S+)', l).group(1)
		seq[id] = ''
	else:
		seq[id] += l.rstrip()
f.close()
for id, s in seq.items():
	seq[id]= len(s)

# parse the blast
out2 = os.path.join(outdir, '%s_%s.locations' % (sp1, sp2))
match = {}
f = open(out, 'r')
for l in f:
	d = re.split('\t', l.rstrip())
	if d[0] not in match:
		start = int(d[6])
		end = int(d[7])
		chr = d[1]
		pos = sorted([int(d[8]), int(d[9])])

		match[d[0]] = {'chr': chr, 'pos': pos}
f.close()

o = open(out2, 'w')
o.write('locus,chr,start,end\n')
for c in match:
	o.write('%s,%s,%s\n' % (c, match[c]['chr'], ','.join([str(x) for x in match[c]['pos']])))
o.close()
