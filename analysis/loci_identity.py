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

match_len = 0
# parse the blast
out2 = os.path.join(outdir, '%s_%s.bed' % (sp1, sp2))
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
		match_len += (end - start)
f.close()

tot_len = sum(seq.values())
print "match length: " + str(match_len)
print "total length: " + str(tot_len)

# write the first bed file
o = open(out2, 'w')
for c in match:
	o.write('%s\t%s\t%s\n' % (match[c]['chr'], match[c]['pos'][0] - 1, match[c]['pos'][1]))
o.close()

# write the second bed file
out3 = os.path.join(outdir, 'ensExons.bed')
o = open(out3, 'w')
out3c = os.path.join(outdir, 'ensCDS.bed')
o1 = open(out3c, 'w')
f = open(os.path.join(outdir, 'ensGene.txt'), 'r')

for l in f:
	d = re.split('\t', l.rstrip())
	starts = re.split(',', d[9])
	ends = re.split(',', d[10])
	starts = [int(x) for x in starts if x != '']
	ends = [int(x) for x in ends if x != '']

	o1.write('%s\t%s\t%s\n' % (d[2], int(d[6]), int(d[7]) + 1))

	for i, j in zip(starts, ends):
		o.write('%s\t%s\t%s\n' % (d[2], i, j+1))
o.close()
o1.close()

# sort the bed files
out2s = re.sub('.bed', '.sorted.bed', out2)
out3s = re.sub('.bed', '.sorted.bed', out3)
out3cs = re.sub('.bed', '.sorted.bed', out3c)

subprocess.call("/Volumes/heloderma4/sonal/bin/bedtools2/bin/sortBed -i %s > %s" % (out2, out2s), shell=True)
subprocess.call("/Volumes/heloderma4/sonal/bin/bedtools2/bin/sortBed -i %s > %s" % (out3, out3s), shell=True)
subprocess.call("/Volumes/heloderma4/sonal/bin/bedtools2/bin/sortBed -i %s > %s" % (out3c, out3cs), shell=True)

out4 = os.path.join(outdir, '%s_%s_probe_gene.intersect.bed' % (sp1, sp2))
subprocess.call("/Volumes/heloderma4/sonal/bin/bedtools2/bin/intersectBed -a %s -b %s > %s" % (out2s, out3s, out4), shell=True)

out5 = os.path.join(outdir, '%s_%s_probe_CDS.intersect.bed' % (sp1, sp2))
subprocess.call("/Volumes/heloderma4/sonal/bin/bedtools2/bin/intersectBed -a %s -b %s > %s" % (out2s, out3cs, out5), shell=True)

exons = 0
f = open(out4, 'r')
for l in f:
	d = re.split('\t', l.rstrip())
	exons += (int(d[2]) - int(d[1]))
f.close()

print "in exons length: " + str(exons)

cds = 0
f = open(out5, 'r')
for l in f:
        d = re.split('\t', l.rstrip())
        cds += (int(d[2]) - int(d[1]))
f.close()

print "in CDS length: " + str(cds)
