import re
import subprocess
import os
import argparse
import pandas as pd
import gzip

dir = '/Volumes/heloderma3/brazil/'
bwa = '/Volumes/heloderma4/sonal/bin/bwa-0.7.12/bwa'
samtools = '/Volumes/heloderma4/sonal/bin/samtools'
bcftools = '/Volumes/heloderma4/sonal/bin/bcftools/bcftools'
readdir = '/Volumes/heloderma3/brazil/trim_reads/'
refdir = '/Volumes/heloderma3/brazil/mitogenomes/'
threads = 10
min_depth = 5

parser = argparse.ArgumentParser(description="run for the species",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                )
parser.add_argument('--ind', type=str, default=None, help='individual for which to run.')
args = parser.parse_args()
ind = args.ind

# get the reads
r1 = os.path.join(readdir, '%s_R1.final.fq.gz' % ind)
r2 = os.path.join(readdir, '%s_R2.final.fq.gz' % ind)
ru = os.path.join(readdir, '%s_unpaired.final.fq.gz' % ind)

# get the ref seq
seq = os.path.join(refdir, '%s_mitogenome.fa' % ind) 

# index ref seq
subprocess.call("%s index %s" % (bwa, seq), shell=True)

# align reads
out1a = os.path.join(refdir, '%s_1.sam' % ind) 
out1b = os.path.join(refdir, '%s_2.sam' % ind)
out2a = os.path.join(refdir, '%s_1.bam' % ind)
out2b = os.path.join(refdir, '%s_2.bam' % ind)
out3a = os.path.join(refdir, '%s_1.sorted.bam' % ind)
out3b = os.path.join(refdir, '%s_2.sorted.bam' % ind)
final = os.path.join(refdir, '%s.sorted.bam' % ind)
vcf = os.path.join(refdir, '%s.vcf.gz' % ind)
seqout1 = os.path.join(refdir, '%s_mito_consensus.fa' % ind)
tmpdir = os.path.join(refdir, '%s' % ind)
if not os.path.isdir(tmpdir):
 	os.mkdir(tmpdir)

subprocess.call("%s mem -t %s %s %s %s > %s" % (bwa, threads, seq, r1, r2, out1a), shell=True)
subprocess.call("%s mem -t %s %s %s > %s" % (bwa, threads, seq, ru, out1b), shell=True)
# turn into bam
subprocess.call("%s view -@ %s -b %s > %s" % (samtools, threads, out1a, out2a), shell=True)
subprocess.call("%s view -@ %s -b %s > %s" % (samtools, threads, out1b, out2b), shell=True)
# sort bams
subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % (samtools, threads, out3a, tmpdir, out2a), shell=True)
subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % (samtools, threads, out3b, tmpdir, out2b), shell=True)
# merge bams
subprocess.call("%s merge %s %s %s" % (samtools, final, out3a, out3b), shell=True)
# call snps
subprocess.call("%s mpileup -ugf %s %s | %s call -mO z -o %s" % (samtools, seq, final, bcftools, vcf), shell=True)

f = gzip.open(vcf, 'r')
seq = {}
for l in f:
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		pos = int(d[1])
		ref = d[3]
		alt = d[4]
		geno = re.search('(\S/\S)', d[9]).group(1)
		cov = int(re.search('DP=(\d+)', d[7]).group(1))

		seq[pos] = {'ref': ref, 'alt': alt, 'geno': geno, 'cov': cov}
f.close()

o = open(seqout1, 'w')
seqmin = min(seq.keys())
seqmax = max(seq.keys())
new = ''
for i in range(seqmin, seqmax+1):
	if i in seq:
		if seq[i]['cov'] >= min_depth:
			if seq[i]['geno'] == '0/0':
				new += seq[i]['ref']
			elif seq[i]['geno'] == '1/1':
				if seq[i]['alt'] != '.':
					new += seq[i]['alt']
				else:
					new += seq[i]['ref']
			else:
				new += 'N'
		else:
			new += 'N'
	else:
		new += 'N'
o.write('>%s\n%s\n' % (ind, new))
o.close()

# clean things up
[os.remove(x) for x in [final, out1a, out1b, out2a, out2b, out3a, out3b, vcf]]
os.rmdir(tmpdir)
