import pandas as pd
import subprocess
import os
import random

outdir = '/scratch/drabosky_flux/sosi/brazil/'
file = '/scratch/drabosky_flux/sosi/brazil/samples.csv'
d = pd.read_csv(file)

lineages = ['Bothrops_moojeni', 'Colobosaura_modesta']

d = d.ix[d.lineage.isin(lineages)]

for ix, row in d.iterrows():
	for count in ['500000', '1000000', '1500000', '2000000']:
		seed = random.randint(1,1000)
		out1 = '/scratch/drabosky_flux/sosi/brazil_sub/raw_reads/%s_%s_R1.fastq' % (row['sample'], int(count) * 2)
		out2 = '/scratch/drabosky_flux/sosi/brazil_sub/raw_reads/%s_%s_R2.fastq' % (row['sample'], int(count) * 2)
		subprocess.call("/home/sosi/bin/seqtk/seqtk sample -s %s %s %s > %s" % (seed, row['read1'], count, out1), shell=True)
		subprocess.call("gzip %s" % (out1), shell=True)
		subprocess.call("/home/sosi/bin/seqtk/seqtk sample -s %s %s %s > %s" % (seed, row['read2'], count, out2), shell=True)
                subprocess.call("gzip %s" % (out2), shell=True)
	
