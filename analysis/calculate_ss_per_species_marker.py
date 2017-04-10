import argparse
import gzip
import os
import numpy as np
import pandas as pd
import re
import subprocess

'''
Sonal Singhal
created on 29 June 2016
'''


def get_args():
	parser = argparse.ArgumentParser(
		description="Calculate pi for a lineage.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
        parser.add_argument(
                '--lineage',
                type=str,
                default=None,
                help='Lineage for which to make calculations.'
                )

        # sample file
        parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )
        
        # base dir
        parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Base directory as necessary"
                     " when used with pipeline"
                )

	# output dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Directory to output pop gen stats, '
                     'only necessary if not running '
                     'in context of pipeline'
                )

	# vcfdir
	parser.add_argument(
		'--vcfdir',
		type=str,
		default=None,
		help='Directory with VCFs, '
		     'only necessary if not running '
                     'in context of pipeline'
		)

	return parser.parse_args()


def get_diversity(lineage, inds, vcf, outdir):
	all = {}

	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				if d[0] not in all:
					all[d[0]] = {'ss': 0, 'denom': 0}
				all[d[0]]['denom'] += 1
				if d[4] in ['A', 'T', 'C', 'G']:
					all[d[0]]['ss'] += 1
	f.close()
		
	out = os.path.join(outdir, '%s_segregatingsites.csv' % lineage)
	o = open(out, 'w')
	o.write('lineage,ninds,locus,locus_length,ss,ss_nuc\n')
	for locus in all:
		denom = 0
		for i in range(1, len(inds) * 2):
			denom += (1 / float(i))
		ss = all[locus]['ss'] / denom
		ss_nuc = ss / all[locus]['denom']

		o.write('%s,%s,%s,%s,%.6f,%.6f\n' % (lineage, len(inds) * 2, locus, all[locus]['denom'], ss, ss_nuc))

	o.close()
			

def get_data(args):
	lineage = args.lineage

	d = pd.read_csv(args.file)
	inds = d.ix[d.lineage == lineage, 'sample'].tolist()
	inds = sorted(inds)

	vcf = os.path.join(args.dir, 'variants', 
		'%s.qual_filtered.cov_filtered.vcf.gz' % args.lineage)

	outdir = args.outdir

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	return lineage, inds, vcf, outdir


def main():
	args = get_args()
	lineage, inds, vcf, outdir = get_data(args)
	get_diversity(lineage, inds, vcf, outdir)

if __name__ == "__main__":
	main()
