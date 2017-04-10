import re
import gzip
import glob
import random
import numpy as np

out = open('/Users/sonal/publications/SqCL/genetic_diversity/genetic_diversity_loci.csv', 'w')

indir = '/Users/sonal/publications/SqCL/variants/'
vcfs = glob.glob(indir + "*gz")

header = ['ind']
header = header + ['all_denom', 'all_div', 
					'uce_denom', 'uce_div',
					'uce_div_sub_AHE_95l', 'uce_div_sub_AHE_95u', 'uce_div_sub_AHE_median',
					'uce_div_sub_gene_95l', 'uce_div_sub_gene_95u', 'uce_div_sub_gene_median',
					'AHE_denom', 'AHE_div',
					'AHE_div_sub_gene_95l', 'AHE_div_sub_gene_95u', 'AHE_div_sub_gene_median',
					'gene_denom', 'gene_div']
out.write(','.join(header) + '\n')

for vcf in vcfs:
	o = gzip.open(vcf, 'r')
	div = {}
	for l in o:
		l = l.decode('utf-8')
		if re.search('#CHROM', l):
			d = re.split('\t', l.rstrip())
			inds = d[9:]
			for ind in inds:
				div[ind] = {'AHE': [], 
				            'uce': [], 
				            'gene': [] 
				            }
		elif not re.search('#', l):
			d = re.split('\t', l.rstrip())
			genos = d[9:]

			loc_type = re.split('-', d[0])[0]
			for ind, geno in zip(inds, genos):
				allele1 = re.search('^(\S)', geno).group(1)
				allele2 = re.search('/(\S)', geno).group(1)

				if allele1 in ['0', '1'] and allele2 in ['0', '1']:
					div[ind][loc_type].append('%s%s' % (allele1, allele2))
	o.close()

	for ind in div:
		res = dict([(x, 0) for x in header])
		res['ind'] = ind

		res['all_denom'] = sum([len(div[ind][x]) for x in div[ind]])
		res['all_div'] = sum([div[ind][x].count('01') for x in div[ind]]) / float(res['all_denom'])

		for type1 in ['AHE', 'uce', 'gene']:
			res['%s_denom' % type1] =  len(div[ind][type1])
			res['%s_div' % type1] = div[ind][type1].count('01') / float(res['%s_denom' % type1])

		for type1, type2 in [('uce', 'AHE'), ('uce', 'gene'), ('AHE', 'gene')]:
			stem = '%s_div_sub_%s' % (type1, type2)
			subsample = res['%s_denom' % type2]
			if subsample > len(div[ind][type1]):
				subsample = len(div[ind][type1]) - 1
			divs = []
			for i in range(100):
				sites = random.sample(div[ind][type1],  subsample)
				tmpdiv = sites.count('01') / float(len(sites))
				divs.append(tmpdiv)
			res['%s_95l' % stem] = np.percentile(divs, 2.5)
			res['%s_95u' % stem] = np.percentile(divs, 97.5)
			res['%s_median' % stem] = np.median(divs)


		vals = [str(res[x]) for x in header]
		out.write(','.join(vals) + '\n')
	print(vcf)
out.close()

