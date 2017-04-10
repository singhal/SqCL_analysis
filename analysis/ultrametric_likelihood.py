import multiprocessing as mp
import sys
import glob
import re

# so it can find the local version of dendropy
sys.path = ['/home/sosi/.local/lib/python2.7/site-packages'] + sys.path
import dendropy
from dendropy.interop import paup

cpu = 8
indir = '/scratch/drabosky_flux/sosi/brazil/phylogeny/alignments/'

def get_tree(aln):
	locname = re.sub('^.*/', '', aln)
	locname = re.sub('\.aln.*', '', locname)
	data = dendropy.DnaCharacterMatrix.get(path=aln, schema="phylip")
	print(aln)
	try:
		# estimate starting tree
		# (estimating ML tree takes WAY too long)
		treefile = '/scratch/drabosky_flux/sosi/brazil/phylogeny/gene_trees/%s.tre' % locname
		tree = dendropy.Tree.get(path=treefile, schema="newick",
			preserve_underscores=True, taxon_namespace=data.taxon_namespace)

		# estimate ultrametric tree
		ult_tree = paup.estimate_ultrametric_tree(data, topology_tree=tree,
				paup_path="/home/sosi/bin/paup4a152_centos64")

		# estimate likelihood of data given bl tree
		est_tree, est_model = paup.estimate_model(data, tree, num_states=2,
				unequal_base_freqs=True, gamma_rates=False, prop_invar=True,
				paup_path="/home/sosi/bin/paup4a152_centos64")
		lnl1 = est_model['log_likelihood']

		# estimate likelihood of data given ultrametric tree
		est_ult_tree, est_ult_model = paup.estimate_model(data, ult_tree,
			num_states=2, unequal_base_freqs=True, gamma_rates=False, prop_invar=True,
			paup_path="/home/sosi/bin/paup4a152_centos64")
		lnl2 = est_ult_model['log_likelihood']

		res = [locname, lnl1, lnl2, len(tree.taxon_namespace), len(data.sequences()[0])]
	except:
		res = [locname, 'NA', 'NA', 'NA', 'NA']

	return res

def run_trees(alns):       
	pool = mp.Pool(cpu)
	res = pool.map(get_tree, alns)
        
        return res

alns = glob.glob(indir + '*phy')
# alns = alns[0:5]
res = run_trees(alns)
res = []
for aln in alns:
	print(aln)
	tmp = get_tree(aln)
	res.append(tmp)
out = open('/scratch/drabosky_flux/sosi/brazil/phylogeny/ultrametric_filtering.csv', 'w')
out.write('locus,lnl_nonultra,lnl_ultra,ntips,loc_length\n')
for d in res:
	out.write(','.join([str(x) for x in d]) + '\n')
out.close()
