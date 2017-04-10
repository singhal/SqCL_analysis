from rpy2.robjects.packages import importr
import rpy2.robjects as ro

import glob
import multiprocessing as mp
import re
import numpy as np

dir = '/scratch/drabosky_flux/sosi/brazil/phylogeny/gene_trees/'
trees = glob.glob("%s*best*" % dir)
loci = [re.search('.*\/(\S+).best', tree).group(1) for tree in trees]

ape = importr('ape')
phangorn = importr('phangorn')
CPU = 8

compare = []
for ix, locus1 in enumerate(loci):
	print(ix)
	for locus2 in loci[(ix + 1):]:
		compare.append([locus1, locus2])

def tree_dist(comp):
	res = [comp[0], comp[1]]
	t1 = '%s%s.bestTree.tre' % (dir, comp[0])
	t2 = '%s%s.bestTree.tre' % (dir, comp[1])
	tree1 = ape.read_tree(t1)
	tree2 = ape.read_tree(t2)

	setdiff = ro.r('setdiff')
	intersect = ro.r('intersect')
	common = intersect(tree1[2], tree2[2])
	if (len(common) > 4):
		tree1 = ape.drop_tip(tree1, setdiff(tree1[2], common))
		tree2 = ape.drop_tip(tree2, setdiff(tree2[2], common))
		dist = phangorn.RF_dist(tree1, tree2, normalize=True)[0]
	else:
		dist = np.nan
	res.append(dist)
	return res

def run_dist(compare):
	pool = mp.Pool(CPU)
	dist = pool.map(tree_dist, compare)

	return(dist)

dist = run_dist(compare)
out = open('/scratch/drabosky_flux/sosi/brazil/phylogeny/tree_dists.csv', 'w')
out.write('locus1,locus2,distance\n')
for res in dist:
	out.write(','.join([str(x) for x in res]) + '\n')
out.close()
