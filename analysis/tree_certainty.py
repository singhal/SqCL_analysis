import multiprocessing as mp
import sys
import glob
import re
import subprocess

trees = glob.glob("/scratch/drabosky_flux/sosi/brazil/phylogeny/gene_trees/*boot*")
cpu = 8

def run_certainty(tree):
	locname = re.sub('.*\/', '', tree)
	locname = re.sub('\.boot.*', '', locname)

	subprocess.call('raxmlHPC-PTHREADS -L MR -z %s -m GTRCAT -n %s' % (tree, locname), shell=True)

	return locname

def run_trees(trees):       
        pool = mp.Pool(cpu)
        res = pool.map(run_certainty, trees)

run_trees(trees)
