import re
import glob

print('locus,tc')
files = glob.glob('/scratch/drabosky_flux/sosi/brazil/phylogeny/tree_certainty/*')
for file in files:
	f = open(file,'r')
	for l in f:
		if re.search('Relative tree certainty for this tree:', l):
			tc = re.search(': ([0-9|\.]+)', l).group(1)
			locname = re.search('info.(\S+)', file).group(1)
			print('%s,%s' % (locname, tc))
			break
	f.close()
