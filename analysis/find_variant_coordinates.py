import re
import pandas as pd
import gzip

sp = 'Colobosaura_modesta'
c = '/Volumes/sosi/brazil/metadata/%s_to_probe_coords.csv' % sp
v = '/Volumes/sosi/brazil/variants/%s.qual_filtered.cov_filtered.vcf.gz' % sp
out = '/Volumes/sosi/brazil/metadata/%s_variants.csv' % sp

coords = {}
c = pd.read_csv(c)
for c, s, e in zip(c.locus, c.start, c.stop):
	coords[c] = int((int(s) + int(e)) / 2.0)

f = gzip.open(v, 'r')
o = open(out, 'w')
o.write('locus,absgap,relgap,maf\n')
for l in f:
	l = l.decode('ascii')
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		if d[4] != '.':
			genos = d[9:len(d)]
			genos = [re.search('(\S/\S)', geno).group(1) for geno in genos]
			g = []
			for geno in genos:
				geno = re.split('/', geno)
				if geno[0] != '.':
					g += geno
			if len(g) > 0:
				af = g.count('1') / float(len(g))
				if af > 0.5:
					af = 1 - af
					if af > 0:
						o.write('%s,%s,%s,%s\n' % (d[0], d[1], int(d[1]) - coords[d[0]], af))
o.close()
f.close()
