import glob
import re

dir = '/Volumes/heloderma3/brazil/mitogenomes/'
files = glob.glob('/Volumes/heloderma3/brazil/mitogenomes/*mitogenome.fa')

out = dir + 'brazil_mitogenomes.fasta'
o = open(out, 'w')

for file in files:
	f = open(file, 'r')
	id = re.search('([^/]+)_mitog', file).group(1)
	#out1 = re.sub('fasta', 'fa', file)
	#o1 = open(out1, 'w')

	head = f.next()
	seq = ''
	for l in f:
		seq += l.rstrip()
	seq = seq.upper()
	seq = re.sub('^X+', '', seq)
	seq = re.sub('^N+', '', seq)
	seq = re.sub('X+$', '', seq)
	seq = re.sub('N+$', '', seq)

	o.write('>%s\n%s\n' % (id, seq))
	#o1.write('>%s\n%s\n' % (id, seq))

	print '%s %s' % (id, seq.count('N'))
	
	#o1.close()
	f.close()
o.close()
