import re
import glob
import os
import multiprocessing as mp
import subprocess

files = glob.glob('/scratch/drabosky_flux/sosi/brazil/phylogeny/alignments/*gb')
outdir = '/scratch/drabosky_flux/sosi/brazil/phylogeny/informativeness/'
mafft = 'mafft'
CPU = 24
num_miss = 40

def get_seq(locus_file):
        f = open(locus_file, 'r')
	locus = re.sub('^.*/', '', locus_file)
	locus = re.sub('\.fasta.*$', '', locus)

        seq = {}
        id = ''
        for l in f:
                if re.search('>', l):
                        id = re.search('>(\S+)', l.rstrip()).group(1)
                        if re.search('^_R_', id):
                                id = re.sub('^_R_', '', id)
                        seq[id] = ''
                else:
                        seq[id] += l.rstrip()
        f.close()

	if 'Gallus_gallus' in seq:
		del seq["Gallus_gallus"]

	if len(seq) >= num_miss:
		out = os.path.join(outdir, '%s.fasta' % locus)
	        o = open(out, 'w')
		for id, s in seq.items():
			o.write('>%s\n%s\n' % (id, re.sub('\s+', '', s)))
		o.close()

		return out
	else:
		return False


def align(params):
        file, mafft = params

        aln_out = file.replace('.fasta', '.fasta.aln')
        proc = subprocess.call("%s --maxiterate 1000 --globalpair "
                               "--adjustdirection --quiet %s > %s" %
                               (mafft, file, aln_out), shell=True)

        os.remove(file)
        return aln_out


def run_alignments(files):
        params = zip(files, [mafft] * len(files))

	pool = mp.Pool(CPU)
	alns = pool.map(align, params)
        
        return alns


# alns = []
# for file in files:
# 	aln = get_seq(file)
#	if aln:
#		alns.append(aln)
#alns2 = run_alignments(alns)

alns2 = glob.glob('/scratch/drabosky_flux/sosi/brazil/phylogeny/informativeness/*aln')
for aln in alns2:
	out = re.sub('.fasta.aln', '.nex', aln)
	if not os.path.isfile(out):
		subprocess.call("python /home/sosi/SqCL_revisions/seq.file.converter.py -i %s -inf FASTA -outf NEXUS > %s" % (aln, out), shell=True)
