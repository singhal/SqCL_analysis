import re
import os

nodes = 1
cpu = 2
mem = 8

names = ['AHE1', 'AHE2', 'AHE3', 'AHE4', 'AHE5', 'uce1',  'uce2', 'uce3', 'uce4', 'uce5', 'gene1']
for name in names:
	sh_out = name + '.sh'
        o = open(sh_out, 'w')

        o.write("#PBS -N %s\n" % (name))
        o.write("#PBS -M sosi@umich.edu\n")
        o.write("#PBS -A drabosky_flux\n")
        o.write("#PBS -l qos=flux\n")
        o.write("#PBS -q flux\n")
        o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
        o.write("#PBS -l walltime=10:00:00:00\n")
        o.write("#PBS -j oe\n")
        o.write("#PBS -V\n")
        # o.write("module load bowtie2/2.1.0 trinity/2.3.2\n")
        o.write("\n")

	xml = '/scratch/drabosky_flux/sosi/brazil/phylogeny/div_dating/%s/%s.xml' % (name, name)
	dir = '/scratch/drabosky_flux/sosi/brazil/phylogeny/div_dating/%s/' % name
	o.write('cd %s\n' % dir)
	o.write('java -Xmx8g -jar /home/sosi/bin/beast/lib/beast.jar %s\n' % (xml))
	o.close()
