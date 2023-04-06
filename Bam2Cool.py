import pysam
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Convert Hi-C bam into mcool') #change

parser.add_argument('-I', '--BAM', required=True, help='Type input Hi-C bam file')
parser.add_argument('-r', '--RES', required=True, help='Type input resolution (5kb, 10kb, 25kb, 100kb)')

args = parser.parse_args()

bam = args.BAM
res = args.RES

if res not in ['5kb', '10kb', '25kb', '100kb']:
	print('invalid resolution!')
	sys.exit()

def bam2cool(bam, res):
	ToolDir = sys.path[0]
	contact_file = f'{bam}_contact.txt'
	chrom_size = f'{ToolDir}/bin/hg19.chr.size.txt'
	write_contact =open(contact_file, 'w')
	samfile = pysam.AlignmentFile(bam, "rb")
	for read in samfile:
		if read.rnext != -1:
			chrom1= read.reference_name
			chrom2= read.next_reference_name
			pos1= read.pos
			pos2= read.pnext
			write_contact.write(f'{chrom1}\t{pos1}\t{chrom2}\t{pos2}\n')
	samfile.close()
	write_contact.close()

	sorted_contact_file = f'{contact_file}.sort'
	csort = ['cooler', 'csort', '-c1 1', '-p1 2', '-c2 3', '-p2 4', '-o', sorted_contact_file, contact_file, chrom_size]
	os.system(' '.join(csort))
	os.remove(contact_file)
	bin_file = f'{ToolDir}/bin/hg19.{res}_bin.bed'
	cool_file = f'{bam}_{res}.cool'
	cload = ['cooler', 'cload', 'pairix', bin_file, sorted_contact_file, cool_file]
	os.system(' '.join(cload))
	os.remove(sorted_contact_file)
	os.remove(sorted_contact_file + '.px2')

	zoomify = ['cooler', 'zoomify', '-o', cool_file.replace('cool', 'mcool'), '-r'] 
	if res == '5kb':
		zoomify.append('5000,10000,25000,50000,100000,200000,500000,1000000,2000000,5000000')
	elif res == '10kb':
		zoomify.append('10000,20000,50000,100000,200000,500000,1000000,2000000,5000000')
	elif res == '25kb':
		zoomify.append('25000,50000,100000,200000,500000,1000000,2000000,5000000')
	elif res == '100kb':
		zoomify.append('100000,200000,500000,1000000,2000000,5000000')
	else:
		print('invalid resolution!')
		sys.exit()
	zoomify.append('--balance --balance-args "--mad-max 5 --min-nnz 10 --min-count 0 --ignore-diags 2"')
	zoomify.append(cool_file)
	os.system(' '.join(zoomify))
	os.remove(cool_file)

bam2cool(bam, res)




