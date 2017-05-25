# interleave_prediction_with_CEGMA_reference_data.py	
# 
# Description:
# 	Take a file of sequences with KOG ids and add to an alignment with KOG IDs in title
#	Note: assumes input is aligned already...
#
# Author:
#	@lonelyjoeparker / Joe Parker / RBG Kew / 2016
#
# Usage:
# 
# 	interleave_prediction_with_CEGMA_reference_data.py  --cegma <input existing CEGMA alignment> --fasta <input fasta> --output <output fasta> --name <name>
#
# Where args:
#
#	<fasta>:			FASTA file to be appended
#	<ouput data>:		output file
#	<name>:				Name to give sequences appended
#	<cegma>:			existing CEGMA alignment
#
# Sample FASTA input:
'''

Reminder, taxlabels:

7291732	Drosophila melanogaster
At3g02190	Arabidopsis thaliana
CE06883	C. elegans
Hs4506647	Homo sapiens
SPCC663	Saccharomyces pombe
YJL189w	Saccharomyces cerevisiae
'''

from Bio import SeqIO
import re
import argparse

# parse input args
parser = argparse.ArgumentParser(description='Take a file of sequences with KOG ids and add to an alignment with KOG IDs in title')
parser.add_argument('--fasta',	help='FASTA file to be appended')
parser.add_argument('--output',	help='output PHYLIP file')
parser.add_argument('--name',	help='Name to give sequences appended')
parser.add_argument('--cegma',	help='Existing cegma alignment')
args = parser.parse_args()

# first find which KOG this is
KOG = re.findall('KOG[0-9]{4}',args.cegma)[0]
if not KOG:
	quit()
	
# load in the cegma KOG alignment
aln = []
sites = 0
for seq in SeqIO.parse(args.cegma,'fasta'):
	aln.append(seq)
	sites = len(seq)
	
# load in the new genes
for seq in SeqIO.parse(args.fasta,'fasta'):
	if re.findall(KOG,seq.id):
		seq.id = args.name + '_' +  KOG
		SeqIO.write(seq,args.output+'single-appended','fasta')
		if len(seq) > sites:
			seq.seq = seq.seq[0:sites]
		elif len(seq) < sites:
			while len(seq) < sites:
				seq.seq = seq.seq + '-'
		aln.append(seq)
		
# write out
SeqIO.write(aln,args.output,'fasta')
