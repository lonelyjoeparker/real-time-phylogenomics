# interleave_prediction_with_Wickett2014_reference_data.py	
# 
# Description:
# 	Take a file of sequences with Wickett2014 ids  ('FNA.7418') and add to an alignment with Wickett2014 ids  ('FNA.7418') in title
#	Note: assumes input is aligned already...
#
# Author:
#	@lonelyjoeparker / Joe Parker / RBG Kew / 2016
#
# Usage:
# 
# 	interleave_prediction_with_Wickett2014_reference_data.py  --Wickett2014 <input existing Wickett2014 alignment> --fasta <input fasta> --output <output fasta> --name <name>
#
# Where args:
#
#	<fasta>:			FASTA file to be appended
#	<ouput data>:		output file
#	<name>:				Name to give sequences appended
#	<Wickett2014>:			existing Wickett2014 alignment
#
# Sample FASTA input:

from Bio import SeqIO
import re
import argparse

# parse input args
parser = argparse.ArgumentParser(description='Take a file of sequences with Wickett2014 ids  (\'FNA.7418\') and add to an alignment with Wickett2014 ids (\'FNA.7418\') in title')
parser.add_argument('--fasta',	help='FASTA file to be appended')
parser.add_argument('--output',	help='output FASTA file')
parser.add_argument('--name',	help='Name to give sequences appended')
parser.add_argument('--Wickett2014',	help='Existing Wickett2014 alignment')
args = parser.parse_args()

# first find which FNA this is
FNA = re.findall('FNA\.([0-9]{4})',args.Wickett2014)[0]
if not FNA:
	quit()
	
# load in the Wickett2014 FNA alignment
aln = []
sites = 0
for seq in SeqIO.parse(args.Wickett2014,'fasta'):
	aln.append(seq)
	sites = len(seq)
	
# load in the new genes
for seq in SeqIO.parse(args.fasta,'fasta'):
	if re.findall(FNA,seq.id):
		single = seq
		seq.id = args.name + '_' +  FNA
		SeqIO.write(single,args.output+'single-appended','fasta')
		if len(seq) > sites:
			seq.seq = seq.seq[0:sites]
		elif len(seq) < sites:
			while len(seq) < sites:
				seq.seq = seq.seq + '-'
		aln.append(seq)
		
# write out
SeqIO.write(aln,args.output,'fasta')
