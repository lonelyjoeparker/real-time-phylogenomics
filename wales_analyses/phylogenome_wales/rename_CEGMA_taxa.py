# rename_CEGMA_taxa.py	
# 
# Description:
# 	Renames taxa in a CEGMA alignment to be human-readable
#	Note: assumes input is aligned already...
#
# Author:
#	@lonelyjoeparker / Joe Parker / RBG Kew / 2016
#
# Usage:
# 
# 	rename_CEGMA_taxa.py  --fasta <input fasta>
#
# Where args:
#
#	<fasta>:			FASTA file to be renamed
#
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
import argparse
import re

# regexps for species - we will be matchin against first two chars on the seq id.
D_mel = '^[0-9]{2}'
A_tha = 'At'
C_ele = 'CE'
H_sap = 'Hs'
S_pom = 'SP'
S_cer = '^Y'

# parse input args
parser = argparse.ArgumentParser(description='Renames taxa in a CEGMA alignment to be human-readable')
parser.add_argument('--fasta',	help='FASTA file to be renamed')
args = parser.parse_args()
aln = []
for seq in SeqIO.parse(args.fasta,'fasta'):
	short = seq.id[:2]
	if re.findall(D_mel,short):
		seq.id = 'D_melanogaster'
	elif short == A_tha:
		seq.id = 'A_thaliana'
	elif short == C_ele:
		seq.id = 'C_elegans'
	elif short == H_sap:
		seq.id = 'H_sapiens'
	elif short == S_pom:
		seq.id = 'S_pombe'
	elif re.findall(S_cer,short):
		seq.id = 'S_cerevisiae'	
	aln.append(seq)
	
SeqIO.write(aln,args.fasta + '.renamed','fasta')