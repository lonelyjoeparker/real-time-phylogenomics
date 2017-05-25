# fasta_to_phylip.py	
# 
# Description:
# 	Converts a fasta file to phylip format
#	Note: assumes input is aligned already...
#
# Author:
#	@lonelyjoeparker / Joe Parker / RBG Kew / 2016
#
# Usage:
# 
# 	fasta_to_phylip.py  --fasta <input fasta> --output <output phylip> [--trim <characters to retain from names]
#
# Where args:
#
#	<fasta>:			FASTA file to be appended
#	<ouput data>:		output file
#	<trim>:				Integer number of characters [1:10] to retain from start of names
#
# Sample FASTA input:
'''
>7291732___KOG0002
maahksfrikqklakklkqnrsvpqwvrlrtgntirynakrrhwrrtklkl
>At3g02190___KOG0002
-------mikkklgkkmrqnrpipnwirlrtdnkirynakrrhwrrtklgf
>SPCC663.04___KOG0002
mpshksfrtkqklakaarqnrpipqwirlrtgntvhynmkrrhwrrtklni
>YJL189w___KOG0002
maaqksfrikqkmakakkqnrplpqwirlrtnntirynakrrnwrrtkmni
>CE06883___KOG0002
msalkksfikrklakkqkqnrpmpqwvrmktgntmkynakrrhwrrtklkl
>Hs4506647___KOG0002
msshktfrikrflakkqkqnrpipqwirmktgnkirynskrrhwrrtklgl
'''
# Sample PHYLIP output:
'''
 6 51
CE06883___ MSALKKSFIK RKLAKKQKQN RPMPQWVRMK TGNTMKYNAK RRHWRRTKLKL
SPCC663.04 MPSHKSFRTK QKLAKAARQN RPIPQWIRLR TGNTVHYNMK RRHWRRTKLNI
Hs4506647_ MSSHKTFRIK RFLAKKQKQN RPIPQWIRMK TGNKIRYNSK RRHWRRTKLGL
At3g02190_ -------MIK KKLGKKMRQN RPIPNWIRLR TDNKIRYNAK RRHWRRTKLGF
YJL189w___ MAAQKSFRIK QKMAKAKKQN RPLPQWIRLR TNNTIRYNAK RRNWRRTKMNI
7291732___ MAAHKSFRIK QKLAKKLKQN RSVPQWVRLR TGNTIRYNAK RRHWRRTKLKL

'''
# Output / details:
#	Will read input KOG data, and fasta file, remove reads that appear in filter,
#	and write amended .fasta to output
from Bio import SeqIO
import argparse

# parse input args
parser = argparse.ArgumentParser(description='Converts a fasta file to phylip format. Note: assumes input is aligned already...')
parser.add_argument('--fasta',	help='FASTA file to be appended')
parser.add_argument('--output',	help='output PHYLIP file')
parser.add_argument('--trim',	help='characters to retain from name (start at left)', type=int, choices=range(1,11))
args = parser.parse_args()

# init global variables
aln = []

# read alignment

for seq in SeqIO.parse(args.fasta,'fasta'):
	if args.trim > 1:
		seq.id = seq.id[args.trim:]
	aln.append(seq)

SeqIO.write(aln,args.output,'phylip')