# add_KOG_IDs_to_JGI_data.py	
# 
# Description:
# 	Given a tab-delimited input file of 
#	gene ID / KOG ID data from JGI and 
#	a fasta file, appends KOG IDs to fasta
#	names
#
# Author:
#	@lonelyjoeparker / Joe Parker / RBG Kew / 2016
#
# Usage:
# 
# 	add_KOG_IDs_to_JGI_data.py --kog <input KOG data> --fasta <input fasta> --output <output file>
#
# Where args:
#
#	<fasta>:			FASTA file to be appended
#	<input KOG data>:	JGI output (tab-delimited) with gene IDs in <fasta> and KOG IDs
#	<ouput data>:		output file
#
# Sample KOG input:
'''
AT3G11500	AT3G11500.1	Athaliana	19662802_peptide	KOG1780
AT5G64990	AT5G64990.1	Athaliana	19669958_peptide	KOG0094
AT3G61620	AT3G61620.2	Athaliana	19660986_peptide	KOG1068
'''
# Sample FASTA input:
'''
>ATCG00905|ATCG00905.1|Athaliana|19638011_peptide
ACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACC
CGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGG
GTCAAAAAGCCAAAATAA
>AT2G20790|AT2G20790.3|Athaliana|19638095_peptide
ATGAAGGATTACGTCAAGCTTTGTAGGAGATCTGATTGTGGACCAGCTGTTGGTGAAGAT
TGTATGTCTTTGGCTTCCTTTTGA
>AT2G25210|AT2G25210.1|Athaliana|19638358_peptide
ATGATTAAGAAGAAATTGGGGAAGAAAATGAGGCAAAACAGGCCTATTCCTCACTGGATT
CGTCTTCGTACCGACAACACCATCAGGTACAATGCTAAGCGTAGGCATTGGCGCAGAACC
AAGCTTGGATTCTAA
'''
# Sample FASTA output:
'''
>ATCG00905|ATCG00905.1|Athaliana|19638011_peptide|KOG1750 ATCG00905|ATCG00905.1|Athaliana|19638011_peptide
ACTATCACCCCCAAAAAACCAAACTCTGCTTTACGTAAAGTTGCCAGAGTACGATTAACC
CGAGGAACCCTAGATGCTGTCGGAGTAAAGGATCGTCAACAAGGGCGTTCTAAATATGGG
GTCAAAAAGCCAAAATAA
>AT2G20790|AT2G20790.3|Athaliana|19638095_peptide|KOG0937 AT2G20790|AT2G20790.3|Athaliana|19638095_peptide
ATGAAGGATTACGTCAAGCTTTGTAGGAGATCTGATTGTGGACCAGCTGTTGGTGAAGAT
TGTATGTCTTTGGCTTCCTTTTGA
>AT2G25210|AT2G25210.1|Athaliana|19638358_peptide|KOG0002 AT2G25210|AT2G25210.1|Athaliana|19638358_peptide
ATGATTAAGAAGAAATTGGGGAAGAAAATGAGGCAAAACAGGCCTATTCCTCACTGGATT
CGTCTTCGTACCGACAACACCATCAGGTACAATGCTAAGCGTAGGCATTGGCGCAGAACC
AAGCTTGGATTCTAA
'''
# Output / details:
#	Will read input KOG data, and fasta file, remove reads that appear in filter,
#	and write amended .fasta to output
from Bio import SeqIO
import argparse

# parse input args
parser = argparse.ArgumentParser(description='Will read .fasta input, and tab-delimited BioMart data (containing gene IDs and KOG IDs), then write .fasta  with KOG IDs appended to output.')
parser.add_argument('--kog',	help='JGI output (tab-delimited) with gene IDs in <fasta> and KOG IDs')
parser.add_argument('--fasta',	help='FASTA file to be appended')
parser.add_argument('--output',	help='output .fasta file')
args = parser.parse_args()

# init global variables
kogs = {}
new_seqs = []

# read in the .tdf with KOG IDs and gene IDs
with open(args.kog) as file:
	for line in file:
		data = line.split("\t")
		# add to global KOG IDs hash
		kogs[data[0]] = data[len(data)-1][:7]

for seq in SeqIO.parse(args.fasta,'fasta'):
	seq_data = seq.id.split("|")[0]
	if seq_data in kogs:
		# this sequence gene ID has an entry from the gene / KOG ID hash
		seq.id = seq.id + '|' + kogs[seq_data]
		new_seqs.append(seq)

SeqIO.write(new_seqs,args.output,'fasta')