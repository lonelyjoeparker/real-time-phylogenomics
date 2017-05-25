#####
# select_fasta_based_on_BLAST.py
# 
# author
#	Joe Parker, RBG Kew, 2016
#
# description
#	Uses a list of sequences (e.g. that meet some arbitrary
#	threshold in a BLAST analysis) to filter a .fasta file
# 
# usage: 
#	python select_fasta_based_on_BLAST.py --fasta <fasta input> --good <list of good sequences>
#
# input:
#	--fasta	a file in .fasta format
#	--good a tab-delimited list of BLAST results. 'good'  sequences IDs are expected to be in 
#		the 3rd column (index [2])
#
# output
#	<fasta>.selected - the input file with only sequences in <good> printed out
#####

# testing data
'''
good_kogs = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.CEGMA.hits.out.KOG-matches-only'
fasta = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.fa'
'''

from Bio import SeqIO
import argparse

# parse arguments
parser = argparse.ArgumentParser(description=
'''
# select_fasta_based_on_BLAST.py
# 
# author
#	Joe Parker, RBG Kew, 2016
#
# description
#	Uses a list of sequences (e.g. that meet some arbitrary
#	threshold in a BLAST analysis) to filter a .fasta file
# 
# usage: 
#	python select_fasta_based_on_BLAST.py --fasta <fasta input> --good <list of good sequences>
#
# input:
#	--fasta	a file in .fasta format
#	--good a tab-delimited list of BLAST results. 'good'  sequences IDs are expected to be in 
#		the 1st column (index [0])
#
# output
#	<fasta>.selected - the input file with only sequences in <good> printed out. also appends the subject hits (column 2, index [1]) if present.
'''
)
parser.add_argument('--fasta',	 help='fasta file containing genes/proteins from focal/test dataset.')
parser.add_argument('--good', help='nonempty flat tab-delimited data file (blast results deemed OK to keep)')
args=parser.parse_args()

# globals
good = {}	# hash holding IDs for hits that were good, and their corresponding subject IDs if possible
aln = []	# list to hold the output alignment itself

# read in 'good' list and parse to hash
with open(args.good,'r') as f:
	for line in f:
		data = line.rstrip().split("\t")
		good[data[0]] = data[1]
		
# walk through fasta file, adding only sequences with matches in 'good' to the output buffer 'aln'
for seq in SeqIO.parse(args.fasta,'fasta'):
	for hit in good:
		if hit in seq.id:
			print(seq.id, hit, good[hit])
			seq.id = seq.id + good[hit]
			aln.append(seq)
			
# write filtered ouput as fasta
SeqIO.write(aln,args.fasta + '.selected.' + str(len(aln)),'fasta')