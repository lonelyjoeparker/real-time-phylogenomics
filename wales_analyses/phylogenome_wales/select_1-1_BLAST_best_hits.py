'''
~/Downloads/phylogenome_wales/select_1-1_BLAST_best_hits.py

Author:
	Joe Parker, RBG Kew, 2016

Input assumptions:

	Two blast-able fasta files containing a reference dataset of genes/proteins,
	and a test/experimental one.
		Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa 
		those 6628 SNAP-predicted proteins with length > 100 aa etc in BLASTP vs. TAIR10
		Downloads/phylogenome_wales/CEGMA_reference_data/A.thaliana.cegma.aa.fa those 
		sequences from the '248alignments' CEGMA cores set with 'At' in the sequence ID 
		(via python)

	Make blastdbs for each:
	makeblastdb -dbtype prot -in SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa -out 'SNAP.A.thal'
	makeblastdb -dbtype prot -in A.thaliana.cegma.aa.fa -out 'cegma.A.thal'
	
	Now blast each way with top hits only:
	blastp -db SNAP.A.thal -query ../CEGMA_reference_data/A.thaliana.cegma.aa.fa \
	-num_threads 8 -max_hsps 1 -num_alignments 1 \
	-outfmt "6 qacc sacc length pident evalue" > BLASTP.cegma-on-snapDB.out

Method:
	Uses two (reciprocal) BLAST searches to identify orthologues in a reference and a 
	focal (experimental) data set, then writes those test orthologues out as fasta. 

Output:

- A fasta file with only 1:1 orthologues (hopefully) between test and reference
  dataset; matches from the test dataset including reference IDs
- A lookup file with the mappings
'''

from Bio import SeqIO
import argparse

# parse input args
parser = argparse.ArgumentParser(description=
'''
Uses two (reciprocal) BLAST searches to identify orthologues in a reference and a 
focal (experimental) data set, then writes those test orthologues out as fasta. 
Output will include:

- A fasta file with only 1:1 orthologues (hopefully) between test and reference
  dataset; matches from the test dataset including reference IDs
- A lookup file with the mappings
'''
)
parser.add_argument('--fasta',	 help='fasta file containing genes/proteins from focal/test dataset. orthologues identified by 1:1 blast will be written out')
parser.add_argument('--blast_1', help='nonempty flat tab-delimited data file (blast results for focal/test dataset as queries)')
parser.add_argument('--blast_2', help='nonempty flat tab-delimited data file (blast results for reference dataset as queries)')
parser.add_argument("--verbose", help="increase output verbosity (0=none, 1=everything)[0]")
args = parser.parse_args()

# global args
aln  = [] 		# list: the sequences we'll retain from FASTA as othologues and print out
matched = []	# list: those queries from BLAST 1 focal/test data which have best-hits in BLAST 2
blast_1 = {} 	# dict: the query-subject mappings from BLAST 1 (assumed to be test on reference DB)
blast_2 = {}	# dict: the query-subject mappings from BLAST 2 (assumed to be reference on test DB)
blast_scores_1 = {}	# dict:  the results from BLAST 2 (assumed to be reference on test DB)
blast_scores_2 = {}	# dict:  the results from BLAST 2 (assumed to be reference on test DB)

# load BLAST 1 library
with open(args.blast_1,'r') as f:
	for line in f:
		data = line.rstrip().split("\t")
		blast_1[data[0]] = data[1]
		blast_scores_1[data[0]] = line

# load BLAST 1 library
with open(args.blast_2,'r') as f:
	for line in f:
		data = line.rstrip().split("\t")
		blast_2[data[0]] = data[1]
		blast_scores_2[data[0]] = line

# walk through looking for reciprocal matches		
for putative in blast_1:
	reciprocal = blast_2.get(blast_1[putative])
	if reciprocal == putative:
		# we have a match
		matched.append(putative)
		if args.verbose:
			print('MATCH',putative,blast_1[putative])
	else:
		# no match
		if args.verbose:
			print('drop',putative,blast_1[putative])

# print out the matches
for name in matched:
	print(blast_scores_1[name])
	
# get the sequences
for seq in SeqIO.parse(args.fasta,'fasta'):
	for orthologue in  matched:
		if orthologue in seq.id:		
			new_seq_id = blast_1[orthologue]
			print seq.id, new_seq_id
			seq.id = new_seq_id + '|' + seq.id
			aln.append(seq)
			
# write out 
SeqIO.write(aln,args.fasta + '.reciprocal.fa','fasta')