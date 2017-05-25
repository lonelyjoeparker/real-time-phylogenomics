#### filter fasta by cutoff .py ####
#
# input a gene predictions (fasta) file
# and a set of BLAST validation data
# (gene predictions blasted against DB)
# and a list of KOG <-> gene ID mappings
# from the reference.
#
# output a fasta file
# containing only those predicted genes
# which have hits to the validation DB
# above a threshold
# appending KOG IDs where available 
# for secondary checks.
#
# Author: Joe Parker, RBG Kew, 2016
# joe.parker@kew.org
#
# Usage: python filter_fasta_by_cutoff.ply
# (no-args - inputs hardcoded as below)
#
####################################

from Bio import SeqIO


## hardcoded inputs

# holds the BLAST validation data, column 'filter' contains [TRUE|FALSE]
cutoff_file_aa = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/all_R7R9_thaliana.genetimes.filtered.cutoff.aa.rdata'
cutoff_file_nt = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/all_R7R9_thaliana.genetimes.filtered.cutoff.nt.rdata'

# the predicted genes with TAIR10 gene names
predicted_aa = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.aa.filtered.fa'
predicted_nt = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.dna.filtered.fa'

# file (tab-delimited) from JGI containing KOG IDs for TAIR10 A thaliana genes
kog_file = '/Users/joeparker/Downloads/phylogenome_wales/JGI_data/mart_export_A_thaliana_TAIR10.tdf'

# output files for fasta NT ans AA
output_aa = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa'
output_nt = '/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.fa'


## global collections

# hashes to hold the filters' output
cutoff_aa = {}
cutoff_nt = {}

# hash to hold KOG IDs
kogs = {}

# lists to hold the selected sequences
aln_aa = []
aln_nt = []


## process the data

# add the aa cutoff data as a hash
# keys will be the unique gene IDs
# values will be the TAIR10 assignment (putative)
# only add filter == TRUE

with open(cutoff_file_aa,'r') as f:
	for line in f:
		data = line.rstrip().split("\t")
		if data[16] == 'TRUE':
#			print(data[0],data[7],data[16])
			cutoff_aa[data[7]] = data[0]
	
# add the nt cutoff data as a hash
# keys will be the unique gene IDs
# values will be the TAIR10 assignment (putative)
# only add filter == TRUE

with open(cutoff_file_nt,'r') as f:
	for line in f:
		data = line.rstrip().split("\t")
		if data[len(data)-1] == 'TRUE':
			cutoff_nt[data[7]] = data[0]

# get the KOG IDs for these putative genes 

with open(kog_file) as file:
	for line in file:
		data = line.split("\t")
		# add to global KOG IDs hash
#		print(data[0],data[len(data)-1][:7])
		kogs[data[0]] = data[len(data)-1][:7]
	

# walk through AA sequences
for seq in SeqIO.parse(predicted_aa,'fasta'):
	if seq.id in cutoff_aa:
		this_kog = 'no_kog_value'
		if cutoff_aa[seq.id] in kogs:
			this_kog = kogs[cutoff_aa[seq.id]]
		seq.id = seq.id + '|' + cutoff_aa[seq.id] + '|' + this_kog
#		print(seq.id, seq.name, cutoff_aa[seq.id],this_kog,len(seq))
		aln_aa.append(seq)

# walk through NT sequences
for seq in SeqIO.parse(predicted_nt,'fasta'):
	if seq.id in cutoff_nt:
		this_kog = 'no_kog_value'
		if cutoff_nt[seq.id] in kogs:
			this_kog = kogs[cutoff_nt[seq.id]]
		seq.id = seq.id + '|' + cutoff_nt[seq.id] + '|' + this_kog
		aln_nt.append(seq)

# write output		
SeqIO.write(aln_aa,output_aa,'fasta')
SeqIO.write(aln_nt,output_nt,'fasta')
