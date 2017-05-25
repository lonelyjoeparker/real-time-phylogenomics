# To implement this:
# 		Load the list of TAIR genes to use as keys
# 		For each BLAST hit:
# 			parse subject gene
# 			parse query up to _Basecall as read_UID, compare to current_read_UID + subject gene
# 			parse length (ignore if length < X) for e.g. X = 500?
# 			if we've not seen this gene subject/read combo, increment the read count for that gene
# 			if we have already seen this read **AND** subject gene is the same then SNAP has (probably) just failed to stitch the exons together; do nothing
# 			if we have seen this read BUT the gene subject is different then a single read probably contains multiple genes; increment each gene subject
# 		Then we can print for each TAIR gene:
# 			name of TAIR
# 			count of reads (int)
# 			list of read IDs (csv)
# 
# This will let us lookup read IDs by time (eventually)

# Pseudocode:
#
# Load TAIR
# Foreach TAIR:
#	parse TID from input e.g.:
#		input='ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010'
#		parse '=|;' for tokens
#		assert (token 1 == token n)
#		add to list
#
# Initialise current_GID, TABLE
# Load SNAP
# Foreach SNAP:
#	parse GID, length
#	ignore length < x
#	if !TID{read} increment TABLE
#	else 

import argparse

# parse input args

parser = argparse.ArgumentParser(description='Match BLASTN (read SNAP predicted genes vs. TAIR10 genes) results to  read name/time info (from poretools.fasta + poretools.index via parse_poretools.index.py) and output to stdout.')
parser.add_argument('--times',help='tab-delimited nonempty flat data file (times from poretools.fasta + poretools.index via parse_poretools.index.py)')
parser.add_argument('--genes',help='nonempty flat data file (genes from TAIR10 annotation)')
parser.add_argument('--blast',help='nonempty flat data file (blastn results)')
parser.add_argument("--verbose", help="increase output verbosity (0=none, 1=everything)[0]")
args = parser.parse_args()

genes = []
times = {}

# Load times
with open(args.times,'r') as f:
	for line in f:
		times_data = line.rstrip().split("\t")
		if len(times_data) == 3:
			times[times_data[0]] = (times_data[1],times_data[2])


# Load TAIR
tair = open(args.genes,'r') # e.g. /Users/joeparker/Downloads/wales_nelumbolab_output/TAIR10_GFF3_names.txt
# Foreach TAIR:
for t in tair:
	tokens = t.split('=')
	gene = tokens[1].split(';')[0] 
	if args.verbose:
		print gene # print if verbose
	genes.append(gene)
	
#	parse TID from input e.g.:
#		input='ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010'
#		parse '=|;' for tokens
#		assert (token 1 == token n)
#		add to list

"""
# bit of testing collections
test_gene = raw_input('enter a gene ID ')
if test_gene in genes:
	print 'found it'
else:
	print 'NO'
"""	

if args.verbose:
	print 'TAIR10 genes list is', len(genes), 'long'

# Initialise current_GID, TABLE
# Load SNAP
blasts = open(args.blast,'r') # e.g 'export.blast.results'
occurence = {}
reads_by_gene = {}
blast_results = {}
genes_matched = 0
genes_missed = 0
blast_total = 0
# Foreach SNAP:
for b in blasts:
#	parse GID, length
	subject_gene = b[0:9]
	read_UID = b[12:48]
	read_alignment = subject_gene + read_UID
	blast_results[read_alignment] = b
	if subject_gene in genes:
		genes_matched += 1
		if (subject_gene in occurence) and (subject_gene in reads_by_gene):
			occurence[subject_gene] = occurence[subject_gene] + 1
			reads_by_gene[subject_gene].add(read_UID)
		else:
			occurence[subject_gene] = 1
			reads_by_gene[subject_gene] = set([read_UID])
		#print subject_gene
	else:
		genes_missed += 1
	blast_total +=1

if args.verbose:
	print "matched: ", genes_matched
	if blast_total != (genes_matched + genes_missed):
		print "ERROR: BLAST key checksum"
	
	print 'total BLAST lines processed:',blast_total
#	ignore length < x
#	if !TID{read} increment TABLE
#	else 

	print 'genes hit with reads =', len(occurence.keys())

# now, finally, we can print out each TAIR10 gene,
#	with a list of reads that hit that gene
#	and a time in POSIX for each of those reads
# print header first
print "TAIR10_gene\tread_UID\ttimes\tread_info\tread\tchannel\tgene_check\tread_SNAP.pred\tlength\tpident\tevalue\tgaps\tnum_gap_open"
for k in sorted(occurence.keys()):
	if args.verbose:
		print k, occurence[k], '[_]', reads_by_gene[k], '[_]',	len(reads_by_gene[k])
	# recall that reads_by_gene is a list of read UIDs e.g. set(['2f32fe5f-7f51-4455-b5c7-71a8dfcf8880', '3df2b437-f55b-4d11-9187-f50728cc10c5']) etc
	for read in reads_by_gene[k]:
		read_alignment = k + read
		if read in times and read_alignment in blast_results:
			read_channel_info = times[read][1].split("_")	# may as well output read and channel #s as clean ints while we're here
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (k,read,times[read][0],times[read][1],read_channel_info[0][2:],read_channel_info[1][4:],blast_results[read_alignment])


"""
# OK this basically works
Now to make the subsampling work all we
have to do is also read in a list/dict of 
times-reads, then for a given dt simply
test reads against that time for existence

and hence get curve, hey presto!

ALSO:

If we do a bit of r on the output from this:
> TAIR = read.table("SNAP_TAIR.tdf",sep="\t",header=T)
> str(TAIR)
'data.frame':	9644 obs. of  4 variables:
 $ locus      : Factor w/ 9644 levels "AT1G01040","AT1G01060",..: 5936 3960 5558 9027 7195 3508 3833 1423 2318 6516 ...
 $ hits_total : int  2 2 1 1 4 1 4 2 1 2 ...
 $ set        : Factor w/ 8734 levels "set([000590ee-35d9-4921-b6ef-295c1a6c3275])",..: 4519 7702 8045 511 3687 6539 1239 190 1716 5531 ...
 $ hits_unique: int  1 1 1 1 2 1 2 1 1 1 ...
> hist(TAIR$hits_total,ylim=c(0,30))
> hist(TAIR$hits_unique,ylim=c(0,30))
> summary(TAIR$hits_total)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   1.000    1.000    2.000    2.981    3.000 1615.000 
> summary(TAIR$hits_unique)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.000   1.000   1.000   1.914   1.000 880.000 

etc...

"""