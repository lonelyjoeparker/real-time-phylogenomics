# python bodges
SCRIPT_INTERLEAVE=/Users/joeparker/Downloads/phylogenome_wales/interleave_prediction_with_CEGMA_reference_data.py

# RAxML binary
BINARY_RAXML=/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS

# fasta file containing the SNAP gene predictions (selected genes)
SNAP=/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa.reciprocal.fa

# trimal alignment threshold (required % taxa with data at each site)
TRIM_FACTOR=0.7

#loops, yay
for f in `cat singles`; 
do 
	# name the files
	SINGLE=$f.aln.renamed.added.snap-appended.fullTaxonSet.KOG.fasingle-appended
	PROFILE=$f.aln.renamed.added
	ALIGNMENT=$f.profile.aligned.fa
	TRIMMED=$f.profile.aligned.trimal.$TRIM_FACTOR.phy
	
	#first interleave the sequence we need
	python $SCRIPT_INTERLEAVE --fasta $SNAP --output $PROFILE.snap-appended.fullTaxonSet.KOG.fa --name A_thal_AA_____ --cegma $PROFILE
	#first align main set
	muscle -in $PROFILE -out tmp 
	# then profile-align new sequence to main set
	muscle -profile -in1 tmp -in2 $SINGLE -out $ALIGNMENT 
	# convert to phylip for RAxML and trim at the same time
	trimal -in $ALIGNMENT -out $TRIMMED -phylip_m10 -gt $TRIM_FACTOR
	# RAxML, dayhoff matrix, pthreads=8
	$BINARY_RAXML -T 8 -s $TRIMMED -m PROTCATDAYHOFF -n $f  
	mv RAxML_bestTree.$f $f.tre; 
	rm RAxML* tmp; 
done 
