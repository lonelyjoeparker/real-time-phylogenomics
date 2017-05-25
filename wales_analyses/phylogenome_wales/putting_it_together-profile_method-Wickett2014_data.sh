# python bodges
SCRIPT_FILTER=/Users/joeparker/Downloads/phylogenome_wales/Wickett-eA-2014_data/FNA/filter_alignment_by_taxon_name.py
SCRIPT_INTERLEAVE=/Users/joeparker/Downloads/phylogenome_wales/interleave_prediction_with_Wickett2014_reference_data.py

# RAxML binary
BINARY_RAXML=/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS

# list of selected genes to process
SELECTED_GENES=/Users/joeparker/Downloads/phylogenome_wales/Wickett-eA-2014_data/FNA/selected_alignments

# fasta file containing the SNAP gene predictions (selected genes)
SNAP=/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.fa.selected.56

# check to see if we really want to do the subset (above) or in fact all genes (if so change SELECTED_GENES)
# implies this command invoked with (e.g. `bash ../../putting_it_together-profile_method-Wickett2014_data.sh all`)
if [ $1 == 'all' ]; then
	echo '+=+= ANALYSING ALL GENES =+=+'
	# list of all BLASTN-validated genes to process
	SELECTED_GENES=/Users/joeparker/Downloads/phylogenome_wales/Wickett-eA-2014_data/FNA/all_207_alignments
	
	# fasta file containing the SNAP gene predictions (all genes)
	SNAP=/Users/joeparker/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.fa
fi;
  
# taxon sets for filtering script
FILTER_10='thaliana|Populus|Vitis|Zea|Oryza|Juniperus|Gingko|Equisetum|Bryum|parkeae'
FILTER_7='thaliana|Vitis|Zea|Oryza|Juniperus|Equisetum'
FILTER_6='thaliana|Vitis|Zea|Oryza|Juniperus'


#loops, for trees, yay
for f in `cat $SELECTED_GENES`; 
do 

# name the files for this gene
DOWNSAMPLED=$f.downsampled
SINGLE=$f.interleavedsingle-appended
ALIGNMENT=$f.profile.aligned.fa
TRIM_FACTOR=0.5
TRIMMED=$f.profile.aligned.trimal.$TRIM_FACTOR.phy

# filter main alignment for selected sequences
python $SCRIPT_FILTER --input $f --keep 'thaliana|Vitis|Zea|Oryza|Juniperus|Equisetum' > $DOWNSAMPLED
# interleave / split this gene for SNAP
python $SCRIPT_INTERLEAVE --fasta $SNAP --output $f.interleaved --name A_tha_SNAP____ --Wickett2014 $DOWNSAMPLED
# then profile-align new sequence to main set
muscle -profile -in1 $DOWNSAMPLED -in2 $SINGLE -out $ALIGNMENT
# convert to phylip for RAxML and trim at the same time
trimal -in $ALIGNMENT -out $TRIMMED -phylip_m10 -gt $TRIM_FACTOR
# RAxML, dayhoff matrix, pthreads=8
#/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS -T 8 -s $TRIMMED -m PROTCATDAYHOFF -n $f  
$BINARY_RAXML -T 8 -s $TRIMMED -m GTRCAT -n $f -N 10
mv RAxML_bestTree.$f $f.tre; 
rm RAxML* $f.interleaved; 
done 

## afterwards: work out how many taxa in each alignment to concat treefile:
#	for i in *profile.aligned.fa; do echo $i;grep '>' $i|wc; done>taxon-counts
#	edit taxon-counts to select only some > selected-trees
#	for i in `cat selected-trees`; do cat $i >> FNA.selected.trees.10; done
#	for i in thaliana Populus Vitis Zea Oryza Juniperus Gingko Equisetum Bryum parkeae; do echo $i; grep '>' FNA.*profile*fa| grep $i |wc;done
#
# **SEE**: ~/Downloads/phylogenome_wales/Wickett-eA-2014_data/FNA/attempt-11-taxa-selected-56-genes/select-trees.sh 
#
###
