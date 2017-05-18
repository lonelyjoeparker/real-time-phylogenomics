#! /usr/bin/perl
#
# compareBlastHits.pl
#
# Synopsis:
# 	Compare queries' blast alignments between two blast searches
#
# 	NB assumes only one hit/HSP per query (-max_hits 1 -max_hsps 1)
#	NB assumes blast output format 6 (tab delimited) with fields {name, length, identity %, e-value}
#	NB e.g. `blastn <db> <query> -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1`
#
# Usage:
#	perl compareBlastHits.pl <blast_search_1_output> <blast_search_2_output>
#
# Output:
#	Tab-separated data, one per query with:
#	query	length_search_1	evalue_search_1	length_search_2	evalue_search_2
#
# Changelog
#	Diffs 16/01/2017:
#	 - number of idents now CORRECTLY counted as [length * (% idents / 100)]
#	Diffs 15/11/2016:
#	 - output data lines for *all* hits including those in *only* one DB not both


open(SEARCH1,$ARGV[0]);
open(SEARCH2,$ARGV[1]);

%search_1_hits = {};
%search_2_hits = {};

#read in the hits

while(<SEARCH1>){
	chomp($line=$_);
	@pieces = split(/\t/,$line);
	$search_1_hits{$pieces[0]} = $line;
	#print "token [$pieces[0]]\n";
}

while(<SEARCH2>){
	chomp($line=$_);
	@pieces = split(/\t/,$line);
	$search_2_hits{$pieces[0]} = $line;
}

close(SEARCH1);
close(SEARCH2);

# initialise stats for output
$hits_1_not_2 = 0;	# number of query sequences present in BLAST 1 alignments but not BLAST 2
$hits_2_not_1 = 0;	# number of query sequences present in BLAST 2 alignments but not BLAST 1
$hits_both = 0;		# number of query sequences present both BLAST 1 alignments and BLAST 2. their scores will be compared.
$cumulative_length_1 = 0;	# cumulative length of alignments in BLAST 1 that are also present in BLAST 2
$cumulative_length_2 = 0;	# cumulative length of alignments in BLAST 2 that are also present in BLAST 2
$cumulative_idents_1 = 0;	# cumulative count of identities in alignments in BLAST 1 that are also present in BLAST 2
$cumulative_idents_2 = 0;	# cumulative count of identities in alignments in BLAST 2 that are also present in BLAST 2
$cumulative_pident_1 = 0;	# cumulative proportion of identities of alignments in BLAST 1 that are also present in BLAST 2
$cumulative_pident_2 = 0;	# cumulative proportion of identities of alignments in BLAST 2 that are also present in BLAST 2
$cumulative_evalue_1 = 0;	# cumulative evalue of alignments in BLAST 1 that are also present in BLAST 2
$cumulative_evalue_2 = 0;	# cumulative evalue of alignments in BLAST 2 that are also present in BLAST 2

# print header
print "key\tBLAST_1length\tBLAST_1pident\tBLAST_1evalue\tBLAST_2length\tBLAST_2pident\tBLAST_2evalue\n";

# iterate through both sets of BLAST results, looking for those present in both BLASTs
foreach $key(keys(%search_1_hits)){
	if(exists($search_2_hits{$key})){
		#	Diffs 16/01/2017:
		#	 - number of idents now CORRECTLY counted as [length * (% idents / 100)]
		#
		# reminder: blast data in tab-delimited, fields {qacc length pident evalue}
		@data_1 = split(/\t/,$search_1_hits{$key});
		@data_2 = split(/\t/,$search_2_hits{$key});
		$length_1 = $data_1[1];	# alignment length BLAST 1
		$length_2 = $data_2[1];	# alignment length BLAST 2
		$pident_1 = $data_1[2];	# alignment % identity BLAST 1
		$pident_2 = $data_2[2];	# alignment % identity BLAST 2
		$evalue_1 = $data_1[3];	# alignment e-value BLAST 1
		$evalue_2 = $data_2[3];	# alignment e-value BLAST 2
		$idents_1 = $length_1 * ($pident_1 / 100.0); # number of identities BLAST 1
		$idents_2 = $length_2 * ($pident_2/ 100.0); # number of identities BLAST 2
		print "$key\t$data_1[1]\t$data_1[2]\t$data_1[3]\t$data_2[1]\t$data_2[2]\t$data_2[3]\n";
		$hits_both++;
		$cumulative_length_1 += $length_1;
		$cumulative_length_2 += $length_2;
		$cumulative_idents_1 += $idents_1;
		$cumulative_idents_2 += $idents_2;
		$cumulative_pident_1 += $pident_1;
		$cumulative_pident_2 += $pident_2;
		$cumulative_evalue_1 += $evalue_1;
		$cumulative_evalue_2 += $evalue_2;
	}else{
		$hits_1_not_2++;
		# Changelog
		#	Diffs 16/01/2017:
		#	 - number of idents now CORRECTLY counted as [length * (% idents / 100)]
		#	Diffs 15/11/2016:
		#	 - output data lines for *all* hits including those in *only* one DB not both
		#	 - null/NA values for missing hits will now be added with the following:
		#		length: -1
		#		pident:	-1
		#		evalue: 999
		#	  	...although this is non-standard from a statistical point of view, we're doing
		#	  it in this case since these values are obviously biologically wrong (so can be 
		#	  coerced to NAs in R later if needed) but will still allow the ROC stats to be
		#	  calculated.
		#		The resultant plots should be correct, except that in the case of x-axis 
		#	  (cutoff) plots the xlim parameter may need to be set to zoom in on the correct
		#	  plot region. However curves for TP and FP rates will be exaggerated (see for e.g.
		#	  /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/wales_nelumbolab_output/ROC_effect_of_missing_observations_as_extreme_values.pdf
		#	  which shows this for dummy data with NA present or coded as extreme values) because
		#	  proportions classified at cutoff points will be in error (too high) ... so the 
		#	  best thing to do is probably to repeat ROC plots recoding as well to check effects.
		@data_1 = split(/\t/,$search_1_hits{$key});
		$length_1 = $data_1[1];	# alignment length BLAST 1
		$pident_1 = $data_1[2];	# alignment % identity BLAST 1
		$evalue_1 = $data_1[3];	# alignment e-value BLAST 1
		$idents_1 = $length_1 * ($pident_1 / 100.0); # number of identities BLAST 1
		print "$key\t$data_1[1]\t$data_1[2]\t$data_1[3]\t-1\t-1\t999\n";
		$cumulative_length_1 += $length_1;
		$cumulative_idents_1 += $idents_1;
		$cumulative_pident_1 += $pident_1;
		$cumulative_evalue_1 += $evalue_1;
	}
}

foreach $key(keys(%search_2_hits)){
	unless(exists($search_1_hits{$key})){
			$hits_2_not_1++;
	# Changelog
	#	Diffs 16/01/2017:
	#	 - number of idents now CORRECTLY counted as [length * (% idents / 100)]
	#	Diffs 15/11/2016:
	#	 - output data lines for *all* hits including those in *only* one DB not both
	#	 - for more see above
	@data_2 = split(/\t/,$search_2_hits{$key});
	$length_2 = $data_2[1];	# alignment length BLAST 2
	$pident_2 = $data_2[2];	# alignment % identity BLAST 2
	$evalue_2 = $data_2[3];	# alignment e-value BLAST 2
	$idents_2 = $length_2 * ($pident_2 / 100.0); # number of identities BLAST 2
	print "$key\t-1\t-1\t999\t$data_2[1]\t$data_2[2]\t$data_2[3]\n";
	$cumulative_length_2 += $length_2;
	$cumulative_idents_2 += $idents_2;
	$cumulative_pident_2 += $pident_2;
	$cumulative_evalue_2 += $evalue_2;
	}
}

# calculate summary stats for those reads which had hits for both BLASTs
# first stats on cumulative sums
$cumulative_length_bias = sprintf ("%+d", $cumulative_length_1 - $cumulative_length_2);	# diff sum(lengths_1) - sum(lengths_2)
$cumulative_idents_bias = sprintf ("%+f", $cumulative_idents_1 - $cumulative_idents_2);	# diff sum(idents_1) - sum(idents_2)
$cumulative_pident_bias = sprintf ("%+f", $cumulative_pident_1 - $cumulative_pident_2);	# diff sum(pidents_1) - sum(pidents_2)
$cumulative_evalue_bias = sprintf ("%+e", $cumulative_evalue_2 - $cumulative_evalue_1);	# diff sum(evalues_1) - sum(evalues_2) NOTE that sign is reversed as LOW evalues are significant!!

# now calculate means
$mean_length_1 = ($cumulative_length_1 / $hits_both);	
$mean_idents_1 = ($cumulative_idents_1 / $hits_both);	
$mean_pident_1 = ($cumulative_pident_1 / $hits_both);	
$mean_evalue_1 = ($cumulative_evalue_1 / $hits_both);	
$mean_length_2 = ($cumulative_length_2 / $hits_both);	
$mean_idents_2 = ($cumulative_idents_2 / $hits_both);	
$mean_pident_2 = ($cumulative_pident_2 / $hits_both);	
$mean_evalue_2 = ($cumulative_evalue_2 / $hits_both);	

# finally stats on means
$mean_length_bias = sprintf ("%+f", $mean_length_1 - $mean_length_2);	# diff sum(lengths_1) - sum(lengths_2)
$mean_idents_bias = sprintf ("%+f", $mean_idents_1 - $mean_idents_2);	# diff sum(idents_1) - sum(idents_2)
$mean_pident_bias = sprintf ("%+f", $mean_pident_1 - $mean_pident_2);	# diff sum(pidents_1) - sum(pidents_2)
$mean_evalue_bias = sprintf ("%+e", $mean_evalue_2 - $mean_evalue_1);	# diff sum(evalues_1) - sum(evalues_2) NOTE that sign is reversed as LOW evalues are significant!!


warn "BLAST hits comparison\n(assumes IDENTICAL queries to each BLAST database, IDENTICAL parameters, and max_hits=1, max_hsps=1)\n";
warn "BLAST input 1:\t$ARGV[0]\n";
warn "BLAST input 2:\t$ARGV[1]\n";
warn "Total hits report:\nNumber of query sequences with hits present in both:\t$hits_both\n";
warn "Stats\n";
warn "\tCumulative length bias ($ARGV[0]: $cumulative_length_1 - $ARGV[1]: $cumulative_length_2)\t$cumulative_length_bias\t$cumulative_length_1\t$cumulative_length_2\n";
warn "\tCumulative identities bias ($ARGV[0]: $cumulative_idents_1 - $ARGV[1]: $cumulative_idents_2)\t$cumulative_idents_bias\t$cumulative_idents_1\t$cumulative_idents_2\n";
warn "\tCumulative % identities bias ($ARGV[0]: $cumulative_pident_1 - $ARGV[1]: $cumulative_pident_2)\t$cumulative_pident_bias\t$cumulative_pident_1\t$cumulative_pident_2\n";
warn "\tCumulative evalues bias ($ARGV[0]: $cumulative_evalue_1 - $ARGV[1]: $cumulative_evalue_2)\t$cumulative_evalue_bias\t$cumulative_evalue_1\t$cumulative_evalue_2\n";
warn "\tmean length bias ($ARGV[0]: $mean_length_1 - $ARGV[1]: $mean_length_2)\t$mean_length_bias\t$mean_length_1\t$mean_length_2\n";
warn "\tmean identities bias ($ARGV[0]: $mean_idents_1 - $ARGV[1]: $mean_idents_2)\t$mean_idents_bias\t$mean_idents_1\t$mean_idents_2\n";
warn "\tmean % identities bias ($ARGV[0]: $mean_pident_1 - $ARGV[1]: $mean_pident_2)\t$mean_pident_bias\t$mean_pident_1\t$mean_pident_2\n";
warn "\tmean evalues bias ($ARGV[0]: $mean_evalue_1 - $ARGV[1]: $mean_evalue_2)\t$mean_evalue_bias\t$mean_evalue_1\t$mean_evalue_2\n";
warn "Unique hits\n$ARGV[0]\t$hits_1_not_2\n$ARGV[1]\t$hits_2_not_1\n";
