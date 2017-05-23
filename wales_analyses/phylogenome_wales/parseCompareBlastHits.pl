#! /usr/bin/perl

# parseCompareBlastHits.pl
#
# Author @lonelyjoeparker / Joe Parker, RBG Kew, 2017
#
# A perl script to parse the *out.summary output files from a compareBlastHits.pl analysis
# of two comparative blast runs into a single line of STDOUT (fields tab-delimited). 
# Typically this would be concatenated into a file.
#
# Usage:
# 	perl parseCompareBlastHits.pl <input.example.out.summary>
#
# Input:
#	*out.summary from compareBlastHits.pl (analysis of pairwise comparative BLAST)
#
# Output:
# 	A single line of tab-delimited data with the following information:
#		Tab number	Field
#		[0]	input filename
#		[1]	input sample (assumed)
#		[2]	pores | miseq (which sequencing platform)
#		[3]	tp-database (assumed true-positive sample, e.g. should be same as [1]
#		[4]	fp-database (assumed incorrect/false-positive)
#		[5]	hits_2-way (number of reads hitting both BLAST DBs)
#		[6]	hits_1-way_T (number of reads hitting only correct / TP DB)
#		[7]	hits_1-way_F (number of reads hitting only incorrect / FP DB)
#		[8]	mean_bias_len (mean value of 'length bias' ID stat amongst all reads)
#		[9]	mean_bias_ident (mean value of 'num. identities bias' ID stat amongst all reads)

# initialise output variables with some sensible defaults that allow parse errors to be noticed
# NB not much explicit exception handling in this code

$input_filename		= 'NULL';
$input_sample		= 'NULL';
$platform			= 'NULL';
$tp_database		= 'NULL';	
$fp_database		= 'NULL';
$hits_2_way 		= "NaN";
$hits_1_way_T		= "NaN";
$hits_1_way_F		= "NaN";
$mean_bias_length	= "NaN";
$mean_bias_nident	= "NaN";

# parse the input
$input_filename = $ARGV[0];
$unique_hits_flag = 0;
open(IN,$input_filename);
while(<IN>){
	# tidy up
	chomp($line = $_);
	@line_fields = split(/\t/,$line);
	
	# rule for sample
	if($line =~ /_query\-([A-Za-z0-9]+)/){
		$input_sample = substr($&,7);
	}
	
	# rule for platform
	if($line =~ 'pores|miseq'){
		$platform = $&;
	}
	
	# rule to work out which DB it is
	if($line =~ /_db\-([A-Za-z0-9\-\_\.]+)_query/){
		$DB_guess = substr($&,4,-6);
		if($line =~ /BLAST\ input\ 1/){
			$tp_database = $DB_guess;
		}else{
			if($line =~ /BLAST\ input\ 2/){
				$fp_database = $DB_guess;
			}
		}
	}
	
	# rule for 2-way hits
	if($line =~ 'both'){
		if($line =~ /[0-9]+/){
			$hits_2_way = $&;
		}
	}
	
	# rule for unique hits
	if($unique_hits_flag>0){
		# 	if flag is set:
		if($line =~ /\t[0-9]+/){
			$nums = substr($&,1);
			if($unique_hits_flag == 1){
				$hits_1_way_T = $nums;
				$unique_hits_flag = 2;
			}else{
				if($unique_hits_flag == 2){
					$hits_1_way_F = $nums;
				}
			}
		}
	}else{
		# 	if flag is not currently set:
		if($line =~ 'Unique hits'){
			# set flag; unique hits on the next two lines
			$unique_hits_flag = 1;
		}
	}
	
	# rule for mean length bias
	if($line =~ "mean length bias"){
		$mean_bias_length = $line_fields[2];
	}

	# rule for mean nident bias
	if($line =~ "mean identities bias"){
		$mean_bias_nident = $line_fields[2];
	}
}

# print the output

print join("\t",
$input_filename,	
$input_sample,	
$platform,		
$tp_database,	
$fp_database,	
$hits_2_way, 	
$hits_1_way_T,	
$hits_1_way_F,	
$mean_bias_length,	
$mean_bias_nident
)."\n";
