#! /usr/bin/perl
# summarise a BLASTN hit table
# e.g. input
#	8b00d8cd-c29d-4525-8bca-72e367caff03_Basecall_Alignment_template	812	82.39	0
#	28d8698f-d88d-4d58-837e-77824f06de7d_Basecall_Alignment_template	158	81.65	4.00E-26
#	5b63a214-d521-4613-a86f-f76351b4b2c7_Basecall_Alignment_template	669	79.97	1.00E-118
#	
# Usage: cat <infile> | sum_blast.pl
# Output: tab delimited:
#	[0]	cumulative sum of alignment lengths
#	[1]	mean of alignment pidents
#	[2]	cumulative sum of alignment identities

# init vars
$cumsum_length = 0;
$cumsum_ident = 0;
$cumsum_pident = 0;
$mean_pident = 0;
$N = 0;

# read STDIN while it's there
while(<>){
	chomp($line=$_);
	#print $line;
	@fields = split("\t",$line);
	if(scalar(@fields)>3){
		$cumsum_length += $fields[1];
		$pident = $fields[2] / 100;
		$cumsum_ident += ($pident * $fields[1]);
		$cumsum_pident += $pident;
		$N++;
	}
}

if($N > 0){$mean_pident = $cumsum_pident / $N;}
print "$cumsum_length\t$cumsum_ident\t$mean_pident\t$N\n";