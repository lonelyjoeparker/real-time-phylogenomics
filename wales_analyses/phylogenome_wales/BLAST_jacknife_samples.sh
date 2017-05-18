# script to generate unholy amounts of pairwise BLASTN data
# uses a very low e-value threshold (1.0)
# so that ROC curves can be calculated post-hoc

# Author
# Dr. Joe Parker, RBG Kew / github @lonelyjoeparker

# !VERSION INFO:
# Diffs 15/11/2016:
#	- ONT thaliana reads set to all_R7R9_thaliana.cleaned.filtered.fasta as this has lambda-phage removed
#	- evalue cutoff is 1.0 for a permissive blast search - filtering / ROC plots are the objective here
#	- A. lyrata ssp. petraea sample TP match will be against *combined* (A.lyrata + A.lyrata ssp. petraea) DBs as neither's that great
#	- compareBlastHits.pl will now run later and output all hits e.g incl those hitting 1 or other DB, not just both
#	- num_threads increased from 4 to 8

# Diffs 2017-05-18:
# - Uncommented all the lines below to enable the full analsys to be re-run *IF* the pathnames below have been edited
# to give real files.
#
# Otherwise will exit with the message below..

# Diffs 2017-05-18 02:
# - Copied this file from massive_BLAST_ROC-onlyCompareBlast.sh to implement a quick analysis of the jacknifed reference data
# - This analysis implements analysis of these jacknife samples via BLAST, orthogonally, e.g.
#
# query petraea, DB lyrata:	 	blastn -db subsample_A.lyra_0000100_draws.fasta -query ../../wales_analyses/phylogenome_wales/inputs/all_R7R9_petraea.fasta  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1|wc
# query petraea, DB thaliana:	blastn -db subsample_A.thal_0000100_draws.fasta -query ../../wales_analyses/phylogenome_wales/inputs/all_R7R9_petraea.fasta  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1|wc
# query thaliana, DB lyrata:	blastn -db subsample_A.lyra_0000100_draws.fasta -query ../../wales_analyses/phylogenome_wales/inputs/all_R7R9_thaliana.cleaned.filtered.fasta  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1|wc
# query thaliana, DB thaliana:	blastn -db subsample_A.thal_0000100_draws.fasta -query ../../wales_analyses/phylogenome_wales/inputs/all_R7R9_thaliana.cleaned.filtered.fasta  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1|wc
#
# /VERSION INFO!

#inputs dir
# /home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/phylogenome_wales/inputs
# lrwxrwxrwx  1 joe joe   69 May 18 20:35 AT2a_S2_L001_all.trimmed.fa -> /media/joe/BiSlDi/WALES-NELUMBOLAB/inputs/AT2a_S2_L001_all.trimmed.fa
# lrwxrwxrwx  1 joe joe   69 May 18 20:34 AL1a_S3_L001_all.trimmed.fa -> /media/joe/BiSlDi/WALES-NELUMBOLAB/inputs/AL1a_S3_L001_all.trimmed.fa
# -rw-r--r--  1 joe joe  66M May 18 15:32 all_R7R9_petraea.fasta
# -rw-r--r--  1 joe joe 251M May 18 15:32 all_R7R9_thaliana.fasta
# -rw-r--r--  1 joe joe 214M May 18 15:32 all_R7R9_thaliana.cleaned.filtered.fasta

# jacknifed databases
# /home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/in-silico-reference-genome-digest
# -rw-rw-r-- 1 joe joe 8.1M May 18 15:26 subsample_A.thal_0100000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe 821K May 18 15:26 subsample_A.thal_0010000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe  83K May 18 15:26 subsample_A.thal_0001000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe 8.4K May 18 15:26 subsample_A.thal_0000100_draws.fasta.nin
# -rw-rw-r-- 1 joe joe  964 May 18 15:26 subsample_A.thal_0000010_draws.fasta.nin
# -rw-rw-r-- 1 joe joe 8.1M May 18 15:26 subsample_A.lyra_0100000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe 821K May 18 15:25 subsample_A.lyra_0010000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe  83K May 18 15:25 subsample_A.lyra_0001000_draws.fasta.nin
# -rw-rw-r-- 1 joe joe 8.4K May 18 15:25 subsample_A.lyra_0000100_draws.fasta.nin
# -rw-rw-r-- 1 joe joe  964 May 18 15:25 subsample_A.lrya_0000010_draws.fasta.nin



# Hardcoded variables: will need to be edited!
# BASEPATH
BASEPATH=/home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/
# DBs - A.thaliana jacknifes
BLAST_DB_THALIANA_10E2=$BASEPATH/in-silico-reference-genome-digest/subsample_A.thal_0000100_draws.fasta
BLAST_DB_THALIANA_10E3=$BASEPATH/in-silico-reference-genome-digest/subsample_A.thal_0001000_draws.fasta
BLAST_DB_THALIANA_10E4=$BASEPATH/in-silico-reference-genome-digest/subsample_A.thal_0010000_draws.fasta
BLAST_DB_THALIANA_10E5=$BASEPATH/in-silico-reference-genome-digest/subsample_A.thal_0100000_draws.fasta
# DBs - A.lyrata jacknifes
BLAST_DB_LYRATA_10E2=$BASEPATH/in-silico-reference-genome-digest/subsample_A.lyra_0000100_draws.fasta
BLAST_DB_LYRATA_10E3=$BASEPATH/in-silico-reference-genome-digest/subsample_A.lyra_0001000_draws.fasta
BLAST_DB_LYRATA_10E4=$BASEPATH/in-silico-reference-genome-digest/subsample_A.lyra_0010000_draws.fasta
BLAST_DB_LYRATA_10E5=$BASEPATH/in-silico-reference-genome-digest/subsample_A.lyra_0100000_draws.fasta
# inputs (queries)
READS_ONT_THALIANA=$BASEPATH/phylogenome_wales/inputs/all_R7R9_thaliana.cleaned.filtered.fasta
READS_ONT_LYRATA=$BASEPATH/phylogenome_wales/inputs/all_R7R9_petraea.fasta
READS_MISEQ_THALIANA=$BASEPATH/phylogenome_wales/inputs/AT2a_S2_L001_all.trimmed.fa
READS_MISEQ_LYRATA=$BASEPATH/phylogenome_wales/inputs/AL1a_S3_L001_all.trimmed.fa
# other
COMPARE_BLAST=$BASEPATH/phylogenome_wales/compareBlastHits.pl
BLAST_OUT=$BASEPATH
BLASTN=/usr/bin/blastn
BLAST_PARAMS=" -num_threads 4 -evalue 1.0 -outfmt \"6 qacc length pident evalue\" -max_target_seqs 1  -max_hsps 1 "

if ! test -e "$COMPARE_BLAST"
then
  echo Did not find one or more variables [e.g. path $COMPARE_BLAST].
  echo You MUST edit the hardcoded variables in this script to use it!
  exit 1
fi


## ORDER of THINGS ##
# For each jacknife from 100:100,000:

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
#	1b ONT thaliana vs DB lyrata (FP)
#	1c ONT lyrata vs DB thaliana (TP)
#	1d ONT lyrata vs DB lyrata (FP)

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
#	2b MiSeq thaliana vs DB lyrata (FP)
#	2c MiSeq lyrata vs DB thaliana (TP)
#	2d MiSeq lyrata vs DB lyrata (FP)

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)


## running it ##

## JACKNIFE: 10E2 ##

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E2 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thaliana-10E2_query-thaliana.out
#	1b ONT tliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E2 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrata-10E2_query-thaliana.out
#	1c ONT lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E2 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thaliana-10E2_query-lyrata.out
#	1d ONT lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E2 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrata-10E2_query-lyrata.out

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E2 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thaliana-10E2_query-thaliana.out
#	2b MiSeq thaliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E2 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrata-10E2_query-thaliana.out
#	2c MiSeq lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E2 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thaliana-10E2_query-lyrata.out
#	2d MiSeq lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E2 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrata-10E2_query-lyrata.out

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
$COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thaliana-10E2_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrata-10E2_query-thaliana.out > $BLAST_OUT/pairwise_10E2_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E2_type-pores_tp-thaliana_fp-lyrata.out.summary
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
$COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrata-10E2_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thaliana-10E2_query-lyrata.out > $BLAST_OUT/pairwise_10E2_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E2_type-pores_tp-lyrata_fp-thaliana.out.summary

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
$COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thaliana-10E2_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrata-10E2_query-thaliana.out > $BLAST_OUT/pairwise_10E2_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E2_type-miseq_tp-thaliana_fp-lyrata.out.summary
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
$COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrata-10E2_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thaliana-10E2_query-lyrata.out > $BLAST_OUT/pairwise_10E2_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E2_type-miseq_tp-lyrata_fp-thaliana.out.summary


## JACKNIFE: 10E3 ##

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E3 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thaliana-10E3_query-thaliana.out
#	1b ONT tliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E3 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrata-10E3_query-thaliana.out
#	1c ONT lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E3 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thaliana-10E3_query-lyrata.out
#	1d ONT lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E3 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrata-10E3_query-lyrata.out

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E3 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thaliana-10E3_query-thaliana.out
#	2b MiSeq thaliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E3 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrata-10E3_query-thaliana.out
#	2c MiSeq lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E3 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thaliana-10E3_query-lyrata.out
#	2d MiSeq lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E3 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrata-10E3_query-lyrata.out

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
$COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thaliana-10E3_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrata-10E3_query-thaliana.out > $BLAST_OUT/pairwise_10E3_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E3_type-pores_tp-thaliana_fp-lyrata.out.summary
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
$COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrata-10E3_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thaliana-10E3_query-lyrata.out > $BLAST_OUT/pairwise_10E3_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E3_type-pores_tp-lyrata_fp-thaliana.out.summary

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
$COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thaliana-10E3_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrata-10E3_query-thaliana.out > $BLAST_OUT/pairwise_10E3_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E3_type-miseq_tp-thaliana_fp-lyrata.out.summary
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
$COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrata-10E3_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thaliana-10E3_query-lyrata.out > $BLAST_OUT/pairwise_10E3_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E3_type-miseq_tp-lyrata_fp-thaliana.out.summary

## JACKNIFE: 10E4 ##

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E4 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thaliana-10E4_query-thaliana.out
#	1b ONT tliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E4 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrata-10E4_query-thaliana.out
#	1c ONT lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E4 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thaliana-10E4_query-lyrata.out
#	1d ONT lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E4 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrata-10E4_query-lyrata.out

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E4 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thaliana-10E4_query-thaliana.out
#	2b MiSeq thaliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E4 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrata-10E4_query-thaliana.out
#	2c MiSeq lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E4 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thaliana-10E4_query-lyrata.out
#	2d MiSeq lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E4 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrata-10E4_query-lyrata.out

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
$COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thaliana-10E4_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrata-10E4_query-thaliana.out > $BLAST_OUT/pairwise_10E4_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E4_type-pores_tp-thaliana_fp-lyrata.out.summary
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
$COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrata-10E4_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thaliana-10E4_query-lyrata.out > $BLAST_OUT/pairwise_10E4_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E4_type-pores_tp-lyrata_fp-thaliana.out.summary

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
$COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thaliana-10E4_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrata-10E4_query-thaliana.out > $BLAST_OUT/pairwise_10E4_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E4_type-miseq_tp-thaliana_fp-lyrata.out.summary
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
$COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrata-10E4_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thaliana-10E4_query-lyrata.out > $BLAST_OUT/pairwise_10E4_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E4_type-miseq_tp-lyrata_fp-thaliana.out.summary

## JACKNIFE: 10E5 ##

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E5 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thaliana-10E5_query-thaliana.out
#	1b ONT tliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E5 -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrata-10E5_query-thaliana.out
#	1c ONT lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E5 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thaliana-10E5_query-lyrata.out
#	1d ONT lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E5 -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrata-10E5_query-lyrata.out

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E5 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thaliana-10E5_query-thaliana.out
#	2b MiSeq thaliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E5 -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrata-10E5_query-thaliana.out
#	2c MiSeq lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA_10E5 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thaliana-10E5_query-lyrata.out
#	2d MiSeq lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA_10E5 -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrata-10E5_query-lyrata.out

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
$COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thaliana-10E5_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrata-10E5_query-thaliana.out > $BLAST_OUT/pairwise_10E5_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E5_type-pores_tp-thaliana_fp-lyrata.out.summary
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
$COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrata-10E5_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thaliana-10E5_query-lyrata.out > $BLAST_OUT/pairwise_10E5_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E5_type-pores_tp-lyrata_fp-thaliana.out.summary

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
$COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thaliana-10E5_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrata-10E5_query-thaliana.out > $BLAST_OUT/pairwise_10E5_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_10E5_type-miseq_tp-thaliana_fp-lyrata.out.summary
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
$COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrata-10E5_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thaliana-10E5_query-lyrata.out > $BLAST_OUT/pairwise_10E5_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_10E5_type-miseq_tp-lyrata_fp-thaliana.out.summary
