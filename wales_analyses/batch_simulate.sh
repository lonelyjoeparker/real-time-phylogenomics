#! /usr/bin/bash

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
# Diffs 2017-05-22:
# - Copied / forked from BLAST_repeat_reference_genome.sh
# - Now takes good/middling genomes (A. thaliana TAIR10 and A.lyrata 1.0), then chops them into smaller contig pieces using Chop-finished-reference-genome-into-contigs-runnable.py
# - Then builds blast DBs, then BLASTs
# /VERSION INFO!

# directories to simulate into
#drwxrwxr-x  2 joe joe 4.0K May 22 15:14 simulate_rubber
#drwxrwxr-x  2 joe joe 4.0K May 22 15:14 simulate_ash
#drwxrwxr-x  2 joe joe 4.0K May 22 15:14 simulate_monkey
#drwxrwxr-x  2 joe joe 4.0K May 22 15:14 simulate_pepper
#drwxrwxr-x  2 joe joe 4.0K May 22 15:14 simulate_nelumbo
#drwxrwxr-x  2 joe joe 4.0K May 22 15:13 simulate_potato

# Hardcoded variables: will need to be edited!
# BASEPATH
BASEPATH=/home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/
# inputs (queries)
READS_ONT_THALIANA=$BASEPATH/phylogenome_wales/inputs/all_R7R9_thaliana.cleaned.filtered.fasta
READS_ONT_LYRATA=$BASEPATH/phylogenome_wales/inputs/all_R7R9_petraea.fasta
READS_MISEQ_THALIANA=$BASEPATH/phylogenome_wales/inputs/AT2a_S2_L001_all.trimmed.fa
READS_MISEQ_LYRATA=$BASEPATH/phylogenome_wales/inputs/AL1a_S3_L001_all.trimmed.fa
# other
COMPARE_BLAST=$BASEPATH/phylogenome_wales/compareBlastHits.pl
BLASTN=/usr/bin/blastn
BLAST_PARAMS=" -num_threads 8 -evalue 1.0 -outfmt \"6 qacc length pident evalue\" -max_target_seqs 1  -max_hsps 1 "

# attempt to check the required, hardcoded variables are present
if ! test -e "$COMPARE_BLAST"
then
  echo Did not find one or more variables [e.g. path $COMPARE_BLAST].
  echo You MUST edit the hardcoded variables in this script to use it!
  exit 1
fi

# delete old outputs
rm -r simulate_*/*fasta
rm -r simulate_*/*fasta.n*

function simulate_and_analyse {
  N50=$1    # target N50 in bp
  simdir=$2 # target output directory
  echo $N50, $simdir

  ## set up ##

  # simulate poor quality genome assemblies
  python in-silico-reference-genome-digest/Chop-finished-reference-genome-into-contigs-runnable.py $N50 $simdir
  output_lyrata_sim=`ls $simdir/*lyra*.fasta`
  output_thalia_sim=`ls $simdir/*thal*.fasta`
  # build a BLASTN database
  makeblastdb -dbtype 'nucl' -in $output_lyrata_sim
  makeblastdb -dbtype 'nucl' -in $output_thalia_sim

  # try and create a BLAST output directory (may fail if it already exists, e.g. from a previous analysis. we're not auto deleting at the start as these babies take AGES to run)
  BLAST_OUT=$BASEPATH/$2/blasts
  mkdir $BLAST_OUT

  ## analyse the poo quality simulation ##

  # 1. Blast ONT reads:
  #	1a ONT thaliana vs DB thaliana (TP)
  $BLASTN -db $output_thalia_sim -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thalianaSIM_query-thaliana.out
  #	1b ONT tliana vs DB lyrata (FP)
  $BLASTN -db $output_lyrata_sim -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrataSIM_query-thaliana.out
  #	1c ONT lyrata vs DB thaliana (TP)
  $BLASTN -db $output_thalia_sim -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thalianaSIM_query-lyrata.out
  #	1d ONT lyrata vs DB lyrata (FP)
  $BLASTN -db $output_lyrata_sim -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrataSIM_query-lyrata.out

  # 2. Blast MiSeq reads:
  #	2a MiSeq thaliana vs DB thaliana (TP)
  $BLASTN -db $output_thalia_sim -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thalianaSIM_query-thaliana.out
  #	2b MiSeq thaliana vs DB lyrata (FP)
  $BLASTN -db $output_lyrata_sim -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrataSIM_query-thaliana.out
  #	2c MiSeq lyrata vs DB thaliana (TP)
  $BLASTN -db $output_thalia_sim -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thalianaSIM_query-lyrata.out
  #	2d MiSeq lyrata vs DB lyrata (FP)
  $BLASTN -db $output_lyrata_sim -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrataSIM_query-lyrata.out

  # 3. Pairwise comparisons - ONT
  #	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
  $COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thalianaSIM_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrataSIM_query-thaliana.out > $BLAST_OUT/pairwise_SIM_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_SIM_type-pores_tp-thaliana_fp-lyrata.out.summary
  #	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
  $COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrataSIM_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thalianaSIM_query-lyrata.out > $BLAST_OUT/pairwise_SIM_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_SIM_type-pores_tp-lyrata_fp-thaliana.out.summary

  # 4. Pairwise comparisons - MiSeq
  #	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
  $COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thalianaSIM_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrataSIM_query-thaliana.out > $BLAST_OUT/pairwise_SIM_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_SIM_type-miseq_tp-thaliana_fp-lyrata.out.summary
  #	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
  $COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrataSIM_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thalianaSIM_query-lyrata.out > $BLAST_OUT/pairwise_SIM_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_SIM_type-miseq_tp-lyrata_fp-thaliana.out.summary

}


simulate_and_analyse 3400000 simulate_nelumbo
simulate_and_analyse 1400000 simulate_potato
simulate_and_analyse 2470000 simulate_pepper
simulate_and_analyse  104000 simulate_ash
simulate_and_analyse 1280000 simulate_rubber
