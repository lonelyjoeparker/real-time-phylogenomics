# script to generate unholy amounts of pairwise BLASTN data
# uses a very low e-value threshold (1.0)
# so that ROC curves can be calculated post-hoc

# VERSION:
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

# Hardcoded variables: will need to be edited!
BLAST_DB_THALIANA=/mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_thaliana/A_thaliana
BLAST_DB_LYRATA=/mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_both_ssp_merged/A_lyrata_ssp
READS_ONT_THALIANA=/mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/all_R7R9_thaliana.cleaned.filtered.fasta
READS_ONT_LYRATA=/mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/all_R7R9_petraea.fasta
READS_MISEQ_THALIANA=/mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/AT2a_S2_L001_all.trimmed.fa
READS_MISEQ_LYRATA=/mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/AL1a_S3_L001_all.trimmed.fa
COMPARE_BLAST=/mnt/HDD_2/joe/WALES-NELUMBOLAB/compareBlastHits.pl
BLAST_OUT=/mnt/HDD_2/joe/WALES-NELUMBOLAB/BLAST_ROC
BLASTN=/home/joe/Downloads/ncbi-blast-2.4.0+/bin/blastn
BLAST_PARAMS=" -num_threads 4 -evalue 1.0 -outfmt \"6 qacc length pident evalue\" -max_target_seqs 1  -max_hsps 1 "

if ! test -e "$COMPARE_BLAST"
then
  echo Did not find one or more variables [e.g. path $COMPARE_BLAST].
  echo You MUST edit the hardcoded variables in this script to use it!
  exit 1
fi

ls -l $BLAST_DB_THALIANA	# this will fail as the blast DB path has to be the name roots without *nin *nsq *nhr suffix
ls -l $BLAST_DB_LYRATA	# this will fail as the blast DB path has to be the name roots without *nin *nsq *nhr suffix
ls -l $READS_ONT_THALIANA
ls -l $READS_ONT_LYRATA
ls -l $READS_MISEQ_THALIANA
ls -l $READS_MISEQ_LYRATA
ls -l $COMPARE_BLAST
echo $BLAST_PARAMS

## ORDER of THINGS ##

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

# 1. Blast ONT reads:
#	1a ONT thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1a_type-pores_db-thaliana_query-thaliana.out
#	1b ONT tliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA -query $READS_ONT_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1b_type-pores_db-lyrata_query-thaliana.out
#	1c ONT lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1c_type-pores_db-thaliana_query-lyrata.out
#	1d ONT lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA -query $READS_ONT_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/1d_type-pores_db-lyrata_query-lyrata.out

# 2. Blast MiSeq reads:
#	2a MiSeq thaliana vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2a_type-miseq_db-thaliana_query-thaliana.out
#	2b MiSeq thaliana vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA -query $READS_MISEQ_THALIANA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2b_type-miseq_db-lyrata_query-thaliana.out
#	2c MiSeq lyrata vs DB thaliana (TP)
$BLASTN -db $BLAST_DB_THALIANA -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2c_type-miseq_db-thaliana_query-lyrata.out
#	2d MiSeq lyrata vs DB lyrata (FP)
$BLASTN -db $BLAST_DB_LYRATA -query $READS_MISEQ_LYRATA  -num_threads 8 -evalue 1.0 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > $BLAST_OUT/2d_type-miseq_db-lyrata_query-lyrata.out

# 3. Pairwise comparisons - ONT
#	3a ONT thaliana, TP (thaliana)(1a) vs FP (lyrata)(1b)
$COMPARE_BLAST $BLAST_OUT/1a_type-pores_db-thaliana_query-thaliana.out $BLAST_OUT/1b_type-pores_db-lyrata_query-thaliana.out > $BLAST_OUT/pairwise_type-pores_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_type-pores_tp-thaliana_fp-lyrata.out.summary
#	3b ONT lyrata, TP (lyrata)(1d) vs FP (thaliana)(1c)
$COMPARE_BLAST $BLAST_OUT/1d_type-pores_db-lyrata_query-lyrata.out $BLAST_OUT/1c_type-pores_db-thaliana_query-lyrata.out > $BLAST_OUT/pairwise_type-pores_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_type-pores_tp-lyrata_fp-thaliana.out.summary

# 4. Pairwise comparisons - MiSeq
#	4a MiSeq thaliana, TP (thaliana)(2a) vs FP (lyrata)(2b)
$COMPARE_BLAST $BLAST_OUT/2a_type-miseq_db-thaliana_query-thaliana.out $BLAST_OUT/2b_type-miseq_db-lyrata_query-thaliana.out > $BLAST_OUT/pairwise_type-miseq_tp-thaliana_fp-lyrata.out 2> $BLAST_OUT/pairwise_type-miseq_tp-thaliana_fp-lyrata.out.summary
#	4b MiSeq lyrata, TP (lyrata)(2d) vs FP (thaliana)(2c)
$COMPARE_BLAST $BLAST_OUT/2d_type-miseq_db-lyrata_query-lyrata.out $BLAST_OUT/2c_type-miseq_db-thaliana_query-lyrata.out > $BLAST_OUT/pairwise_type-miseq_tp-lyrata_fp-thaliana.out 2> $BLAST_OUT/pairwise_type-miseq_tp-lyrata_fp-thaliana.out.summary
