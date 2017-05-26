#commands from tues night 20160524

# housekeeping...
#
#first rebasecall all R7 in one go using nanocall:

docker run -v /Users/joeparker/Public/:/drop -v /Users/joeparker/Public/basecalled/:/basecalled --rm nanocall /drop/R7-original -o /basecalled/R7-all-nanocall.fasta -t 4 2>/dev/null

#move all to minion dir on LaCie
cp basecalled/R7-all-nanocall.fasta /Volumes/LaCie/minION_analyses/wales_20160518/analyses

## BLAST ##
# get some stats and extract fasta from R9
#
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py fasta main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/pass 2> /dev/null > analyses/R9-partial-metrichor.fasta
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py fasta main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/fail 2> /dev/null >> analyses/R9-partial-metrichor.fasta
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py qualdist main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/pass/ 2> /dev/null > analyses/read_stats/R9-partial-metrichor.qualdist.tdf
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py qualdist main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/fail/ 2> /dev/null >> analyses/read_stats/R9-partial-metrichor.qualdist.tdf
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py index main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/pass/ 2> /dev/null > analyses/read_stats/R9-partial-metrichor.index.tdf
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py index main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/fail/ 2> /dev/null >> analyses/read_stats/R9-partial-metrichor.index.tdf


#blast not very good... maybe try other params and/or BLAT...
#
blastn -db ../Arabidopsis_lyrata_petraea/A_lyrata -query analyses/R7-all-nanocall.fasta -outfmt "6 length pident evalue qacc"|wc
#a few, maybe 30
blastn -db ../Arabidopsis_thaliana/A_thaliana -query analyses/R7-all-nanocall.fasta -outfmt "6 length pident evalue qacc" -max_target_seqs 1
#a few, maybe 30
blastn -db ../Arabidopsis_lyrata_petraea/A_lyrata -query analyses/R9-partial-metrichor.fasta -outfmt "6 length pident evalue qacc" -max_target_seqs 1|wc
#      45     180    3649
blastn -db ../Arabidopsis_thaliana/A_thaliana -query analyses/R9-partial-metrichor.fasta -outfmt "6 length pident evalue qacc" -max_target_seqs 1|wc
#      29     116    2348


## MAPPING ##
# index the input genomes to map 
cd /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/
bwa index Arabidopsis_lytata_ADBK01.1.fsa_nt
cd /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/
bwa index arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna

# map EVERYTHING (all R7 or R9 runs, *NOT* sorted by source organism!!!) to lyrata/thaliana
# http://angus.readthedocs.io/en/stable/analyzing_nanopore_data.html
cd /Volumes/LaCie/minION_analyses/wales_20160518/analyses/mapping
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna ../R7-all-nanocall.fasta | samtools view -bS - | samtools sort - R7-A_thaliana-mapped-reads.sorted
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna ../R9-partial-metrichor.fasta | samtools view -bS - | samtools sort - R9-A_thaliana-mapped-reads.sorted
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt ../R7-all-nanocall.fasta | samtools view -bS - | samtools sort - R7-A_lyrata-mapped-reads.sorted
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt ../R9-partial-metrichor.fasta | samtools view -bS - | samtools sort - R9-A_lyrata-mapped-reads.sorted

# don't forget to index for BAMView etc
samtools index R7-A_thaliana-mapped-reads.sorted.bam 
samtools index R7-A_lyrata-mapped-reads.sorted.bam 
samtools index R9-A_thaliana-mapped-reads.sorted.bam 
samtools index R9-A_lyrata-mapped-reads.sorted.bam 

# calculate avg. coverage
# https://www.biostars.org/p/5165/#67920 NB this will OMIT bases that aren't covered - so almost CERTAIN TO BE TOO HIGH
samtools depth R7-A_thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#Average =  0.897763
samtools depth R9-A_thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#Average =  2.52333
samtools depth R7-A_lyrata-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#Average =  0.903082
samtools depth R9-A_lyrata-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
#Average =  2.36002
#
# CORRECT genome sizes:
# A lyrata ssp. petraea ~ 203Mb http://www.ncbi.nlm.nih.gov/genome/genomes/493
# A thaliana ~ 120Mb http://www.ncbi.nlm.nih.gov/genome/genomes/4
samtools depth R9-A_lyrata-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/203000000}'
#Average =  0.00604587
samtools depth R7-A_lyrata-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/203000000}'
#Average =  0.00026701
samtools depth R9-A_thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/120000000}'
#Average =  0.0099659
samtools depth R7-A_thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/120000000}'
#Average =  0.000710983

#these aren't very high... maybe trim an adapter-sized chunk (?20-30bp) from 5' and 3' ..?

## BY SAMPLE ##

# split by name http://stackoverflow.com/questions/17368872/how-to-move-or-copy-files-listed-by-find-command-in-unix
# to find /path/to/search/ -type f -name "regular-expression-to-find-files" | xargs cp -t /target/path/
# TURNS OUT R9-partial are ALL LYRATA
# eg
find R9-original-PARTIAL_TUES_PM/reads/downloads/fail/ -type f -name "*lyrata*"|wc
#    5037    5037  550606
find R9-original-PARTIAL_TUES_PM/reads/downloads/fail/ -type f -name "*fast5*"|wc
#    5037    5037  550606
find R7-original/reads/uploaded/ -type f -name "*thaliana*" -exec cp {} R7-original/by-name-at \;
find R7-original/reads/uploaded/ -type f -name "*lyr*" -exec cp {} R7-original/by-name-ap \;
#blast R7 reads.. there aren't any R7 lyrata nanocalled...
blastn -db ../../Arabidopsis_thaliana/A_thaliana -query R7-thaliana-nanocall.fasta -outfmt "6 length pident evalue qacc"|wc
#      30     120    3144
blastn -db ../../Arabidopsis_lyrata_petraea/A_lyrata -query R7-thaliana-nanocall.fasta -outfmt "6 length pident evalue qacc"|wc
#      19      76    1988
# more reads blast thaliana than lyrata, better evals too.
#
# map R7 thaliana against thaliana AND petraea
#map
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna ../R7- | samtools view -bS - | samtools sort - R7-A_thaliana-thaliana-mapped-reads.sorted
bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/ ../R7-thaliana-nanocall.fasta | samtools view -bS - | samtools sort - R7-A_lyrata-thaliana-mapped-reads.sorted
#index
samtools index R7-A_thaliana-thaliana-mapped-reads.sorted.bam 
samtools index R7-A_lyrata-thaliana-mapped-reads.sorted.bam 
##coverage of thaliana R7 on thaliana
samtools depth R7-A_thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/120000000}'
#	Average =  0.000684383
##coverage of thaliana R7 on lyrata
samtools depth R7-A_lyrata-thaliana-mapped-reads.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/203000000}'
#	Average =  0.000251778

# hmm... coverage of thaliana on thaliana is ~2x 3x better than thaliana on lyrata. not conclusive but..

# all looks good except why is it so poo!?
# how many of the R9 failing reads are short?
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py fasta main_backup_copies_from_raw/R9-original-PARTIAL_TUES_/reads/downloads/fail 2> /dev/null|grep -e ">"|wc
#    4995   14985 1276825
# 5000 R9 in fail dir
python2.7 ~/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py fasta main_backup_copies_from_raw/R9-original-PARTIAL_TUES_PM/reads/downloads/fail --min-length 500 2> /dev/null|grep -e ">"|wc
#    3126    9378  798144
# of which > 3000 have 0.5kb or longer, surely not adapters...
# issues with R9 https://wiki.nanoporetech.com/display/FMDP/2016/05/10/R9+on+the+Plain?focusedCommentId=36992745#comment-36992745
# adapters for R9 rapid at https://wiki.nanoporetech.com/display/FMDP/Rapid+sequencing+protocol+with+SQK-RAD001
# check and blast for adapters
makeblastdb -in adapters.fasta -dbtype nucl
blastn -db adapters.fasta -query R9-partial-metrichor.fasta -outfmt "6 qacc sacc"
# <20 hits from >5000 sequences.. clearly blast is not working well with default params, or maybe this is garbage