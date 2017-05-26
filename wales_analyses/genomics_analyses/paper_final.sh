### final paper analyses

## collect all inputs

cd /Users/joeparker/Downloads/wales_nelumbolab_output
# lyrata
cat _Volumes_LaCie_wales_primary_data_collection_R7_A_lyrata_.poretools.fasta _Volumes_LaCie_wales_primary_data_collection_R9_A_lyrata_.poretools.fasta > all_R7R9_petraea.fasta 
cp all_R7R9_petraea.fasta /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_petraea.fasta
# thaliana
cat R7_A_thaliana.poretools.fasta _Volumes_SCI-FEST-A_R9_A_thaliana.poretools.fasta > all_R7R9_thaliana.fasta
cp all_R7R9_thaliana.fasta /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta


## map with BWA

# ONT petraea on lyrata
# /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt
# don't forget to indext reference first if it isn't done:
# bwa index /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt
bwa mem -t6 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_petraea.fasta > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_lyrata.mapped-reads.bam
samtools view -bh  /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_lyrata.mapped-reads.bam | samtools sort - /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_lyrata.mapped-reads.sorted
samtools index /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_lyrata.mapped-reads.sorted.bam 
samtools depth /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_lyrata.mapped-reads.sorted.bam  |  awk '{sum+=$3} END { print "Average mapping A. petraea on A. lyrata = ",sum/NR}' 
# Average mapping A. petraea on A. lyrata =  4.07007

# ONT petraea on petraea
# /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt
# don't forget to indext reference first if it isn't done:
# bwa index /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt
bwa mem -t7 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_petraea.fasta > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_petraea.mapped-reads.bam
samtools view -bh  /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_petraea.mapped-reads.bam | samtools sort - /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_petraea.mapped-reads.sorted
samtools index /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_petraea.mapped-reads.sorted.bam 
samtools depth /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_petraea.mapped-reads.sorted.bam  |  awk '{sum+=$3} END { print "Average mapping A. petraea on A. petraea  = ",sum/NR}' 
# Average mapping A. petraea on A. petraea  =  4.35829


# ONT thaliana
# don't forget to indext reference first if it isn't done:
# bwa index Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna 
#bwa mem -t3 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta | samtools view -bS - | samtools sort - /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.sorted
bwa mem -t6 -x ont2d /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.bam
samtools view -bh  /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.bam | samtools sort - /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.sorted
samtools index /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.sorted.bam 
samtools depth /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9-A_thaliana.mapped-reads.sorted.bam  |  awk '{sum+=$3} END { print "Average mapping A. thaliana on A. thaliana  = ",sum/NR}' 
# Average =  1.81726

# MiSeq reads/reports are in: /mnt/HDD_2/joe/WALES-NELUMBOLAB
# MiSeq trimmed reads to actually USE are /mnt/HDD_2/joe/WALES-NELUMBOLAB/arabidopsis_trimmed 
# 645M Sep 13 16:19 AL1a_S3_L001_1P.fq.gz	- forward, paired trimmed reads (see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
# 211M Sep 13 16:13 AL1a_S3_L001_1U.fq.gz	- forward, unpaired trimmed reads
# 731M Sep 13 16:14 AL1a_S3_L001_2P.fq.gz	- reverse, paired trimmed reads
#  13M Sep 13 16:14 AL1a_S3_L001_2U.fq.gz	- reverse, unpaired trimmed reads
# 560M Sep 13 16:15 AL2a_S4_L001_1P.fq.gz
# 177M Sep 13 16:13 AL2a_S4_L001_1U.fq.gz
# 631M Sep 13 16:13 AL2a_S4_L001_2P.fq.gz
#  11M Sep 13 16:13 AL2a_S4_L001_2U.fq.gz
# 697M Sep 13 16:18 AT1a_S1_L001_1P.fq.gz
# 192M Sep 13 16:18 AT1a_S1_L001_1U.fq.gz
# 791M Sep 13 16:17 AT1a_S1_L001_2P.fq.gz
#  11M Sep 13 16:14 AT1a_S1_L001_2U.fq.gz
# 192M Sep 13 16:19 AT2a_S2_L001_1U.fq.gz
# 9.6M Sep 13 16:15 AT2a_S2_L001_2U.fq.gz

# README info says: 
#	T1 and AT2 are thaliana. Differnet individuals. AT2 was used for minION.
#	AL1 and AL2 are Lyrata. Different rosettes but quite likely to be the same plant (not definite). AL1 was used for MinION.
# MiSeq petraea on lyrata
# MiSeq petraea on petraea
# MiSeq thaliana

# first interleave paired reads
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/arabidopsis_trimmed 
interleave-reads.py AT1a_S1_L001_1P.fq.gz AT1a_S1_L001_2P.fq.gz --gzip -o AT1a_S1_L001_paired.interleaved.fq.gz
interleave-reads.py AT2a_S2_L001_1P.fq.gz AT2a_S2_L001_2P.fq.gz --gzip -o AT2a_S2_L001_paired.interleaved.fq.gz 
interleave-reads.py AL1a_S3_L001_1P.fq.gz AL1a_S3_L001_2P.fq.gz -o AL1a_S3_L001_paired.interleaved.fq.gz --gzip
interleave-reads.py AL1a_S3_L001_1P.fq.gz AL1a_S3_L001_2P.fq.gz -o AL1a_S3_L001_paired.interleaved.fq.gz --gzip
# then concatenate unpaired reads
cat AT1a_S1_L001_*U*gz > AT1a_S1_L001_unpaired.concatenated.fq.gz
cat AT2a_S2_L001_*U*gz > AT2a_S2_L001_unpaired.concatenated.fq.gz
cat AL1a_S3_L001_*U*gz > AL1a_S3_L001_unpaired.concatenated.fq.gz
cat AL2a_S4_L001_*U*gz > AL2a_S4_L001_unpaired.concatenated.fq.gz 
# then get quick stats
for i in ls *P*gz; do echo $i ;gunzip -c $i | grep -e "@M" - | wc; done > paired.read.counts.tdf
for i in ls *U*gz; do echo $i ;gunzip -c $i | grep -e "@M" - | wc; done > unpaired.read.counts.tdf
for i in ls *interleaved*gz; do echo $i ;gunzip -c $i | grep -e "@M" - | wc; done > paired.interleaved.read.counts.tdf
for i in ls *concatenated*gz; do echo $i ;gunzip -c $i | grep -e "@M" - | wc; done > unpaired.concatenated.read.counts.tdf
# concatenate all for convenience (we'll use them for BLAST / LAST)
cat AL1a_S3_L001_paired.interleaved.fq.gz AL1a_S3_L001_unpaired.concatenated.fq.gz > AL1a_S3_L001_all.trimmed.fq.gz
cat AL2a_S4_L001_paired.interleaved.fq.gz AL2a_S4_L001_unpaired.concatenated.fq.gz > AL2a_S4_L001_all.trimmed.fq.gz
cat AT1a_S1_L001_paired.interleaved.fq.gz AT1a_S1_L001_unpaired.concatenated.fq.gz > AT1a_S1_L001_all.trimmed.fq.gz
cat AT2a_S2_L001_paired.interleaved.fq.gz AT2a_S2_L001_unpaired.concatenated.fq.gz > AT2a_S2_L001_all.trimmed.fq.gz
# extract .fasta for blast too
gunzip -c AL1a_S3_L001_all.trimmed.fq.gz | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' - > AL1a_S3_L001_all.trimmed.fa
gunzip -c AL2a_S4_L001_all.trimmed.fq.gz | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' - > AL2a_S4_L001_all.trimmed.fa
gunzip -c AT1a_S1_L001_all.trimmed.fq.gz | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' - > AT1a_S1_L001_all.trimmed.fa
gunzip -c AT2a_S2_L001_all.trimmed.fq.gz | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' - > AT2a_S2_L001_all.trimmed.fa

# now we can actually map
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/mapping_illumina/
# paired first

#AL1a, lyrata
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt -p ../arabidopsis_trimmed/AL1a_S3_L001_paired.interleaved.fq.gz > AL1a.lyrata.paired.sam
samtools view -bS AL1a.lyrata.paired.sam -o AL1a.lyrata.paired.bam
samtools sort AL1a.lyrata.paired.bam AL1a.lyrata.paired.sorted
samtools depth AL1a.lyrata.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL1a paired-end on A. lyrata  = ",sum/NR}' > AL1a.paired.lyrata.depth.txt
#AL1a, petraea
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt -p ../arabidopsis_trimmed/AL1a_S3_L001_paired.interleaved.fq.gz > AL1a.petraea.paired.sam
samtools view -bS AL1a.petraea.paired.sam -o AL1a.petraea.paired.bam
samtools sort AL1a.petraea.paired.bam AL1a.petraea.paired.sorted
samtools depth AL1a.petraea.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL1a paired-end on A. lyrata ssp. petraea = ",sum/NR}' > AL1a.paired.petraea.depth.txt
#AL2a, lyrata
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt -p ../arabidopsis_trimmed/AL2a_S3_L001_paired.interleaved.fq.gz > AL2a.lyrata.paired.sam
samtools view -bS AL2a.lyrata.paired.sam -o AL2a.lyrata.paired.bam
samtools sort AL2a.lyrata.paired.bam AL2a.lyrata.paired.sorted
samtools depth AL2a.lyrata.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL2a paired-end on A. lyrata  = ",sum/NR}' > AL2a.paired.lyrata.depth.txt
#AL2a, petraea
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt -p ../arabidopsis_trimmed/AL2a_S3_L001_paired.interleaved.fq.gz > AL2a.petraea.paired.sam
samtools view -bS AL2a.petraea.paired.sam -o AL2a.petraea.paired.bam
samtools sort AL2a.petraea.paired.bam AL2a.petraea.paired.sorted
samtools depth AL2a.petraea.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL2a paired-end on A. lyrata ssp. petraea = ",sum/NR}' > AL2a.paired.petraea.depth.txt
#AT1a
bwa mem -t 6 ../../REFERENCE-GENOMES/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna -p ../arabidopsis_trimmed/AT1a_S1_L001_paired.interleaved.fq.gz > AT1a.paired.sam
samtools view -bS AT1a.paired.sam -o AT1a.paired.bam
samtools sort AT1a.paired.bam AT1a.paired.sorted
samtools depth AT1a.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. thaliana AT1a paired-end on A. thaliana  = ",sum/NR}' > AT1a.paired.thaliana.depth.txt
#AT2a
bwa mem -t 6 ../../REFERENCE-GENOMES/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna -p ../arabidopsis_trimmed/AT2a_S2_L001_paired.interleaved.fq.gz > AT2a.paired.sam
samtools view -bS AT2a.paired.sam -o AT2a.paired.bam
samtools sort AT2a.paired.bam AT2a.paired.sorted
samtools depth AT2a.paired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. thaliana AT2a paired-end on A. thaliana  = ",sum/NR}' > AT2a.paired.thaliana.depth.txt

# then unpaired
#AL1a, lyrata
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt -p ../arabidopsis_trimmed/AL1a_S3_L001_unpaired.concatenated.fq.gz > AL1a.lyrata.unpaired.sam
samtools view -bS AL1a.lyrata.unpaired.sam -o AL1a.lyrata.unpaired.bam
samtools sort AL1a.lyrata.unpaired.bam AL1a.lyrata.unpaired.sorted
samtools depth AL1a.lyrata.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL1a unpaired reads on A. lyrata  = ",sum/NR}' > AL1a.unpaired.lyrata.depth.txt
#AL1a, petraea
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt -p ../arabidopsis_trimmed/AL1a_S3_L001_unpaired.concatenated.fq.gz > AL1a.petraea.unpaired.sam
samtools view -bS AL1a.petraea.unpaired.sam -o AL1a.petraea.unpaired.bam
samtools sort AL1a.petraea.unpaired.bam AL1a.petraea.unpaired.sorted
samtools depth AL1a.petraea.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL1a unpaired reads on A. lyrata ssp. petraea = ",sum/NR}' > AL1a.unpaired.petraea.depth.txt
#AL2a, lyrata
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt -p ../arabidopsis_trimmed/AL2a_S3_L001_unpaired.concatenated.fq.gz > AL2a.lyrata.unpaired.sam
samtools view -bS AL2a.lyrata.unpaired.sam -o AL2a.lyrata.unpaired.bam
samtools sort AL2a.lyrata.unpaired.bam AL2a.lyrata.unpaired.sorted
samtools depth AL2a.lyrata.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL2a unpaired reads on A. lyrata  = ",sum/NR}' > AL2a.unpaired.lyrata.depth.txt
#AL2a, petraea
bwa mem -t 6  ../../REFERENCE-GENOMES/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt -p ../arabidopsis_trimmed/AL2a_S3_L001_unpaired.concatenated.fq.gz > AL2a.petraea.unpaired.sam
samtools view -bS AL2a.petraea.unpaired.sam -o AL2a.petraea.unpaired.bam
samtools sort AL2a.petraea.unpaired.bam AL2a.petraea.unpaired.sorted
samtools depth AL2a.petraea.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. lyrata ssp. petraea AL2a unpaired reads on A. lyrata ssp. petraea = ",sum/NR}' > AL2a.unpaired.petraea.depth.txt
#AT1a
bwa mem -t 6 ../../REFERENCE-GENOMES/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna -p ../arabidopsis_trimmed/AT1a_S1_L001_unpaired.concatenated.fq.gz > AT1a.unpaired.sam
samtools view -bS AT1a.unpaired.sam -o AT1a.unpaired.bam
samtools sort AT1a.unpaired.bam AT1a.unpaired.sorted
samtools depth AT1a.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. thaliana AT1a unpaired reads on A. thaliana  = ",sum/NR}' > AT1a.unpaired.thaliana.depth.txt
#AT2a
bwa mem -t 6 ../../REFERENCE-GENOMES/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna -p ../arabidopsis_trimmed/AT2a_S2_L001_unpaired.concatenated.fq.gz > AT2a.unpaired.sam
samtools view -bS AT2a.unpaired.sam -o AT2a.unpaired.bam
samtools sort AT2a.unpaired.bam AT2a.unpaired.sorted
samtools depth AT2a.unpaired.sorted.bam |  awk '{sum+=$3} END { print "Average mapping A. thaliana AT2a unpaired reads on A. thaliana  = ",sum/NR}' > AT2a.unpaired.thaliana.depth.txt


## align with LAST

# ONT petraea on lyrata Arabidopsis_lyrata_petraea/A_lyrata-LAST
echo Mapping petraea to lyrata > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_lyrata.maf.report.txt
lastal /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata-LAST /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_petraea.fasta | last-split > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_lyrata.maf
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/countErrorsFromMaf.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_lyrata.maf >> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_lyrata.maf.report.txt
#	Total aligned length:	980358
#	Total aligned errors:	220907
#	Nominal aligned error rate:	0.225332990601393


# ONT petraea on petraea
echo Mapping petraea to petraea > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_petraea.maf.report.txt
lastal /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/A_petraea-LAST /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_petraea.fasta | last-split > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_petraea.maf
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/countErrorsFromMaf.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_petraea.maf >> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_petraea.maf.report.txt
#	Total aligned length:	811232
#	Total aligned errors:	190965
#	Nominal aligned error rate:	0.23540121691452

# ONT thaliana (DB is at /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana-LAST)
echo Mapping thaliana to thaliana > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_thaliana.maf.report.txt
lastal /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana-LAST /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta | last-split > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_thaliana.maf
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/countErrorsFromMaf.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_thaliana.maf >> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_thaliana.maf.report.txt
# see /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/mapping/all_R7R9_thaliana.maf.report.txt
#	Total aligned length:	53056319
#	Total aligned errors:	11106856
#	Nominal aligned error rate:	0.20934087040603


## LAST align illumina on haemodorum
# see also /mnt/HDD_2/joe/WALES-NELUMBOLAB/LAST_illumina/LAST-mapping.sh
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/LAST_illumina/

# MiSeq thaliana (AT2a) on thaliana
echo Mapping AT2a thaliana to thaliana > all_AT2a_thaliana.maf.report.txt
lastal /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_thaliana/A_thaliana-LAST /mnt/HDD_2/joe/WALES-NELUMBOLAB/arabidopsis_trimmed/AL1a_S3_L001_all.trimmed.fa | last-split > all_AT2a_thaliana.maf
perl /mnt/HDD_2/joe/WALES-NELUMBOLAB/countErrorsFromMaf.pl all_AT2a_thaliana.maf >> all_AT2a_thaliana.maf.report.txt

# MiSeq petraea (AL1a) on lyrata
echo Mapping AL1a petraea to lyrata >  all_AL1a_lyrata.maf.report.txt
lastal /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/A_lyrata-LAST /mnt/HDD_2/joe/WALES-NELUMBOLAB/arabidopsis_trimmed/AL1a_S3_L001_all.trimmed.fa | last-split > all_AL1a_lyrata.maf
perl /mnt/HDD_2/joe/WALES-NELUMBOLAB/countErrorsFromMaf.pl all_AL1a_lyrata.maf >> all_AL1a_lyrata.maf.report.txt

# MiSeq petraea (AL1a) on petraea
echo Mapping AL1a petraea to petraea > all_AL1a_petraea.maf.report.txt
lastal /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_petraea_contigs/A_petraea-LAST /mnt/HDD_2/joe/WALES-NELUMBOLAB/arabidopsis_trimmed/AL1a_S3_L001_all.trimmed.fa | last-split > all_AL1a_petraea.maf
perl /mnt/HDD_2/joe/WALES-NELUMBOLAB/countErrorsFromMaf.pl all_AL1a_petraea.maf >> all_AL1a_petraea.maf.report.txt


## Check that field-sequenced and lab-sequenced data are comparable

## Predict genes from sequencing

## de novo assembly
# dirs on both boletus and senecio:
#	/mnt/HDD_2/joe/WALES-NELUMBOLAB/SPAdes-petraea
#	/mnt/HDD_2/joe/WALES-NELUMBOLAB/SPAdes-thaliana
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/SPAdes-petraea
nohup spades.py -o ./hybrid_petraea  --12 ../arabidopsis_trimmed/AL1a_S3_L001_paired.interleaved.fq.gz -s ../arabidopsis_trimmed/AL1a_S3_L001_unpaired.concatenated.fq.gz --nanopore ../inputs/all_R7R9_petraea.fasta -m 210 -t 20 > hybrid_petraea.out 2> hybrid_petraea.err &
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/SPAdes-thaliana
nohup spades.py -o ./hybrid_thaliana --12 ../arabidopsis_trimmed/AT2a_S2_L001_paired.interleaved.fq.gz -s ../arabidopsis_trimmed/AT2a_S2_L001_unpaired.concatenated.fq.gz --nanopore ../inputs/all_R7R9_thaliana.fasta -m 210 -t 20 > hybrid_thaliana.out 2> hybrid_thaliana.err &

## BLAST and LAST TP/FP rates

# e-value 0.0001, ONT
#blast, R7R9 thaliana ONT queries, one HSP per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.0001-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.0001-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.0001-DB.out;
#blast, R7R9 thaliana ONT queries, keep all HSPs per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-allHits-evalue0.0001-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-allHits-evalue0.0001-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-allHits-evalue0.0001-DB.out;
#compare hits
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.0001-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.0001-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.0001.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.0001.summary;
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.0001-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.0001-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.0001.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.0001.summary;

# e-value 0.01, ONT
#blast, R7R9 thaliana ONT queries, one HSP per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.01-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.01-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.01-DB.out;
#blast, R7R9 thaliana ONT queries, keep all HSPs per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-allHits-evalue0.01-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-allHits-evalue0.01-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.01 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-allHits-evalue0.01-DB.out;
#compare hits
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.01-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.01-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.01.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.01.summary;
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.01-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.01-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.01.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.01.summary;

# e-value 0.1, ONT
#blast, R7R9 thaliana ONT queries, one HSP per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.1-DB.out;
#blast, R7R9 thaliana ONT queries, keep all HSPs per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-allHits-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-allHits-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-allHits-evalue0.1-DB.out;
#compare hits
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.1-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.1.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.1.summary;
perl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/LondonCalling-lastMinute/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.1-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.1.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.1.summary;

# e-value 0.1, Illumina
#blast, all thaliana Illumina queries, one HSP per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1 > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.1-DB.out;
#blast, all thaliana Illumina, keep all HSPs per query
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_thaliana/A_thaliana 													-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-allHits-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_lyrata_petraea/A_lyrata 												-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-allHits-evalue0.1-DB.out;
blastn -db /Volumes/LaCie/minION_analyses/Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt 	-query /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/inputs/all_R7R9_thaliana.fasta -num_threads 8 -evalue 0.1 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-allHits-evalue0.1-DB.out;
#compare hits
perl /mnt/HDD_2/joe/WALES-NELUMBOLAB/BLAST_illumina/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-petraea-singleHit-evalue0.1-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.1.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-petraea-evalue0.1.summary;
perl /mnt/HDD_2/joe/WALES-NELUMBOLAB/BLAST_illumina/compareBlastHits.pl /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-thaliana-singleHit-evalue0.1-DB.out /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-on-lyrata-singleHit-evalue0.1-DB.out > /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.1.out 2> /Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/blast/all_R7R9_thaliana-vs-lyrata-evalue0.1.summary;

# we'll repeat blasting but first merge both A. lyrata subspecies as neither DB seems to be that good...
cd /mnt/HDD_2/joe/REFERENCE-GENOMES
mkdir Arabidopsis_both_ssp_merged
cat Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt Arabidopsis_petraea_contigs/Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt | makeblastdb -in - -input_type fasta -title Arabidopsis_lyrata_merged_ssp -out Arabidopsis_both_ssp_merged/A_lyrata_ssp -dbtype nucl
blastdbcmd -db Arabidopsis_both_ssp_merged/A_lyrata_ssp -info # check the DB has built as expected
# now repeat blasting but with evalue thresh = 1.0 (pretty massively lax), still only best hit from each query
# this will then be used to calculate ROC curves for:
#	1.0 < evalue ≤ 0
#	0.75 < pident % ≤ 100
#	150 < length < 100,000 (this last will obviously only be complete for ONT)
#
# in each case we can do simple pairwise A. thaliana vs A. lyrata for ONT / Illumina


## de novo assembly of hybrid ONT / Illumina data:

# thaliana
# SPAdes assembly:
nohup spades.py -o ./hybrid_thaliana --12 ../arabidopsis_trimmed/AT2a_S2_L001_paired.interleaved.fq.gz -s ../arabidopsis_trimmed/AT2a_S2_L001_unpaired.concatenated.fq.gz --nanopore ../inputs/all_R7R9_thaliana.fasta -m 210 -t 20 > hybrid_thaliana.out 2> hybrid_thaliana.err &
# Quast 4.3 stats vs TAIR thaliana genome:
quast.py hybrid_thaliana/contigs.fasta hybrid_thaliana/scaffolds.fasta -o hybrid_thaliana_quast_reports/ -R /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna

# petraea
nohup spades.py -o ./hybrid_petraea  --12 ../arabidopsis_trimmed/AL1a_S3_L001_paired.interleaved.fq.gz -s ../arabidopsis_trimmed/AL1a_S3_L001_unpaired.concatenated.fq.gz --nanopore ../inputs/all_R7R9_petraea.fasta -m 210 -t 20 > hybrid_petraea.out 2> hybrid_petraea.err &
# Quast 4.3 stats vs lyrata ADBK 0.1.1 genome
quast.py hybrid_petraea/contigs.fasta hybrid_petraea/scaffolds.fasta -o hybrid_petraea_quast_reports/ -R /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt

# Quast 4.3 stats vs petraea contics
##quast.py hybrid_petraea/contigs.fasta hybrid_petraea/scaffolds.fasta -o hybrid_petraea_quast_reports/ -R /mnt/HDD_2/joe/REFERENCE-GENOMES/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt

#try to run CEGMA
cd /mnt/HDD_2/joe/WALES-NELUMBOLAB/CEGMA/CEGMA_v2-master
nohup cegma -T 1 --genome /mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/all_R7R9_petraea.fasta --output CEGMA_lyrata > CEGMA_lyrata.nohup.out 2>CEGMA_lyrata.nohup.err &
nohup cegma -T 1 --genome /mnt/HDD_2/joe/WALES-NELUMBOLAB/inputs/all_R7R9_thaliana.fasta --output CEGMA_thaliana > CEGMA_thaliana.nohup.out 2>CEGMA_thaliana.nohup.err &
