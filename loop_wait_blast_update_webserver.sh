#! /bin/bash
# run from 
#~/Public/KewSciFest/mini-metrichor					: raw reads pulled from minion-1 by dir watcher
#~/Public/KewSciFest/mini-metrichor/uploaded		: uploaded by metrichor
#~/Public/KewSciFest/mini-metrichor/downloads		: fast5 d/l by metrichor
#~/Public/KewSciFest/data							: fasta files, unknown_sample.fa and the blast results

while true;
do
echo `date` > /Library/WebServer/Documents/realtimephylogenomics/data/SF.txt;

# get stats in poretools 
python2.7 /Users/joeparker/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py stats mini-metrichor/downloads/ > data/unknown_sample.stats 
# poretools fasta from the local metrichor d/l to the unknown_reads.fa
python2.7 /Users/joeparker/Documents/all_work/programming/repo-git/poretools/poretools/poretools_main.py fasta ~/Public/KewSciFest/mini-metrichor/downloads > ~/Public/KewSciFest/data/unknown_sample.fa

# blast the new reads from the local unknown_sample.fa:
 for db in {'erycina','sorbus','beta','napenthes','silene','lyrata'};
 #for db in {'erycina','sorbus','beta','napenthes','silene','thaliana','lyrata'};#lyrata only for now
 do
	blastn -db  /Volumes/LaCie/minION_analyses/sci-fest/presequenced_data/$db     -query /Users/joeparker/Public/KewSciFest/data/unknown_sample.fa -num_threads 8 -evalue 0.0001 -outfmt "6 qacc length pident evalue" -max_target_seqs 1  -max_hsps 1  | /Users/joeparker/Documents/all_work/programming/repo-git/real-time-phylogenomics/sum_blast.pl - > data/agg.blasts.$db;
	echo $db >> data/agg.blasts.$db;
 done;

# copy the blast results and read stats to rtp
cp data/agg.blasts.* /Library/WebServer/Documents/realtimephylogenomics/data/
cp data/unknown_sample.stats /Library/WebServer/Documents/realtimephylogenomics/data/

#sleep
sleep 10
done;


