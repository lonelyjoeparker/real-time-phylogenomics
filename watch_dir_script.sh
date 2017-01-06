#!/bin/bash
#get the args we need
#invoke with e.g 
#java -jar ~/Documents/all_work/programming/metrichor/SimplifiedWatchDir.jar /Volumes/F-sci-fest-B-minion-1/data/20160805_sample_01/ /bin/bash /Users/joeparker/Public/KewSciFest/watch_dir_script.sh /Users/joeparker/Public/KewSciFest/mini-metrichor
#nb CHECK PATHS

# /Users/joeparker/Public/KewSciFest								: scripts etc executed from here
# /Users/joeparker/Public/KewSciFest/mini-metrichor					: raw reads pulled from minion-1 by dir watcher
# /Users/joeparker/Public/KewSciFest/mini-metrichor/uploaded		: uploaded by metrichor
# /Users/joeparker/Public/KewSciFest/mini-metrichor/downloads		: fast5 d/l by metrichor
# /Users/joeparker/Public/KewSciFest/data							: fasta files, unknown_sample.fa and the blast results

raw=$1
output=$2 #not used at the moment but there if we want it..

echo [nanocall-docker]: try to basecall $raw

# copy files from /Volumes/Minion/reads/whereever to /Users/joeparker/Public/reads etc
cp  /Volumes/F-sci-fest-B-minion-1/data/20160805_sample_01/reads/$raw /Users/joeparker/Public/KewSciFest/mini-metrichor/$raw	
