#! /bin/bash

cd ~/Documents/all_work/programming/metrichor/downloads
../extract_fastq_v2.sh ../converted_reads/out.fastq

#try and add a random number to out.fastq... 
somenum=`awk -v min=0 -v max=999 'BEGIN{srand(); print int(min+rand()*(max-min+1))}'`
echo -e "@${somenum}_taxon\nCTTTTTAATTG\n+\nHVQF7#%0(*," >> ../converted_reads/out.fastq