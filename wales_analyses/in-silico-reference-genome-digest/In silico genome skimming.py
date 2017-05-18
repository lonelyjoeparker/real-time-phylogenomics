
# coding: utf-8

# In[2]:

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import random
from math import floor
from random import randint
import warnings


# In[39]:

# subsampling a genome according to the exponential, with mean=~ 1kbp


# In[11]:

def make_subsample(genome, label, num_draws):
    partial_subsample = []
    assembly  = list(SeqIO.parse(genome,'fasta')) # cast to a list; we need to address individual elements more than iterate reliably
    samples = num_draws # this is how many reads from each scaffold - nb should be a global increment but... meh
    outfile = 'subsample_'+label+'_'+str(samples).zfill(7)+'_draws.fasta'
    print('making '+str(samples)+' draws, output to '+outfile)
    for i in range(0,samples):
        # randomly pick (uniform distribution) one of the scaffolds/contigs
        # to draw this sample from 
        which_scaffold = randint(0,len(assembly)-1)
        scaffold = assembly[which_scaffold]
        #print(scaffold.id)
        some_new_seq = False
        nonzero_length_check = False
        while not nonzero_length_check:
            new_length = int(random.exponential(1500))
            start_pos = int(random.uniform(1,(len(scaffold)-new_length)))
            some_new_seq = scaffold[start_pos:(start_pos+new_length)]
            some_new_seq.id = some_new_seq.id + "_sampled_" + str(start_pos) + "_" + str(new_length)
            # WARNINGS if short (verbose)
            #if not (len(some_new_seq) >0):
            #warnings.warn('short sequence: '+str(len(some_new_seq))+' id '+some_new_seq.id)
            nonzero_length_check = (len(some_new_seq) >0)
        partial_subsample.append(some_new_seq)

    if(SeqIO.write(partial_subsample,outfile,'fasta')):
        print('wrote out '+str(len(partial_subsample))+' sequences to '+outfile)
    else:
        warnings.warn('COULD NOT write '+str(len(partial_subsample))+' sequences to '+outfile)

    del partial_subsample[:]
    pass


# input dir: /media/joe/BiSlDi/genomes/Arabidopsis_thaliana
# input file: arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna

genome = "/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna"

make_subsample(genome,'A.thal',10)
make_subsample(genome,'A.thal',100)
make_subsample(genome,'A.thal',1000)
make_subsample(genome,'A.thal',10000)
make_subsample(genome,'A.thal',100000)

# should now have input files with sub/resampled reads, exponentially distributed read lengths with median 1500,
# uniform picks from all seven A. thaliana TAIR10 chromosomes (e.g. mitochondrial and chloroplast genomes at same
# frequency as nuclear DNA)
# e.g. (example - random process means exact sizes will vary)
#
#joe-Tower:in-silico-reference-genome-digest (master*) joe$ lll *fasta
#-rw-rw-r-- 1 joe joe 1.1G May 18 12:19 subsample_100000_draws.fasta
#-rw-rw-r-- 1 joe joe 109M May 18 12:18 subsample_10000_draws.fasta
#-rw-rw-r-- 1 joe joe  11M May 18 12:18 subsample_1000_draws.fasta
#-rw-rw-r-- 1 joe joe 1.2M May 18 12:18 subsample_100_draws.fasta
#-rw-rw-r-- 1 joe joe 104K May 18 12:18 subsample_10_draws.fasta



# In[12]:

# repeat for Arabidopsis lyrata 

genomic_data = "/media/joe/BiSlDi/genomes/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt"

make_subsample(genomic_data,'A.lrya',10)
make_subsample(genomic_data,'A.lyra',100)
make_subsample(genomic_data,'A.lyra',1000)
make_subsample(genomic_data,'A.lyra',10000)
make_subsample(genomic_data,'A.lyra',100000)


# In[30]:

# then make a blast DB for each subsampled file, e.g.
#
# if: 
#-rw-rw-r-- 1 joe joe  102K May 18 14:56 subsample_A.thal_0000010_draws.fasta
#
# then run:
# makeblastdb -dbtype nucl -out A.thal.subsample_10 -title 'subsampled A.thaliana genome, 10 draws per chromosome' -in subsample_A.thal_0000010_draws.fasta 
#
# to give:
#-rw-rw-r-- 1 joe joe   11K May 18 15:02 A.thal.subsample_10.nhr
#-rw-rw-r-- 1 joe joe   964 May 18 15:02 A.thal.subsample_10.nin
#-rw-rw-r-- 1 joe joe   24K May 18 15:02 A.thal.subsample_10.nsq
#
# (the database will be referenced by 'A.thal.subsample_10' in blastn etc.)
#
# to automate this run:

#for i in *fasta;
#do
#echo $i; 
#makeblastdb -dbtype nucl -out $i -title '$i: subsampled A.thaliana genome, 10 draws per chromosome' -in $i;
#done

# to crunch through, etc:
#makeblastdb -dbtype nucl -out A.thal.subsample_10E1 -title 'subsampled A.thaliana genome, 10E1 draws per chromosome' -in subsample_A.thal_0000010_draws.fasta 
#makeblastdb -dbtype nucl -out A.thal.subsample_10E2 -title 'subsampled A.thaliana genome, 10E2 draws per chromosome' -in subsample_A.thal_0000100_draws.fasta 
#makeblastdb -dbtype nucl -out A.thal.subsample_10E3 -title 'subsampled A.thaliana genome, 10E3 draws per chromosome' -in subsample_A.thal_0001000_draws.fasta 


# In[ ]:




# In[ ]:



