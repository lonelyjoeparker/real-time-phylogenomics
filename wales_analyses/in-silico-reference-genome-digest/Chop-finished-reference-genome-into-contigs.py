
# coding: utf-8

# In[1]:

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import random, mean, median
from math import floor
import warnings
import matplotlib.pyplot as plt
import numpy as np

screwy_sample = SeqIO.parse('/home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/in-silico-reference-genome-digest/subsample_A.lyra_0001000_draws.fasta','fasta')
thal_genome = SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna",'fasta')
alyr_genomic_data = SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt",'fasta')


seq_lens = []

for seq in screwy_sample:
    seq_lens.append(len(seq))
    pass


seq_lens_t = []

for seq in thal_genome:
    seq_lens_t.append(len(seq))
    pass


seq_lens_l = []

for seq in alyr_genomic_data:
    seq_lens_l.append(len(seq))
    pass

    
len(seq_lens)
print(seq_lens[1:10])
print(len(seq_lens))
print(mean(seq_lens))
print(mean(seq_lens_t))
print(mean(seq_lens_l))
print(median(seq_lens_l))


# In[ ]:


plt.hist(seq_lens_l,bins=range(min(seq_lens_l), max(seq_lens_l) + 500, 500))
plt.show()


# In[10]:

thal_genome = SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna",'fasta')
for chromosome in thal_genome:
    print(chromosome.id, str(len(chromosome)))


# In[3]:

# divide a sequence up 

all_contigs_bin = []
all_contigs_lengths = []

for scaffold in SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna",'fasta'):
    # chop up this scaffold
    contigs_lengths = []
    start_pos=0
    total_length = len(scaffold)
    remaining_length = total_length - start_pos
    # we'll propose new length pieces, exponentially distributed, from the 5' until no scaffold left
    while remaining_length > 0:
        valid_length = False
        while not valid_length:
            new_length = int(random.exponential(35000))
            valid_length = (new_length>300)
        end_pos = start_pos + new_length
        if end_pos < total_length:
            new_piece=scaffold[start_pos:end_pos]
            start_pos = end_pos+1
            remaining_length = total_length - start_pos
            contigs_lengths.append(len(new_piece))
            all_contigs_lengths.append(len(new_piece))
            all_contigs_bin.append(new_piece)
        else:
            new_piece=scaffold[start_pos:]
            contigs_lengths.append(len(new_piece))
            all_contigs_lengths.append(len(new_piece))
            all_contigs_bin.append(new_piece)
            print('last piece ',str(len(new_piece)))
            remaining_length=0 # run out of sequence
    
    print(scaffold.id,' mean ',str(mean(contigs_lengths)),' N ',str(len(contigs_lengths)))
print('overall mean ',str(mean(all_contigs_lengths)),' N ',str(len(all_contigs_lengths)))

SeqIO.write(all_contigs_bin,'A.thal.downsample_reference.fa','fasta')

#plt.hist(all_contigs_lengths,bins=range(min(all_contigs_lengths), max(all_contigs_lengths) + 500, 500))
#plt.show()


# In[3]:

sentence='this is a valid_le'
sentence[4:]

