
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
import argparse

# set up an argument parser
N50 = -1
parser = argparse.ArgumentParser(description='Chop a genome assembly up into pooer-quality chunks')
parser.add_argument('N50',help='target output N50 (Kbp, approximate)',nargs=1,type=int)
parser.add_argument('output_dir',help='target output directory',nargs='?', type=argparse.FileType('d'))
screwy_sample = SeqIO.parse('/home/joe/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/in-silico-reference-genome-digest/subsample_A.lyra_0001000_draws.fasta','fasta')
thal_genome = SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna",'fasta')
alyr_genomic_data = SeqIO.parse("/media/joe/BiSlDi/genomes/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt",'fasta')

N50 = parser.parse_args().N50
output_dir = parser.parse_args().output_dir

print(N50,output_dir)

quit()

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


# In[2]:

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
