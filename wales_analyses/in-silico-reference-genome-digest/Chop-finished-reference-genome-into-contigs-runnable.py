
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
import argparse, os
#import os.path as path

# subclass argparse to check dir path:
# code from https://codereview.stackexchange.com/questions/28608/checking-if-cli-arguments-are-valid-files-directories-in-python
# Author @aclewis182 https://github.com/aclewis182
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname
# /end https://codereview.stackexchange.com/questions/28608/checking-if-cli-arguments-are-valid-files-directories-in-python

# set up an argument parser
N50 = -1
parser = argparse.ArgumentParser(description='Chop a genome assembly up into pooer-quality chunks')
parser.add_argument('N50',help='target output N50 (Kbp, approximate)',nargs=1,type=int)
parser.add_argument('output_dir',help='target output directory',type=is_dir, action=FullPaths)

# hardcoded genome paths
thal_genome_path = "/media/joe/BiSlDi/genomes/Arabidopsis_thaliana/arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna"
lyra_genome_path = "/media/joe/BiSlDi/genomes/Arabidopsis_lyrata_petraea/Arabidopsis_lytata_ADBK01.1.fsa_nt"

# parse args to get target N50 and output directory
inputs = parser.parse_args()
N50 = inputs.N50
output_dir = inputs.output_dir

print('writing to ',inputs.output_dir)

def chop_up_genomes(starting_genome, target_N50, target_output,target_label):
    # divide a sequence up

    all_contigs_bin = []
    all_contigs_lengths = []

    for scaffold in SeqIO.parse(starting_genome,'fasta'):
        # chop up this scaffold
        new_mean = int(int(target_N50[0]) * (35.0/55.0))  # fudge factor for exponential distribution mean -> N50 stat
        contigs_lengths = []
        start_pos=0
        total_length = len(scaffold)
        remaining_length = total_length - start_pos
        # we'll propose new length pieces, exponentially distributed, from the 5' until no scaffold left
        while remaining_length > 0:
            valid_length = False
            while not valid_length:
                new_length = int(random.exponential(new_mean))
                valid_length = (new_length>300) # minimum acceptable length
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
    print('overall mean :',str(mean(all_contigs_lengths)),' N ',str(len(all_contigs_lengths)))

    new_output = os.path.join(target_output,str('downsample_'+target_label+'['+str(int(mean(all_contigs_lengths)))+'].fasta'))
    print('Writing to: ',type(new_output),new_output)
    SeqIO.write(all_contigs_bin,new_output,'fasta')

    pass

print(type(output_dir[0]))
chop_up_genomes(thal_genome_path,N50,output_dir,'some_thaliana_junk')
chop_up_genomes(lyra_genome_path,N50,output_dir,'some_alyrata_junk')
