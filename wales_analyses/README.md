# Parker et al (2017 in prep) analyses

## Introduction

This repository contains code, instructions and listed dependencies needed to re-run the analyses presented in our field-sequencing paper (preprint: http://biorxiv.org/content/early/2017/03/14/107656).

We are in the process of fully documenting this code to make these scripts easier to run on your own system, please note however that producing a multifunctional analysis tool is not our intention (rather enhancing the reproducibility of this paper). This work partly involves altering where certain executables and input files are located to make porting the scripts simpler; as a result later [commits](https://github.com/lonelyjoeparker/real-time-phylogenomics/commits/master) may alter the various component scripts slightly from those used for the paper.

## Instructions

There are a number of steps in this paper, in approximately this order:

1. Input preparation - Oxford Nanopore 1D reads
    - Basecalling
    - Filtering
2. Input preparation - Illumina MiSeq 300bp paired-end
    - Trimming / pairing

3. Sample ID via BLAST
    - BLAST database preparation
    - Pairwise alignment with BLASTN
    - Parse BLAST output
    - Analyse BLAST performance

4. Genomics applications
    - Read matching / mapping with BWA and LAAST
    - De novo genome assembly with Abyss (MiSeq data only)
    - De novo genome assembly with Abyss (MiSeq + MinION data)
    - De novo genome assembly with CANU (MinION data only)
    - Genome assembly summary statistics with QUAST
    - Genome completeness estimation with CEGMA

5. (Phylogenomics analyses)[https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/in-silico-reference-genome-digest]:
    - Annotation of raw reads to putative CDS using SNAP
    - Identification of putative CDS (othologue identification) to: TAIR10 reference loci; Wickett *et al.* (2014) phylogenomics CDS; Korf (2004) CEGMA genes (BLAST, Python, perl)
    - Multiple sequence alignment (MUSCLE, Trimal)
    - Phylogenetic inference (Star-BEAST, RAxML, TreeAnnotator)

6. Other analyses to investigate robustness of BLAST-based ID-by-sequencing
    - (*in silico* reference genome digest)[https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/in-silico-reference-genome-digest] to investigate likely pairwise performance on incomplete or poor-contiguity reference genomes
    - Jacknife sampling of MinION reads to investigate relationship between sequencing yield and identification confidence

## Dependencies

There are a number of external bioinformatics / data analysis tools required for this code to run. Which (with tested versions and citations) are laid out in [Dependencies.md](Dependencies.md)

Note that many of these tools (e.g. raxmlHPC; QUAST; CEGMA for instance) may need additional steps to integrate them into your ```$PATH```.
