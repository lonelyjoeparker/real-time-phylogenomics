# Parker et al (2017 in prep) analyses

## Introduction

This repository contains code, instructions and listed dependencies needed to re-run the analyses presented in our field-sequencing paper (preprint: http://biorxiv.org/content/early/2017/03/14/107656).

We are in the process of fully documenting this code to make these scripts easier to run on your own system, please note however that producing a multifunctional analysis tool is not our intention (rather enhancing the reproducibility of this paper). This work partly involves altering where certain executables and input files are located to make porting the scripts simpler; as a result later [commits](https://github.com/lonelyjoeparker/real-time-phylogenomics/commits/master) may alter the various component scripts slightly from those used for the paper.

## Instructions

There are a number of steps in this paper, in approximately this order:

1. Input preparation - Oxford Nanopore 1D reads
    - [Basecalling](genomics_analyses/extract_fasta.sh) (Metrichor / EPI2ME and nanocall)
    - Filtering (BLASTN)
2. Input preparation - Illumina MiSeq 300bp paired-end
    - Trimming / pairing (TRimmomatic)

3. Sample ID via BLAST
    - BLAST database preparation
    ```bash
    # Arabidopsis_thaliana
    makeblastdb -in arabidopsis_thaliana_GCF_000001735.3_TAIR10_genomic.fna -dbtype nucl -out A_thaliana
    # A. lyrata
    makeblastdb -in Arabidopsis_lytata_ADBK01.1.fsa_nt -dbtype nucl
    # A. lyrata ssp. petraea
    makeblastdb -in Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt -dbtype nucl
    # A. lyrata ssp. petraea and A. lyrata combined
    makeblastdb -in Arabidopsis_lyrata_petraea_contigs_BASP01.1.fsa_nt  Arabidopsis_lytata_ADBK01.1.fsa_nt -dbtype nucl -out A_lyrata_ssp
    # Bacteriophage lambda
    makeblastdb -in lambda_phage_assembly_sequence.fasta -dbtype nucl -out lambda_genome
    ```
    - [Pairwise alignment with BLASTN](https://github.com/lonelyjoeparker/real-time-phylogenomics/blob/master/wales_analyses/phylogenome_wales/massive_BLAST_ROC-onlyCompareBlast.sh)
    - [Parse BLAST output](https://github.com/lonelyjoeparker/real-time-phylogenomics/blob/master/wales_analyses/phylogenome_wales/parseCompareBlastHits.pl)
    - [Analyse BLAST performance](https://github.com/lonelyjoeparker/real-time-phylogenomics/blob/master/wales_analyses/phylogenome_wales/BLAST_ROC-final-mod.r) via R ROCR package
    - See also [paper_final.sh](genomics_analyses/paper_final.sh) Badly named(!) since it is mainly input/filtering steps

4. [Genomics applications](https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/genomics_analyses)
    - [Read matching / mapping with BWA and LAST](https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/genomics_analyses/commands_mapping.sh)
    - De novo genome assembly with Abyss (MiSeq data only)
    - De novo genome assembly with Abyss (MiSeq + MinION data)
    - De novo genome assembly with CANU (MinION data only)
    - Genome assembly summary statistics with QUAST
    - Genome completeness estimation with CEGMA (see also [paper_final.sh](genomics_analyses/paper_final.sh) Badly named(!) since it is mainly input/filtering steps...)


5. [Phylogenomics analyses](https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/phylogenome_wales):
    - Annotation of raw reads to putative CDS using SNAP
    - Identification of putative CDS (othologue identification) to: TAIR10 reference loci; Wickett *et al.* (2014) phylogenomics CDS; Korf (2004) CEGMA genes (BLAST, Python, perl)
    - Multiple sequence alignment (MUSCLE, Trimal)
    - Phylogenetic inference (Star-BEAST, RAxML, TreeAnnotator)

6. Other analyses to investigate robustness of BLAST-based ID-by-sequencing
    - [*in silico* reference genome digest](https://github.com/lonelyjoeparker/real-time-phylogenomics/tree/master/wales_analyses/in-silico-reference-genome-digest) to investigate likely pairwise performance on incomplete or poor-contiguity reference genomes
    - [Jacknife sampling of MinION reads](https://github.com/lonelyjoeparker/real-time-phylogenomics/blob/master/wales_analyses/phylogenome_wales/jacknife_sample_realtime_ID.r) to investigate relationship between sequencing yield and identification confidence

## Dependencies

There are a number of external bioinformatics / data analysis tools required for this code to run. Which (with tested versions and citations) are laid out in [Dependencies.md](Dependencies.md)

Note that many of these tools (e.g. raxmlHPC; QUAST; CEGMA for instance) may need additional steps to integrate them into your ```$PATH```.
