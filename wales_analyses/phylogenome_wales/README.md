# Parker et al (2017 in prep) RTnS phylogenomics analyses

## Introduction

This repository contains code, instructions and listed dependencies needed to re-run the analyses presented in our field-sequencing paper (preprint: http://biorxiv.org/content/early/2017/03/14/107656).

We are in the process of fully documenting this code to make these scripts easier to run on your own system, please note however that producing a multifunctional analysis tool is not our intention (rather enhancing the reproducibility of this paper). This work partly involves altering where certain executables and input files are located to make porting the scripts simpler; as a result later [commits](https://github.com/lonelyjoeparker/real-time-phylogenomics/commits/master) may alter the various component scripts slightly from those used for the paper.

For the other analyses in the paper see [../README.md](../README.md)

## Instructions

This folder contains code to perform the following steps in the phylogenomics analyses:

* Annotation of individual ONT reads using SNAP, and assignment to TAIR10 coding sequences (CDS) via BLAST
* Indexing of the TAIR10 loci to both Wickett *et al.* (2014) phylogenomics loci, and the CEGMA (Korf, 2004) loci via BLAST and database information
* Compilation of the above into multiple sequence alignments (using outgroup and backbone taxa) using MUSCLE and Trimal
* Phylogeny inference via RAxML (ML) and Star-BEAST (multispecies coalescent)
* Summary tree construction using TreeAnnotator.

## Dependencies

There are a number of external bioinformatics / data analysis tools required for this code to run. Which (with tested versions and citations) are laid out in [Dependencies.md](../Dependencies.md)

Note that many of these tools (e.g. raxmlHPC; QUAST; CEGMA for instance) may need additional steps to integrate them into your ```$PATH```.
