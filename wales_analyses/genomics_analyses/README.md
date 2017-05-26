# Parker et al (2017 in prep) RTnS and HTS genomics analyses

## Introduction

This repository contains code, instructions and listed dependencies needed to re-run the analyses presented in our field-sequencing paper (preprint: http://biorxiv.org/content/early/2017/03/14/107656).

We are in the process of fully documenting this code to make these scripts easier to run on your own system, please note however that producing a multifunctional analysis tool is not our intention (rather enhancing the reproducibility of this paper). This work partly involves altering where certain executables and input files are located to make porting the scripts simpler; as a result later [commits](https://github.com/lonelyjoeparker/real-time-phylogenomics/commits/master) may alter the various component scripts slightly from those used for the paper.

For the other analyses in the paper see [../README.md](../README.md)

## Instructions

* Read matching / mapping with BWA and LAAST
* *de novo* genome assembly with Abyss (MiSeq data only)
* *de novo* genome assembly with Abyss (MiSeq + MinION data)
* *de novo* genome assembly with CANU (MinION data only)
* Calculate genome assembly summary statistics with QUAST
* Genome completeness estimation with CEGMA

## Dependencies

There are a number of external bioinformatics / data analysis tools required for this code to run. Which (with tested versions and citations) are laid out in [Dependencies.md](../Dependencies.md)

Note that many of these tools (e.g. raxmlHPC; QUAST; CEGMA for instance) may need additional steps to integrate them into your ```$PATH```.
