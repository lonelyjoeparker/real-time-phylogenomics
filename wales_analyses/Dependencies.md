# Known dependencies

There are a number of external bioinformatics tools these methods depend on. Versions used are listed below. Full platform-specific installation instructions for each are beyond the scope of this documentation.

Note that standard computaional science tools (Python 2.7+; perl 5; POSIX-compatible; bash shell; R; etc) are also assumed.

## Lab / fieldwork dependencies

 - DNA extraction: Qiagen DNeasy plant mini prep kit (Qiagen)
 - DNA quantification: Quantus fluorimeter (Promega)
 - DNA sequencing, field-based
     - Hardware: MinION Mk1b (Oxford Nanopore Technologies)
     - Flowcells: R7.3 and R9 1D MinION Flowcells (ONT)
     - Library kit: Nanopore RAD-001 (ONT)
 - DNA sequencing, lab-based
     - Hardware: MiSeq (Illumina)
     - Library kit: NEBNext Ultra II sequencing library kit (paired-end, 300bp; New England Biolabs)



## Bioinformatics tools

Software | VERSION | Citation / purpose
-------- | ------- | ------------------
Trimmomatic | 0.32 | Bolger et al. (2014) / Illumina read adapter trimming
ncbi-blast suite (BLASTN, BLASTP, makeblastdb) | 2.2.31 | Camacho et al. (2008) /  Sequence matching/mapping
nanocall | 0.6.13 | [https://github.com/mateidavid/nanocall] / Nanopore basecalling
BWA | 0.7.12-r1039 | Li & Durbin, (2009) / read mapping
LAST | 581 | Kielbasa et al. (2011) / read mapping/alignment
ROCR | NA | Sing et al. (2005) / Classifier performance evaluation
ABYSS | 1.9.0 | Simpson et al. (2009) / de novo genome assembly (short-read)
HybridSPAdes | 3.5.0 | Antipov et al. (2016) / de novo genome assembly (hybrid short- and long-read)
Canu | 1.3 | Koren at al. (2016) / de novo genome assembly (long-read)
Quast | 4.3 | Gurevich et al. (2013) / Genome assembly statistics
CEGMA | 2.5 | Parra et al. (2007) / Genome completeness estimation
SNAP | NA | Korf (2004) / Gene annotation
muscle | 3.8.31 | Edgar (2004) / Multiple sequence alignment
Trimal | 1.4.rev15 | Capella-Gutierrez et al. (2009) / Alignment trimming and manipulation
raxmlHPC | 8.2.4 | Stamatakis (2014) / Phylogeny inference via maximum-likelihood
(Star)BEAST | 2.4.4 | Bouckaert et al. (2014) / Beyesian phylogeny inference using the multispecies coalescent
Tracer | 1.5 | Drummond et al. (2012) / Inspect performance of MCMC inference (BEAST)
TreeAnnotator | 1.7.4 | Drummond et al. (2012) / Compile maximum clade-credibility trees from multiple input phylogenies


## References

Antipov D, Korobeynikov A, McLean JS, Pevzner PA. (2016) hybridSPAdes: an algorithm for hybrid
assembly of short and long reads. *Bioinformatics* **32**(7):1009-15. doi: 10.1093/bioinformatics/btv688.

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina
Sequence Data. *Bioinformatics* **30**(15):2114-2120.
Bouckaert, R., Heled, J., Kühnert, D., Vaughan, T., Wu, C-H., Xie, D., Suchard, MA., Rambaut, A., &

Drummond, A. J. (2014). BEAST 2: A software platform for bayesian evolutionary analysis. *PLoS
Computational Biology*, **10**(4):e1003537.

Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008)
BLAST+: architecture and applications. *BMC Bioinformatics* **10**:421.

Drummond, A.J., Suchard, M.A., Xie, D. & Rambaut, A. (2012) Bayesian phylogenetics with BEAUti
and the BEAST 1.7 *Molecular Biology And Evolution* **29**:1969-1973.

Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.
*Nucleic Acids Res* **32**(5):1792-97.

Gureyvich, A., Saveliev, V., Vyahhi, N. & Tesler, G. (2013) QUAST: quality assessment tool for
genome assemblies. *Bioinformatics* **29**(8):1072-1075.

Heled, J. and Drummond, A.J. (2010) Bayesian inference of species trees from multilocus data. *Mol.
Biol. Evol.* **27**(3):570-580.

Kiełbasa SM, Wan R, Sato K, Horton P, Frith MC. (2011) Adaptive seeds tame genomic sequence
comparison. *Genome Res.* **21**(3):487-93.

Koren, S., Walenz, B.P., Berlin, K., Miller, J.R. & Phillippy, A.M. (2016) Canu: scalable and accurate
long-read assembly via adaptive k-mer weighting and repeat separation.
http://dx.doi.org/10.1101/071282.

Korf I. (2004) Gene finding in novel Genomes. *BMC Bioinformatics* **5**:59.

Li, H. & Durbin, R. (2009) Fast and accurate short read alignment with Burrows–Wheeler transform.
*Bioinformatics* **25**(14):1754–1760.

Parra, G., Bradnam, K. & Korf, I. (2007) CEGMA: a pipeline to accurately annotate core genes in
eukaryotic genomes. *Bioinformatics* **23**(9):1061-1067.

Simpson, Jared T., Kim Wong, Shaun D. Jackman, Jacqueline E. Schein, Steven JM Jones, and Inanc
Birol. (2009) ABySS: a parallel assembler for short read sequence data. *Genome research* **19**(6):1117-
1123.

Sing, T., Sander, O., Beerenwinkel, N. & Lengauer, T. (2005) ROCR: Visualizing classifier
performance in R. *Bioinformatics* **21**(20):3940-3941.

Stamatakis, A. (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large
phylogenies. *Bioinformatics* **30**(9):1312-1313.

Wickett, N.J., et al. (2014) Phylotranscriptomic analysis of the origin and early diversification of land
plants. *PNAS* **111**(45):E4859-4868.
