# notes on prerequisites (see also Evernote https://www.evernote.com/shard/s571/nl/2147483647/10d95461-1265-4745-9d32-9afe670ca3ea/):

## extract info from basecalled fast5 files:
#	python2.7 poretools fasta <dir> > output.fasta
#	python2.7 poretools index <dir> > output.index
#
# .fasta for R7 and R9 thaliana reads concatenated and cleaned of whitespace:
# ~/Downloads/wales_nelumbolab_output/all_R7R9_thaliana.cleaned.fasta
# 
## run SNAP to predict genes, using 'thale' .gff gene models 
#	nb. release date for SNAP is 2006, not sure how old those gene models are
#	 ./snap HMM/thale ~/Downloads/wales_nelumbolab_output/all_R7R9_thaliana.fasta > ~/Downloads/wales_nelumbolab_output/all_R7R9_thaliana.fasta.snap.prediction 2>/dev/null
#
## export nucleotide seqs (+1000bp flankers) from prediction:
#	/Applications/Phylogenetics/snap/fathom all_R7R9_thaliana.fasta.ann all_R7R9_thaliana.cleaned.fasta -export 1000  2>/dev/null
#
## verifying the predictions are sensible with blastn:
#	Make a blast DB for TAIR10 genes:
#	makeblastdb -dbtype nucl -in TAIR10_seq_20101214_updated.txt -out TAIR10_genes_nt
#	BLASTN the predicted genes-from-reads against TAIR10:
#	blastn -db TAIR10_genes_nt -query export.dna -num_threads 8 -max_target_seqs 1 -max_hsps 1 -outfmt "6 sacc qacc length evalue" > export.blast.results
#	total hits: 28745
#	hits > 500: 22833
#	hits > 300: 25432
#	hits > 100: 27639
#
## then parse the poretools index / fasta to get times that can be linked to read UID (hash)
# 	python parse_poretools_index.py --input R7_A_thaliana.poretools.index --fasta R7_A_thaliana.poretools.fasta > all_R7R9_thaliana.times.tdf
#	python parse_poretools_index.py --input _Volumes_SCI-FEST-A_R9_A_thaliana.poretools.index --fasta _Volumes_SCI-FEST-A_R9_A_thaliana.poretools.fasta >> all_R7R9_thaliana.times.tdf 
#
#	output: all_R7R9_thaliana.times.tdf
#
## grep the gene names from the TAIR10 .gff into ~/Downloads/wales_nelumbolab_output/TAIR10_GFF3_names.txt
#
## then produce the genetimes data:
#	python parse_SNAP_genes_to_TAIR.py --blast export.blast.results --times all_R7R9_thaliana.times.tdf --genes TAIR10_GFF3_names.txt >all_R7R9_thaliana.genetimes.rdata
#
#	output: all_R7R9_thaliana.genetimes.rdata
#
# /prerequisite

# read input
#diff gts = read.table("all_R7R9_thaliana.genetimes.rdata",sep="\t",header=T)
gts = read.table("~/Downloads/phylogenome_wales/SNAP-predicted/all_R7R9_thaliana.genetimes.filtered.nt.rdata",sep="\t",header=T)
ats = read.table("~/Downloads/phylogenome_wales/SNAP-predicted/all_R7R9_thaliana.genetimes.filtered.aa.rdata",sep="\t",header=T)
str(gts)
str(ats)
# loopthrough and plot
i=1
t0 = min(gts$times)
cum_read_counts = c()
cum_gene_counts = c()
cum_read_times = c()
for (t in seq(1,(max(gts$times)-t0),length.out=20)){
	dt = t0 + t
	cum_read_counts[i] = length(gts$times[gts$times<dt])	# number of reads younger than dt
	cum_gene_counts[i] = length(levels(as.factor(as.character(gts[gts$times<dt,]$TAIR10_gene)))) # number of genes (excl. duplicates) hit by reads younger than dt
	cum_read_times[i] = dt
	i = i+1
}

#plot read cumulative count
plot(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts,main='Genes predicted directly from single Nanopore reads',xlab='Sequencing time',ylab='Cumulative num. reads matched against genes (non-unique)',sub='Start: 18:00:04 BST 18 May 2016')
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts)
polygon(
	c(
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))
polygon(
	c(
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))	

#plot gene cumualtive count
plot(
	as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,
	main='Genes predicted directly from single Nanopore reads',
	xlab='Sequencing time',
	ylab='Cumulative num. genes hit by one or more reads',
	sub='Start: 18:00:04 BST 18 May 2016')
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts)
polygon(
	c(
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))
polygon(
	c(
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))	

#plot both, why not
plot(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts,main='Genes predicted directly from single Nanopore reads',xlab='Sequencing time',ylab='Cumulative num. reads matched against genes (non-unique)',sub='Start: 18:00:04 BST 18 May 2016')
polygon(
	c(
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+19600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+16000,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))
polygon(
	c(
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+48400,origin="1970-01-01")
	),
	c(0,0,18000,18000),
	col=rgb(0.5,0.1,0.1,0.1))	
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts)
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,col="red")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,col="red")
legend(t0,17000,col="black","Cumulative num. predicted\ngenes (inc. duplicates)",pch=1)
legend(t0,15000,col="red"  ,"Cumulative num. TAIR10\ngenes hit by >1 reads",pch=1)
legend(t0,13000,fill=rgb(0.5,0.1,0.1,0.1),"MinION paused in transit")


cumulative=data.frame(times=cum_read_times,human_time=as.POSIXct(cumulative$times,origin="1970-01-01"),reads=cum_read_counts,genes=cum_gene_counts)
> cumulative
        times          human_time reads genes
1  1463590805 2016-05-18 18:00:05     1     1
2  1463594842 2016-05-18 19:07:21  3148  2155
3  1463598878 2016-05-18 20:14:38  5551  3559
4  1463602915 2016-05-18 21:21:55  7625  4667
5  1463606952 2016-05-18 22:29:11  7779  4754
6  1463610989 2016-05-18 23:36:28  7779  4754
7  1463615025 2016-05-19 00:43:45  9145  5435
8  1463619062 2016-05-19 01:51:02 11706  6694
9  1463623099 2016-05-19 02:58:18 13443  7482
10 1463627136 2016-05-19 04:05:35 14521  7965
11 1463631172 2016-05-19 05:12:52 15251  8267
12 1463635209 2016-05-19 06:20:09 15838  8509
13 1463639246 2016-05-19 07:27:25 16019  8601
14 1463643283 2016-05-19 08:34:42 16019  8601
15 1463647319 2016-05-19 09:41:59 16019  8601
16 1463651356 2016-05-19 10:49:16 16019  8601
17 1463655393 2016-05-19 11:56:32 16019  8601
18 1463659430 2016-05-19 13:03:49 16019  8601
19 1463663466 2016-05-19 14:11:06 17402  9194
20 1463667503 2016-05-19 15:18:23 18097  9501

## repeat with filtering

#measure	length	pident	evalue         
#Min.: 	28	70.36	0.000e+00  
#1st Qu.:	570	78.68	0.000e+00  
#Median:	980	82.93	0.000e+00  
#Mean:	1070	82.93	3.706e-08  
#3rd Qu.:	1482	87.00	0.000e+00  
#Max.:	6826	100.00	5.000e-05 
 
#all:	length(gts[(gts$pident>0)&(gts$length>0),1])	18098	
#filter:	length(gts[(gts$pident>78.68)&(gts$length>570),1])	10102
gts_filter=gts[(gts$pident>78.68)&(gts$length>570),]
# loopthrough and plot, this time with 100 time-increments for better granularity
i=1
t0 = min(gts_filter$times)
cum_read_counts = c()
cum_gene_counts = c()
cum_read_times = c()
for (t in seq(1,(max(gts_filter$times)-t0),length.out=100)){
	dt = t0 + t
	cum_read_counts[i] = length(gts_filter$times[gts_filter$times<dt])	# number of reads younger than dt
	cum_gene_counts[i] = length(levels(as.factor(as.character(gts_filter[gts_filter$times<dt,]$TAIR10_gene)))) # number of genes (excl. duplicates) hit by reads younger than dt
	cum_read_times[i] = dt
	i = i+1
}

### now replot this as above...

### now the AA data:

## bit more plotting / filtering
ats$idents = ats$length * (ats$pident / 100) 		# number of matches
ats$mismatch = ats$length - ats$idents - ats$gaps 	# number of mismatches
ats$gap_fraction = na.omit(ats$gaps/ats$mismatch)			# normalize num gaps by num total mismatches, e.g. find alignments that are more gappy than expected
ats$gap_open_rate  = ats$num_gap_open/ats$length	# normalize gap openings by total length e.g. find alignments that have fewer gap opening events (e.g. long contiguous chunks missing which look like dropped exons) than expected
# plot 'em
plot(gap_fraction~gap_open_rate,data=ats,col="light grey",ylim=c(-100,100)
points(gap_fraction~gap_open_rate,data=ats[ats$pident>30,],col="yellow")
points(gap_fraction~gap_open_rate,data=ats[ats$pident>35,],col="orange")
points(gap_fraction~gap_open_rate,data=ats[ats$pident>40,],col="red")
points(gap_fraction~gap_open_rate,data=ats[ats$pident>45,],col="purple")

#that doesn't show any *particularly* clear patterns, so just pick the top 5% by length, gap fraction and gap opening rate
length(ats[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1),1]) # 5351
ats_filtered_top_5pc = ats[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1),]
# loopthrough and plot
i=1
t0 = min(ats$times[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1)])
cum_read_counts_aa = c()
cum_gene_counts_aa = c()
cum_read_times_aa = c()
#for (t in seq(1,(max(ats$times[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1)])-t0),length.out=100)){
for (t in seq(1,(max(gts_filter$times)-t0),length.out=100)){		#use the same times as nt filtered data (gts_filtered) so they can be plotted on the same axes
	dt = t0 + t
	cum_read_counts_aa[i] = length(ats$times[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1)&(ats$times<dt)])	# number of reads younger than dt
	cum_gene_counts_aa[i] = length(levels(as.factor(as.character(ats[(ats$length>100)&(ats$gap_fraction>0.16)&(ats$gap_open_rate<0.1)&(ats$times<dt),]$TAIR10_gene)))) # number of genes (excl. duplicates) hit by reads younger than dt
	cum_read_times_aa[i] = dt
	i = i+1
}


#plot both, why not
plot(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts,main='Genes predicted directly from single Nanopore reads',xlab='Sequencing time',ylab='Cumulative num. coding sequences predicted from reads',sub='Start: 18:00:04 BST 18 May 2016')
polygon(
	c(
		as.POSIXct(min(cum_read_times)+12000,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+21600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+21600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+12000,origin="1970-01-01")
	),
	c(0,0,5000,5000),
	col=rgb(0.5,0.1,0.1,0.05))
polygon(
	c(
		as.POSIXct(min(cum_read_times)+46400,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+46400,origin="1970-01-01")
	),
	c(0,0,10000,10000),
	col=rgb(0.5,0.1,0.1,0.05))	
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts,lty=3)
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,col="black")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,col="black")
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts_aa,col="red")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts_aa,col="red",lty=3)
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts_aa,col="red")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts_aa,col="red")
legend(t0,10000,col="black","Cumulative num. predicted\ngenes (inc. duplicates)",pch=1,lty=3)
legend(t0,9000,col="black"  ,"Cumulative num. TAIR10\ngenes hit by >1 reads",pch=1,lty=1)
legend(t0,8000,col="red","Cumulative num. predicted\nproteins sequences (inc. duplicates)",pch=1,lty=3)
legend(t0,7000,col="red"  ,"Cumulative num. TAIR10\nproteins hit by >1 reads",pch=1,lty=1)
legend(t0,6000,fill=rgb(0.5,0.1,0.1,0.05),"MinION paused in transit (approx)")



### final plot for Figure 2 in Papadopulos et al
#plot both, why not
plot(
	as.POSIXct(cum_read_times,origin="1970-01-01"),
	cum_read_counts,	
	col="white",
	mar=c(5,5,5,5),
	cex.axis=2,
	cex.lab=1.4,
	bty="n",
	xlab='Time (UTC+1)',
	ylab='Predicted genes',
	sub='18 May 2016'
	)
polygon(
	c(
		as.POSIXct(min(cum_read_times)+12000,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+22600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+22600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+12000,origin="1970-01-01")
	),
	c(0,0,10101,10101),
	border=rgb(1,1,1,1.00),
	col=rgb(0.9,0.1,0.1,0.1))
polygon(
	c(
		as.POSIXct(min(cum_read_times)+46400,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+68600,origin="1970-01-01"),
		as.POSIXct(min(cum_read_times)+46400,origin="1970-01-01")
	),
	c(0,0,10101,10101),
	border=rgb(1,1,1,1.00),
	col=rgb(0.9,0.1,0.1,0.1))	
// don't plot points 
//points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,col="white")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"), cum_read_counts,lty=3,lwd=10)
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts,lwd=10)
/* don't plot AAs for publication figure
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts_aa,col="red")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_read_counts_aa,col="red",lty=3)
points(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts_aa,col="red")
lines(as.POSIXct(cum_read_times,origin="1970-01-01"),cum_gene_counts_aa,col="red")
 */
/* or legend text 
legend(t0,10000,col="black","Cumulative num. predicted\ngenes (inc. duplicates)",pch=1,lty=3)
legend(t0,9000,col="black"  ,"Cumulative num. TAIR10\ngenes hit by >1 reads",pch=1,lty=1)
legend(t0,8000,col="red","Cumulative num. predicted\nproteins sequences (inc. duplicates)",pch=1,lty=3)
legend(t0,7000,col="red"  ,"Cumulative num. TAIR10\nproteins hit by >1 reads",pch=1,lty=1)
legend(t0,6000,fill=rgb(0.5,0.1,0.1,0.05),"MinION paused in transit (approx)")
 */