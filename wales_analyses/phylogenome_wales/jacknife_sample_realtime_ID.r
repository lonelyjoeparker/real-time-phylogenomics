# assumes ~/Downloads/phylogenome_wales/BLAST_ROC-final-mod.r has 
# been loaded to make the pore and miseq data available
# e.g source('~/Downloads/phylogenome_wales/BLAST_ROC-final-mod.r')

# A script to jacknife sample from A.thaliana read data
# The aim is to simulate real-time ID by BLAST and determine
# likely error rates (e.g. 'how many reads are needed for 
# ID'-type questions.)

library(ggplot2)

# a script to calculate the quantities
calc_ID_stats <- function(subsample){
	# clean NAs out
	subsample=na.omit(subsample)
	N = length(subsample)
	# assuming a subsample of the data has been taken
	# calculate some of the likely ID stats
	# and return them in a vector
	# cutoff-based stats
	num_cutoff_0   = length(subsample[subsample>0])
	num_cutoff_1   = length(subsample[subsample>1])
	num_cutoff_10  = length(subsample[subsample>10])
	num_cutoff_100 = length(subsample[subsample>100])
	# score the cutoffs
	score_cutoff_0 	 = num_cutoff_0 / N
	score_cutoff_1 	 = num_cutoff_1 / N
	score_cutoff_10  = num_cutoff_10 / N
	score_cutoff_100 = num_cutoff_100 / N
	# aggregate stats
	mean_diff=mean(subsample)
	sum_diff=sum(subsample)
	#return the output
	return(c(score_cutoff_0,score_cutoff_1,score_cutoff_10,score_cutoff_100,mean_diff,sum_diff))
}

# a script to repeat the stat on some data 
# returns a row for a data frame
do_replication <- function(vector,label,subsample_count,replicates){
	# vector 	: the observations under test
	# label 	: some useful label of some kind
	# subsample_count : number of elements in vector that should be subsampled with each replicate
	# replicate	: the number of replicates overall
	#
	# set up an empty frame
	empty=data.frame(cutoff_0=NA,cutoff_1=NA,cutoff_10=NA,cutoff_100=NA,mu=NA,agg=NA)
	for(i in 1:replicates){
		empty[i,] = calc_ID_stats(sample(vector,subsample_count))
	}
	means=c()
	vars=c()
	# summary stats
	for(i in 1:6){
		means[i]= mean(empty[,i])
		vars[i]=  sd(empty[,i])
	}
	# set up the return vector
	returns=c(label,subsample_count,replicates,means,vars)
	#print(empty)
	#print(returns)
	return(returns)
}

# a script to iterate through and get results(!)
ID_frame_tp=data.frame(num_reads=NA,labal=NA,ss_count=NA,replicates=NA,mu_cutoff_0=NA,mu_cutoff_1=NA,mu_cutoff_10=NA,mu_cutoff_100=NA,mu_readsAvg=NA,mu_reads_agg=NA,var_cutoff_0=NA,var_cutoff_1=NA,var_cutoff_10=NA,var_cutoff_100=NA,var_readsAvg=NA,var_reads_agg=NA)
subsamples_to_make=c(1:9,seq(10,90,by=10),seq(100,900,by=100),seq(1000,10000,by=1000))
for(i in 1:length(subsamples_to_make)){
	ID_frame_tp[i,] = c(subsamples_to_make[i],do_replication(pores_all$diff_nident[pores_all$TP=='thaliana'],'thaliana',subsamples_to_make[i],1000))
}

four_colours=rainbow(4)
#by now have:
#> str(ID_frame_tp)
#'data.frame':	130 obs. of  16 variables:
# $ num_reads     : chr  "1" "2" "3" "4" ...
# $ labal         : chr  "thaliana" "thaliana" "thaliana" "thaliana" ...
# $ ss_count      : chr  "1" "2" "3" "4" ...
# $ replicates    : chr  "30" "30" "30" "30" ...
# $ mu_cutoff_0   : chr  "1" "0.966666666666667" "0.933333333333333" "0.966666666666667" ...
# $ mu_cutoff_1   : chr  "0.966666666666667" "0.966666666666667" "0.922222222222222" "0.95" ...
# $ mu_cutoff_10  : chr  "0.933333333333333" "0.866666666666667" "0.811111111111111" "0.816666666666667" ...
# $ mu_cutoff_100 : chr  "0.766666666666667" "0.65" "0.6" "0.633333333333333" ...
# $ mu_readsAvg   : chr  "1080.96192633333" "760.679360833333" "613.454452777778" "767.914227166667" ...
# $ mu_reads_agg  : chr  "1080.96192633333" "1521.35872166667" "1840.36335833333" "3071.65690866667" ...
# $ var_cutoff_0  : chr  "0" "0.0160919540229885" "0.0183908045977012" "0.0117816091954023" ...
# $ var_cutoff_1  : chr  "0.0333333333333333" "0.0160919540229885" "0.0282247765006386" "0.0146551724137931" ...
# $ var_cutoff_10 : chr  "0.064367816091954" "0.0505747126436782" "0.0588761174968072" "0.0385057471264368" ...
# $ var_cutoff_100: chr  "0.185057471264368" "0.123275862068966" "0.0873563218390805" "0.0419540229885058" ...
# $ var_readsAvg  : chr  "1213893.96121447" "757306.283869482" "226621.004043959" "372434.481166005" ...
# $ var_reads_agg : chr  "1213893.96121447" "3029225.13547793" "2039589.03639563" "5958951.69865608" ...

plot(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_0,col=four_colours[1],pch=2,ylim=c(0,1),xlab="number of reads (log10)",ylab="correct guesses")
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_1,col=four_colours[2],pch=2)
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_10,col=four_colours[3],pch=2)
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_100,col=four_colours[4],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_0,col=four_colours[1],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_1,col=four_colours[2],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_10,col=four_colours[3],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_cutoff_100,col=four_colours[4],pch=2)
legend(2,0.4,"cutoff>0",col=four_colours[1],pch=2,bty='n')
legend(2,0.3,"cutoff>1",col=four_colours[2],pch=2,bty='n')
legend(2,0.2,"cutoff>10",col=four_colours[3],pch=2,bty='n')
legend(2,0.1,"cutoff>100",col=four_colours[4],pch=2,bty='n')

plot(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_0,col=four_colours[1],pch=2,ylim=c(0,0.5),xlab="number of reads (log10)",ylab="correct guess variance")
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_1,col=four_colours[2],pch=2)
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_10,col=four_colours[3],pch=2)
points(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_100,col=four_colours[4],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_0,col=four_colours[1],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_1,col=four_colours[2],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_10,col=four_colours[3],pch=2)
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$var_cutoff_100,col=four_colours[4],pch=2)
legend(2,0.4,"cutoff>0",col=four_colours[1],pch=2,bty='n')
legend(2,0.3,"cutoff>1",col=four_colours[2],pch=2,bty='n')
legend(2,0.2,"cutoff>10",col=four_colours[3],pch=2,bty='n')
legend(2,0.1,"cutoff>100",col=four_colours[4],pch=2,bty='n')


plot(log10(as.numeric(ID_frame_tp$ss_count)),log10(as.numeric(ID_frame_tp$var_readsAvg)),xlab="reads (log10)",ylab="variance (mean read num. idents. diff)")
plot(log10(as.numeric(ID_frame_tp$ss_count)),log10(as.numeric(ID_frame_tp$var_reads_agg)),xlab="num. reads (log10)",ylab="aggregate num. idents. difference (log10)")
plot(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_readsAvg,ylab="mean read num. idents. diff",xlab="num. reads (log10)")
lines(log10(as.numeric(ID_frame_tp$ss_count)),ID_frame_tp$mu_readsAvg,ylab="mean read num. idents. diff",xlab="num. reads (log10)")


# oddly these vars all now characters. change to numerics
ID_frame_clone=ID_frame_tp
for(i in c(1,3:16)){
	ID_frame_clone[,i]=as.numeric(ID_frame_clone[,i])
}
ID_frame_tp=ID_frame_clone[order(ID_frame_clone$ss_count),]

# to make ggplot work properly (with groups):
ID_frame_tp$cutoff=0
ID_frame_tp$cutoff_mu=ID_frame_tp$mu_cutoff_0
ID_frame_tp$cutoff_var=ID_frame_tp$var_cutoff_0
cutoff_1=ID_frame_tp
cutoff_1$cutoff=1
cutoff_1$cutoff_mu=cutoff_1$mu_cutoff_1
cutoff_1$cutoff_var=cutoff_1$var_cutoff_1
cutoff_10=ID_frame_tp
cutoff_10$cutoff=10
cutoff_10$cutoff_mu=cutoff_10$mu_cutoff_10
cutoff_10$cutoff_var=cutoff_10$var_cutoff_10
cutoff_100=ID_frame_tp
cutoff_100$cutoff=100
cutoff_100$cutoff_mu=cutoff_100$mu_cutoff_100
cutoff_100$cutoff_var=cutoff_100$var_cutoff_100

#bind them to a single table
cutoff_with_groups=rbind(ID_frame_tp[,c(1:4,17:19)],cutoff_1[,c(1:4,17:19)],cutoff_10[,c(1:4,17:19)],cutoff_100[,c(1:4,17:19)])

# convert cutoff to a factor
cutoff_with_groups$cutoff=as.factor(cutoff_with_groups$cutoff)

#plot cutoffs
cutoffPlotLog=ggplot(data=cutoff_with_groups,aes(x=ss_count,y=cutoff_mu,col=cutoff,group=cutoff)) +
	geom_point(size=4) +
	geom_errorbar(
		aes(ymin=cutoff_mu-cutoff_var, ymax=cutoff_mu+cutoff_var),
		size=1,
		width=0
	) +
	xlab("Number of reads (log10)") + 
	ylab(expression(paste("Correct ",italic("A. thaliana")," read ID fraction"))) +
	ggtitle("Simulated ID precision") +
	geom_line(linetype='dashed',size=1.1) +
	theme_bw() + theme(
		legend.key=element_blank(),
		axis.text =element_text(size=15),
		axis.title=element_text(size=13),
		panel.grid.major.y  = element_line(colour="grey",linetype="solid"  ,size=0.5), 
		panel.grid.major.x  = element_line(colour="grey",linetype="solid"  ,size=1), 
		panel.grid.minor.x = element_line(colour="grey"     ,linetype="dashed" ,size=1)
	) +
	scale_color_discrete("Cutoff threshold") +
	coord_trans(x="log10") +
	scale_x_continuous(
		breaks = trans_breaks("log10", function(x) 10^x),
		labels = trans_format("log10", math_format(10^.x))
	)
cutoffPlotLog 

diffPlot=ggplot(data=ID_frame_tp,aes(x=ss_count,y=mu_readsAvg)) +
	geom_point(size=4) +
	geom_errorbar(
		aes(ymin=mu_readsAvg-var_readsAvg, ymax=mu_readsAvg+var_readsAvg),
		size=1,
		width=0.0
	) +
	theme_bw() + theme(
		legend.key=element_blank(),
		axis.text =element_text(size=15),
		axis.title=element_text(size=13),
		panel.grid.major.y  = element_line(colour="grey",linetype="solid"  ,size=0.5), 
		panel.grid.major.x  = element_line(colour="grey",linetype="solid"  ,size=1), 
		panel.grid.minor.x = element_line(colour="grey"     ,linetype="dashed" ,size=1)
	) +
	xlab("Number of reads (log10)") +
	ylab("Mean BLAST num. idents bias") +
	ggtitle("Simulated ID progress") +
	#geom_smooth(linetype='dotted',method="loess") +
	geom_line(linetype='dashed',size=1.1) +
	coord_trans(x="log10") +
	scale_x_continuous(
		breaks = trans_breaks("log10", function(x) 10^x),
		labels = trans_format("log10", math_format(10^.x))
	)
diffPlot 

#ggplot.geom_errorbar() isn't going to like -ve values when we come to log10 plot
#so floor zero
agg_reads_error_y_low=ID_frame_tp$mu_reads_agg - ID_frame_tp$var_reads_agg
#> agg_reads_error_y_low
# [1]    -303.22482      32.92865     427.52462     846.12843    1453.02046
# [6]    2060.48301    2601.60487    3128.99363    3711.82832    4177.02498
#[11]   10454.53865   16757.13514   23730.05952   30779.55138   37093.13080
#[16]   44070.16641   50615.15853   58398.95760   65434.97856  136334.72270
#[21]  207209.38044  282231.46068  353367.21398  426886.49625  502284.30937
#[26]  575593.36228  646459.55693  722922.83800 1466891.98225 2210146.49343
#[31] 2956383.14375 3705818.91622 4457565.33330 5211605.05791 5956099.13163
#[36] 6709209.61479 7454891.12649
agg_reads_error_y_low[1]=1

# OK now should be good to plot
aggPlot=ggplot(data=ID_frame_tp,aes(x=ss_count,y=mu_reads_agg)) +
	theme_bw() + theme(
		legend.key=element_blank(),
		axis.text =element_text(size=15),
		axis.title=element_text(size=13),
		panel.grid.major.y = element_line(colour="grey",linetype="solid"  ,size=0.7), 
		panel.grid.minor.y = element_line(colour="grey",linetype="dashed" ,size=0.3), 
		panel.grid.major.x = element_line(colour="grey",linetype="solid"  ,size=0.5), 
		panel.grid.minor.x = element_line(colour="grey",linetype="dashed" ,size=0.1)
	) +
	geom_point(size=4) +
	geom_errorbar(			#NB using error low values coerced ≥1
		aes(ymin=agg_reads_error_y_low, ymax=mu_reads_agg+var_reads_agg),
		size=1,
		width=(1:37)*0.25	# bodge the error bars' horizontal width
	) + 
	xlab("Number of reads (log10)") +
	ylab("Cumulative BLAST num. idents bias (log10)") +
	ggtitle("Simulated ID progress") +
	#geom_smooth(linetype='dotted',method="loess") +
	geom_line(linetype='dashed',size=1.1) +
	coord_trans(x="log10",y="log10") +
	scale_x_continuous(
		breaks = trans_breaks("log10", function(x) 10^x),
		labels = trans_format("log10", math_format(10^.x))
	) +
	scale_y_continuous(
		breaks = trans_breaks("log10", function(x) 10^x),
		labels = trans_format("log10", math_format(10^.x))
	)
aggPlot 
