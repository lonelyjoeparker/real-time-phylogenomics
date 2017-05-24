# quick R script to analyse the results of simulated genome fragmentation
# on pairwise sample ID by BLASTN
#
# Joe Parker / @lonelyjoeparker, RBG Kew, 2017

# hardcode original sims file
parsed_sims_file = "~/Documents/all_work/programming/repo-git/real-time-phylogenomics/wales_analyses/in-silico-reference-genome-digest/simulate_genome_fragmentation.tdf"

sims = read.table(parsed_sims_file,sep="\t",header=T)

# subset the data by treatment
pores_thaliana = subset(sims, sample.assumed=='thaliana' & platform=='pores')
pores_lyrata   = subset(sims, sample.assumed=='lyrata'   & platform=='pores')
miseq_thaliana = subset(sims, sample.assumed=='thaliana' & platform=='miseq')
miseq_lyrata   = subset(sims, sample.assumed=='lyrata'   & platform=='miseq')

# plot N50 vs ID score (mean length bias) for treatments
# nb. N50 on log(10) scale
quartz()
plot(log10(sims$N50),sims$length.bias.mean,main="Simulated and empirical BLAST ID values by N50",xlab="N50 (approx)",ylab="Mean length bias",pch=3)
lines(log10(pores_thaliana$N50),pores_thaliana$length.bias.mean,col="red",lwd=3)
lines(log10(pores_lyrata$N50),  pores_lyrata$length.bias.mean,  lty=2,    lwd=3)
lines(log10(miseq_thaliana$N50),miseq_thaliana$length.bias.mean,col="red",lwd=3)
lines(log10(miseq_lyrata$N50),  miseq_lyrata$length.bias.mean,  lty=2,    lwd=3)
legend(5.5,1200,"MinION (A.thaliana)",col="red",lwd=3)
legend(5.5, 900,"MinION (A.lyrata)",lty=2    ,lwd=3)
legend(5.5, 350, "MiSeq (A.thaliana)",col="red",lwd=3)
legend(5.5, 250, "MiSeq (A.lyrata)",lty=2    ,lwd=3)

# plot N50 vs ID score (mean num. identities bias) for treatments
# nb. N50 on log(10) scale
quartz()
plot(log10(sims$N50),sims$idents.bias.mean,main="Simulated and empirical BLAST ID values by N50",xlab="N50 (approx)",ylab="Mean num. identities bias",pch=3)
lines(log10(pores_thaliana$N50),pores_thaliana$idents.bias.mean,col="red",lwd=3)
lines(log10(pores_lyrata$N50),  pores_lyrata$idents.bias.mean,  lty=2,    lwd=3)
lines(log10(miseq_thaliana$N50),miseq_thaliana$idents.bias.mean,col="red",lwd=3)
lines(log10(miseq_lyrata$N50),  miseq_lyrata$idents.bias.mean,  lty=2,    lwd=3)
legend(5.5,1150,"MinION (A.thaliana)",col="red",lwd=3)
legend(5.5, 700,"MinION (A.lyrata)",lty=2    ,lwd=3)
legend(5.5, 350, "MiSeq (A.thaliana)",col="red",lwd=3)
legend(5.5, 250, "MiSeq (A.lyrata)",lty=2    ,lwd=3)