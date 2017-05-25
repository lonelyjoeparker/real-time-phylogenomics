library(ROCR)

# read inputs

# ONT data, TP==A.thaliana
pores_tp_thaliana=read.table("/Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/BLAST_ROC/pairwise_type-pores_tp-thaliana_fp-lyrata.out",header=T,sep="\t")
# ONT data, TP==A.lyrata
pores_tp_lyrata=read.table("/Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/BLAST_ROC/pairwise_type-pores_tp-lyrata_fp-thaliana.out",header=T,sep="\t")
# miseq data, TP==A.thaliana
miseq_tp_thaliana=read.table("/Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/BLAST_ROC/pairwise_type-miseq_tp-thaliana_fp-lyrata.out",header=T,sep="\t")
# miseq data, TP==A.lyrata
miseq_tp_lyrata=read.table("/Volumes/LaCie/minION_analyses/wales_20160518/analyses/paper-final/BLAST_ROC/pairwise_type-miseq_tp-lyrata_fp-thaliana.out",header=T,sep="\t")

# number of queries (reads) with zero BLAST hits
zero_pores_tp_thaliana=58629
zero_pores_tp_lyrata  =25660
zero_miseq_tp_thaliana=185107
zero_miseq_tp_lyrata  =470270

# adding blank rows for reads that have had zero hits, laborious way...
# first work out how long the existing matrices are
##length_pores_tp_thaliana=length(pores_tp_thaliana[,1])
##pores_tp_thaliana[(length_pores_tp_thaliana+1):(length_pores_tp_thaliana+zero_pores_tp_thaliana),c(2,3,5,6)]=-1
##pores_tp_thaliana[(length_pores_tp_thaliana+1):(length_pores_tp_thaliana+zero_pores_tp_thaliana),c(4,7)]=999
##pores_tp_thaliana[(length_pores_tp_thaliana+1):(length_pores_tp_thaliana+zero_pores_tp_thaliana),8]='thaliana'
##length_pores_tp_lyrata=length(pores_tp_lyrata[,1])
##pores_tp_lyrata[(length_pores_tp_lyrata+1):(length_pores_tp_lyrata+zero_pores_tp_lyrata),c(2,3,5,6)]=-1
##pores_tp_lyrata[(length_pores_tp_lyrata+1):(length_pores_tp_lyrata+zero_pores_tp_lyrata),c(4,7)]=999
##pores_tp_lyrata[(length_pores_tp_lyrata+1):(length_pores_tp_lyrata+zero_pores_tp_lyrata),8]='a.lyrata'
##length_miseq_tp_thaliana=length(miseq_tp_thaliana[,1])
##miseq_tp_thaliana[(length_miseq_tp_thaliana+1):(length_miseq_tp_thaliana+zero_miseq_tp_thaliana),c(2,3,5,6)]=-1
##miseq_tp_thaliana[(length_miseq_tp_thaliana+1):(length_miseq_tp_thaliana+zero_miseq_tp_thaliana),c(4,7)]=999
##miseq_tp_thaliana[(length_miseq_tp_thaliana+1):(length_miseq_tp_thaliana+zero_miseq_tp_thaliana),8]='thaliana'
##length_miseq_tp_lyrata=length(miseq_tp_lyrata[,1])
##miseq_tp_lyrata[(length_miseq_tp_lyrata+1):(length_miseq_tp_lyrata+zero_miseq_tp_lyrata),c(2,3,5,6)]=-1
##miseq_tp_lyrata[(length_miseq_tp_lyrata+1):(length_miseq_tp_lyrata+zero_miseq_tp_lyrata),c(4,7)]=999
##miseq_tp_lyrata[(length_miseq_tp_lyrata+1):(length_miseq_tp_lyrata+zero_miseq_tp_lyrata),8]='a.lyrata'

# create species flags for TP
pores_tp_thaliana$TP='thaliana'
miseq_tp_thaliana$TP='thaliana'
  pores_tp_lyrata$TP='a.lyrata'
  miseq_tp_lyrata$TP='a.lyrata'

# create species flags for TP
pores_tp_thaliana$BLAST_1nident=(pores_tp_thaliana$BLAST_1pident/100) * pores_tp_thaliana$BLAST_1length
miseq_tp_thaliana$BLAST_1nident=(miseq_tp_thaliana$BLAST_1pident/100) * miseq_tp_thaliana$BLAST_1length
  pores_tp_lyrata$BLAST_1nident=(  pores_tp_lyrata$BLAST_1pident/100) *   pores_tp_lyrata$BLAST_1length
  miseq_tp_lyrata$BLAST_1nident=(  miseq_tp_lyrata$BLAST_1pident/100) *   miseq_tp_lyrata$BLAST_1length
pores_tp_thaliana$BLAST_2nident=(pores_tp_thaliana$BLAST_2pident/100) * pores_tp_thaliana$BLAST_2length
miseq_tp_thaliana$BLAST_2nident=(miseq_tp_thaliana$BLAST_2pident/100) * miseq_tp_thaliana$BLAST_2length
  pores_tp_lyrata$BLAST_2nident=(  pores_tp_lyrata$BLAST_2pident/100) *   pores_tp_lyrata$BLAST_2length
  miseq_tp_lyrata$BLAST_2nident=(  miseq_tp_lyrata$BLAST_2pident/100) *   miseq_tp_lyrata$BLAST_2length

# make inverts, e.g. direction of BLAST database comparison
# is reversed for A. lyrata samples
pores_tp_lyrata_invert = pores_tp_lyrata
pores_tp_lyrata_invert[,c(2:7,9,10)] = -1 * pores_tp_lyrata[,c(2:7,9,10)]
miseq_tp_lyrata_invert = miseq_tp_lyrata
miseq_tp_lyrata_invert[,c(2:7,9,10)] = -1 * miseq_tp_lyrata[,c(2:7,9,10)]

# bind species together
pores_all = rbind(pores_tp_thaliana,pores_tp_lyrata_invert)
miseq_all = rbind(miseq_tp_thaliana,miseq_tp_lyrata_invert)

# calculate the difference stats. values >0 indicate A.thaliana (all are T)
pores_all$diff_length=pores_all$BLAST_1length-pores_all$BLAST_2length
pores_all$diff_pident=pores_all$BLAST_1pident-pores_all$BLAST_2pident 
pores_all$diff_evalue=pores_all$BLAST_2evalue-pores_all$BLAST_1evalue 
pores_all$diff_nident=pores_all$BLAST_1nident-pores_all$BLAST_2nident 
miseq_all$diff_length=miseq_all$BLAST_1length-miseq_all$BLAST_2length
miseq_all$diff_pident=miseq_all$BLAST_1pident-miseq_all$BLAST_2pident
miseq_all$diff_evalue=miseq_all$BLAST_2evalue-miseq_all$BLAST_1evalue
miseq_all$diff_nident=miseq_all$BLAST_1nident-miseq_all$BLAST_2nident

# labels for classifier (in this case binary as above, e.g. positve==T negative==F)
#pores_all$pred_length_labels=pores_all$diff_length>0
#pores_all$pred_pident_labels=pores_all$diff_pident>0
#pores_all$pred_evalue_labels=pores_all$diff_evalue<0
#miseq_all$pred_length_labels=miseq_all$diff_length>0
#miseq_all$pred_pident_labels=miseq_all$diff_pident>0
#miseq_all$pred_evalue_labels=miseq_all$diff_evalue<0

# make a ROCR pred object
pores_all_pred_length = prediction(pores_all$diff_length,pores_all$TP)
pores_all_pred_pident = prediction(pores_all$diff_pident,pores_all$TP)
pores_all_pred_evalue = prediction(pores_all$diff_evalue,pores_all$TP)
pores_all_pred_nident = prediction(pores_all$diff_nident,pores_all$TP)
miseq_all_pred_length = prediction(miseq_all$diff_length,miseq_all$TP)
miseq_all_pred_pident = prediction(miseq_all$diff_pident,miseq_all$TP)
miseq_all_pred_evalue = prediction(miseq_all$diff_evalue,miseq_all$TP)
miseq_all_pred_nident = prediction(miseq_all$diff_nident,miseq_all$TP)



multiplot_ROC <- function(){
#multiplot
#par(mfrow=c(4,3))
par(mfrow=c(4,4))
## tp vs fp rates (classical ROC)
# e-value
plot(performance(pores_all_pred_evalue,measure="tpr",x.measure="fpr"),main="ROC: e-value",lwd=2)
plot(performance(miseq_all_pred_evalue,measure="tpr",x.measure="fpr"),add=T,col="red",lwd=2)
abline(0,1,lty=4)
# length
plot(performance(pores_all_pred_length,measure="tpr",x.measure="fpr"),main="ROC: length",lwd=2)
plot(performance(miseq_all_pred_length,measure="tpr",x.measure="fpr"),add=T,col="red",lwd=2)
abline(0,1,lty=4)
# pidents
plot(performance(pores_all_pred_pident,measure="tpr",x.measure="fpr"),main="ROC: % identities",lwd=2)
plot(performance(miseq_all_pred_pident,measure="tpr",x.measure="fpr"),add=T,col="red",lwd=2)
abline(0,1,lty=4)
# nidents
plot(performance(pores_all_pred_nident,measure="tpr",x.measure="fpr"),main="ROC: num. identities",lwd=2)
plot(performance(miseq_all_pred_nident,measure="tpr",x.measure="fpr"),add=T,col="red",lwd=2)
abline(0,1,lty=4)

## false-positive rate vs cutoff
# e-value
plot(performance(pores_all_pred_evalue,measure="fpr"),main="False-positive rate: e-value",xlim=c(-1e-20,1e-20),lwd=2)
plot(performance(miseq_all_pred_evalue,measure="fpr"),add=T,col="red",lwd=2)
# length
plot(performance(pores_all_pred_length,measure="fpr"),main="False-positive rate: length",xlim=c(-200,10),lwd=2)
plot(performance(miseq_all_pred_length,measure="fpr"),add=T,col="red",lwd=2)
# pidents
plot(performance(pores_all_pred_pident,measure="fpr"),main="False-positive rate: % identities",xlim=c(-30,10),lwd=2)
plot(performance(miseq_all_pred_pident,measure="fpr"),add=T,col="red",lwd=2)
# nidents
plot(performance(pores_all_pred_nident,measure="fpr"),main="False-positive rate: num. identities",xlim=c(-30,10),lwd=2)
plot(performance(miseq_all_pred_nident,measure="fpr"),add=T,col="red",lwd=2)

## true positive rate vs cutoff
# e-value
plot(performance(pores_all_pred_evalue,measure="tpr"),main="True-positive rate: e-value",xlim=c(-1e-20,1e-20),lwd=2)
plot(performance(miseq_all_pred_evalue,measure="tpr"),add=T,col="red",lwd=2)
# length
plot(performance(pores_all_pred_length,measure="tpr"),main="True-positive rate: length",xlim=c(-200,4000),lwd=2)
plot(performance(miseq_all_pred_length,measure="tpr"),add=T,col="red",lwd=2)
# pidents
plot(performance(pores_all_pred_pident,measure="tpr"),main="True-positive rate: % identities",xlim=c(-10,100),lwd=2)
plot(performance(miseq_all_pred_pident,measure="tpr"),add=T,col="red",lwd=2)
# nidents
plot(performance(pores_all_pred_nident,measure="tpr"),main="True-positive rate: num. identities",xlim=c(-10,100),lwd=2)
plot(performance(miseq_all_pred_nident,measure="tpr"),add=T,col="red",lwd=2)

## accuracy () vs cutoff
# e-value
plot(performance(pores_all_pred_evalue,measure="acc"),main="Accuracy: e-value",xlim=c(-1e-20,1e-20),lwd=2)
plot(performance(miseq_all_pred_evalue,measure="acc"),add=T,col="red",lwd=2)
# length
plot(performance(pores_all_pred_length,measure="acc"),main="Accuracy: length",xlim=c(-200,4000),lwd=2)
plot(performance(miseq_all_pred_length,measure="acc"),add=T,col="red",lwd=2)
# pidents
plot(performance(pores_all_pred_pident,measure="acc"),main="Accuracy: % identities",lwd=2)
plot(performance(miseq_all_pred_pident,measure="acc"),add=T,col="red",lwd=2)
# nidents
plot(performance(pores_all_pred_nident,measure="acc"),main="Accuracy: num. identities",lwd=2)
plot(performance(miseq_all_pred_nident,measure="acc"),add=T,col="red",lwd=2)
}

multiplot_ROC_2a <- function(){
# plotting for Figure Supplementary 2a, maximum info
par(oma=c(3,3,0,0),mar=c(5,5,1,1),mfrow=c(3,2),cex.axis=1.8,cex.lab=1.65,bty="n")
## false-positive rate vs cutoff
# length
plot(performance(pores_all_pred_length,measure="fpr"),xlab='',ylab='FP',main="",xlim=c(-200,10),lwd=4,bty="n")
plot(performance(pores_all_pred_length_NA,measure="fpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_length,measure="fpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_length_NA,measure="fpr"),add=T,col="pink",lty=3,lwd=4)
# pidents
plot(performance(pores_all_pred_pident,measure="fpr"),xlab='',ylab='',main="",xlim=c(-30,10),lwd=4,bty="n")
plot(performance(pores_all_pred_pident_NA,measure="fpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_pident,measure="fpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_pident_NA,measure="fpr"),add=T,col="pink",lty=3,lwd=4)

## true positive rate vs cutoff
# length
plot(performance(pores_all_pred_length,measure="tpr"),xlab='',ylab='TP',main="",xlim=c(-200,4000),lwd=4,bty="n")
plot(performance(pores_all_pred_length_NA,measure="tpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_length,measure="tpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_length_NA,measure="tpr"),add=T,col="pink",lty=3,lwd=4)
# pidents
plot(performance(pores_all_pred_pident,measure="tpr"),xlab='',ylab='',main="",xlim=c(-10,100),lwd=4,bty="n")
plot(performance(pores_all_pred_pident_NA,measure="tpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_pident,measure="tpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_pident_NA,measure="tpr"),add=T,col="pink",lty=3,lwd=4)

## accuracy () vs cutoff
# length
plot(performance(pores_all_pred_length,measure="acc"),xlab="Length",ylab='Accuracy',main="",xlim=c(-200,4000),lwd=4,bty="n")
plot(performance(pores_all_pred_length_NA,measure="acc"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_length,measure="acc"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_length_NA,measure="acc"),add=T,col="pink",lty=3,lwd=4)
# pidents
plot(performance(pores_all_pred_pident,measure="acc"),xlab="% Identities",ylab='',main="",lwd=4,bty="n")
plot(performance(pores_all_pred_pident_NA,measure="acc"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_pident,measure="acc"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_pident_NA,measure="acc"),add=T,col="pink",lty=3,lwd=4)

# axes
mtext(text="Cut-off",side=1,line=0,outer=TRUE,cex=1.9)
mtext(text="Rate",side=2,line=0,outer=TRUE,cex=1.9)
}


# plotting for Figure 2b, original style 
## accuracy () vs cutoff
# length
par(mfrow=c(1,1),cex.axis=1.8,cex.lab=1.65,bty="n",mar=c(5,5,1,1))
plot(performance(pores_all_pred_length,measure="acc"),xlab="Length bias",ylab='Accuracy',main="",xlim=c(0,2000),lwd=6,bty="n")
plot(performance(pores_all_pred_length_NA,measure="acc"),add=T,col="grey",lty=3,lwd=6)
plot(performance(miseq_all_pred_length,measure="acc"),add=T,col="red",lwd=6)
plot(performance(miseq_all_pred_length_NA,measure="acc"),add=T,col="pink",lty=3,lwd=6)

plot_2a2b_idents <- function(){
#/* repeated for (a,b), (c,d) format */
#/* plotting for Figure 2b, original style */
# 2a, 2b: %idents
par(mfrow=c(2,1),cex.axis=1.8,cex.lab=1.65,bty="n",mar=c(5,5,1,1))
plot(performance(pores_all_pred_pident,measure="tpr"),xlab='',ylab='TP',main="",xlim=c(-10,100),lwd=4,bty="n")
plot(performance(pores_all_pred_pident_NA,measure="tpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_pident,measure="tpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_pident_NA,measure="tpr"),add=T,col="pink",lty=3,lwd=4)
plot(performance(pores_all_pred_pident,measure="acc"),xlab="% Identities bias",ylab='Accuracy',main="",lwd=4,bty="n",xlim=c(-25,100))
plot(performance(pores_all_pred_pident_NA,measure="acc"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_pident,measure="acc"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_pident_NA,measure="acc"),add=T,col="pink",lty=3,lwd=4)
}


plot_2a2b_length <- function(){
# 2c, 2d: length
par(mfrow=c(2,1),cex.axis=1.8,cex.lab=1.65,bty="n",mar=c(5,5,1,1))
plot(performance(pores_all_pred_length,measure="tpr"),xlab='',ylab='TP',main="",xlim=c(0,2000),lwd=4,bty="n")
plot(performance(pores_all_pred_length_NA,measure="tpr"),add=T,col="grey",lty=3,lwd=4)
plot(performance(miseq_all_pred_length,measure="tpr"),add=T,col="red",lwd=4)
plot(performance(miseq_all_pred_length_NA,measure="tpr"),add=T,col="pink",lty=3,lwd=4)
plot(performance(pores_all_pred_length,measure="acc"),xlab="Length bias",ylab='Accuracy',main="",xlim=c(0,2000),lwd=6,bty="n")
plot(performance(pores_all_pred_length_NA,measure="acc"),add=T,col="grey",lty=3,lwd=6)
plot(performance(miseq_all_pred_length,measure="acc"),add=T,col="red",lwd=6)
plot(performance(miseq_all_pred_length_NA,measure="acc"),add=T,col="pink",lty=3,lwd=6)
}

# finally, run the plotting calls
multiplot_ROC()
break()
multiplot_ROC_2a()
plot_2a2b_idents()
plot_2a2b_length()
