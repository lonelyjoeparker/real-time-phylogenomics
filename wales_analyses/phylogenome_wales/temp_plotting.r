#read in data
rtp=read.table("Downloads/serial_temperature_RTP.log",header=T,sep="\t")
cold=read.table("Downloads/serial_temperature_smallfridgebox.log",header=T,sep="\t")

#approximately convert diode mV to temp celsius
colds=(107-cold$bias_mV)*2.1
warms=(107-rtp$bias_mV)*2.1

#create time-series
cold_ser=ts(colds)
cold_hw=HoltWinters(cold_ser,beta=F,gamma=F)
rtp_ser=ts(warms)
rtp_hw=HoltWinters(rtp_ser,beta=F,gamma=F)

#open new device and enable 2-column plotting
quartz()
par(mfrow=c(1,2))
par(cex=0.8)

#plot both series, managing y-axis (temp) limits manually
plot(rtp_hw,main="Room temperature, Holt-winters trend",xlab=("Time (22hrs total)"),ylim=c(0,16),ylab=c("temp celsius (approx)"))
plot(cold_hw,main="Icebox warming up, Holt-winters trend",xlab=("Time (48hrs total)"),ylim=c(0,16),ylab=c("temp celsius (approx)"))
