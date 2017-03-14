require(caTools)

evol_d005 <- as.matrix(read.table("disp_005/processed_results/evol.txt"))
cnts_d005 <- as.matrix(read.table("disp_005/processed_results/cnts.txt"))
evol_d01  <- as.matrix(read.table("disp_01/processed_results/evol.txt"))
cnts_d01 <- as.matrix(read.table("disp_01/processed_results/cnts.txt"))
evol_d02  <- as.matrix(read.table("disp_02/processed_results/evol.txt"))
cnts_d02 <- as.matrix(read.table("disp_02/processed_results/cnts.txt"))

x11(width=3.75,height=10)

wd <- 10

ts <- c(5,10,15,50)
cols <- gray(seq(0.8,0,len=length(ts)))
cntbreak <- 20

par(mfrow=c(3,1),bty="l",mar=c(1.5,6,2,1.5),lwd=2,cex.axis=1.5,cex.lab=1.75,oma=c(4,0,0,0))

plot(1,1,type="n",xlab="",ylab=expression(paste("mutation rate of ",tau[opt]," (*",10^-4,")",sep="")),xlim=c(0,250),ylim=c(.5,20),xaxt="n",
	 main=expression(d==0.05),cex.main=1.75)
#axis(side=1,labels=F)

for (t in 1:length(ts)) {

	m <- as.numeric(evol_d005[ts[t],][which(cnts_d005[ts[t],]>=cntbreak)])
	if (length(which(m==0))>0) {m <- m[-which(m==0)]}

	lines(runmean(10000*10^-m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.03,"A",xpd=T,cex=2)

legend("topright",col=cols,bty="n",cex=1.35,lwd=2,c("t=500","t=1000","t=1500","t=5000"))

plot(1,1,type="n",xlab="",ylab=expression(paste("mutation rate of ",tau[opt]," (*",10^-4,")",sep="")),xlim=c(0,250),ylim=c(.5,20),xaxt="n",
	 main=expression(d==0.1),cex.main=1.75)
#axis(side=1,labels=F)

for (t in 1:length(ts)) {

	m <- as.numeric(evol_d01[ts[t],][which(cnts_d01[ts[t],]>=cntbreak)])
	if (length(which(m==0))>0) {m <- m[-which(m==0)]}

	lines(runmean(10000*10^-m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.03,"B",xpd=T,cex=2)


plot(1,1,type="n",xlab="",ylab=expression(paste("mutation rate of ",tau[opt]," (*",10^-4,")",sep="")),xlim=c(0,250),ylim=c(.5,20),
	 main=expression(d==0.2),cex.main=1.75)
#axis(side=1,labels=T)

for (t in 1:length(ts)) {

	m <- as.numeric(evol_d02[ts[t],][which(cnts_d02[ts[t],]>=cntbreak)])
	if (length(which(m==0))>0) {m <- m[-which(m==0)]}

	lines(runmean(10000*10^-m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.03,"C",xpd=T,cex=2)

title(xlab="spatial location", outer=T, line=2, adj=0.6)

dev.copy2eps(file="figure_2.eps", title="Cobben & Kubisch - Figure 2")
