require(caTools)

dens <- as.matrix(read.table("processed_results/dens.txt"))
disp <- as.matrix(read.table("processed_results/disp.txt"))
evol <- as.matrix(read.table("processed_results/evol.txt"))
div  <- as.matrix(read.table("processed_results/div.txt"))
neut_div  <- as.matrix(read.table("processed_results/neut_div.txt"))
disp_div <- as.matrix(read.table("processed_results/disp_div.txt"))
cnts <- as.matrix(read.table("processed_results/cnts.txt"))

# figure 5

x11(width=7.5,height=10)

wd <- 20

ts <- c(1,3,5,10,15,50)
cols <- gray(seq(0.8,0,len=length(ts)))
cntbreak <- 10

par(mfrow=c(3,2),bty="l",mar=c(1,6,2,1),oma=c(5,0,0,0),lwd=2,cex.axis=1.75,cex.lab=1.85)

plot(1,1,type="n",xlab="",ylab="population density",xlim=c(0,250),ylim=c(.2,.7),xaxt="n")
axis(side=1,labels=F)

for (t in 1:length(ts)) {
	m <- as.numeric(dens[ts[t],][which(cnts[ts[t],]>=cntbreak)])

	lines(runmean(m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.03,"A",xpd=T,cex=2)

plot(1,1,type="n",xlab="",ylab="dispersal probability",xlim=c(0,250),ylim=c(.1,.2),xaxt="n",yaxt="n")
axis(side=2,at=c(0.1,0.15,0.2))
axis(side=1,labels=F)

for (t in 1:length(ts)) {
	m <- as.numeric(disp[ts[t],][which(cnts[ts[t],]>=cntbreak)])

	lines(runmean(m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.025,"B",xpd=T,cex=2)


plot(1,1,type="n",xlab="",ylab=expression(paste("mutation rate of ",tau[opt]," (*",10^-4,")",sep="")),
     xlim=c(0,250),ylim=c(0.00,15),xaxt="n",yaxt="n")
axis(side=1,labels=F)
axis(side=2,at=c(0,5,10,15,20))

for (t in 1:length(ts)) {
	m <- as.numeric(evol[ts[t],][which(cnts[ts[t],]>=cntbreak)])
  x <- 10000*10^-m
  if (length(which(x==10000))>0) {x <- x[-which(x==10000)]}

	lines(runmean(x,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.05,"C",xpd=T,cex=2)

plot(1,1,type="n",xlab="",ylab=expression(paste("genetic diversity (*",10^-2,")",sep="")),
     xlim=c(0,250),ylim=range(100*div,na.rm=T),xaxt="n")
axis(side=1,labels=F)

for (t in 1:length(ts)) {
	m <- as.numeric(div[ts[t],][which(cnts[ts[t],]>=cntbreak)])

	lines(runmean(100*m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.05,"D",xpd=T,cex=2)

plot(1,1,type="n",xlab="",ylab=expression(paste("disp. locus gen. diversity (*",10^-2,")",sep="")),
     xlim=c(0,250),ylim=c(0.0,0.2))

for (t in 1:length(ts)) {
  m <- as.numeric(disp_div[ts[t],][which(cnts[ts[t],]>=cntbreak)])

  lines(runmean(100*m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.05,"E",xpd=T,cex=2)

plot(1,1,type="n",xlab="",ylab=expression(paste("neutral genetic diversity (*",10^-2,")",sep="")),
     xlim=c(0,250),ylim=c(0,40))

for (t in 1:length(ts)) {
	m <- as.numeric(neut_div[ts[t],][which(cnts[ts[t],]>=cntbreak)])

	lines(runmean(100*m,wd),col=cols[t])
}
text(10,par("usr")[[4]]*1.05,"F",xpd=T,cex=2)
legend("topright",bty="n",lwd=2,col=cols,c("t = 100", "t = 300", "t = 500", "t = 1000", "t = 1500", "t = 5000"),cex=1.75)

#plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")

title(xlab="spatial location",outer=T)

dev.copy2eps(file="figure5.eps",title="Cobben & Kubisch - Figure 5")
