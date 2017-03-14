dens <- as.matrix(read.table("processed_results/dens.txt"))
disp <- as.matrix(read.table("processed_results/disp.txt"))
evol <- as.matrix(read.table("processed_results/evol.txt"))

tmax <- dim(dens)[[1]]

dens_t <- numeric()
disp_t <- numeric()
evol_t <- numeric()

for (t in 1:tmax) {
  dens_t[t] <- mean(dens[t,],na.rm=T)
  disp_t[t] <- mean(disp[t,],na.rm=T)
  evol_t[t] <- 10000*10^-mean(evol[t,],na.rm=T)
}

ts <- 1:50*100

x11(width=10,height=3.75)
par(mfrow=c(1,3),bty="l",mar=c(1,6,3,1),oma=c(5,0,0,0),lwd=2,cex.axis=1.75,cex.lab=1.85)

plot(1,1,type="n",xlim=range(ts),ylim=c(0,0.6),xlab="",ylab="population density",yaxt="n",xaxt="n")
axis(side=1,at=c(0,2500,5000))
axis(side=2,at=c(0,0.3,0.6))
lines(dens_t~ts)
text(250,par("usr")[[4]]*1.005,"A",xpd=T,cex=2)

plot(1,1,type="n",xlim=range(ts),ylim=c(0.15,0.25),xlab="",ylab="dispersal probability",yaxt="n",xaxt="n")
axis(side=1,at=c(0,2500,5000))
axis(side=2,at=c(0.15,0.2,0.25))
lines(disp_t~ts)
text(250,par("usr")[[4]]*1.01,"B",xpd=T,cex=2)

plot(1,1,type="n",xlim=range(ts),ylim=c(0.0,15),xlab="",ylab=expression(paste("mutation rate of ",tau[opt]," (*",10^-4,")",sep="")),yaxt="n",xaxt="n")
axis(side=1,at=c(0,2500,5000))
axis(side=2,at=c(0,5,10,15))
lines(evol_t~ts)
text(250,par("usr")[[4]]*1.015,"C",xpd=T,cex=2)


title(xlab="time (generations)",outer=T)

dev.copy2eps(file="figure4a.eps",title="Cobben & Kubisch - Figure 4a")
