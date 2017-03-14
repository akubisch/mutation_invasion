require(caTools)

occ_d02   <- as.matrix(read.table("disp_02/processed_results/occ.txt"))
occ_mutfix  <- as.matrix(read.table("mut_4/processed_results/occ.txt"))
occ_mutlowfix <- as.matrix(read.table("mut_5/processed_results/occ.txt"))
occ_d02_mutfix <- as.matrix(read.table("disp_02_fixed_mut/processed_results/occ.txt"))
occ_d02_lowmutfix <- as.matrix(read.table("disp_02_fixed_lowmut/processed_results/occ.txt"))
occ  <- as.matrix(read.table("dispevol/processed_results/occ.txt"))


R_mutfix <- numeric()
for (t in 1:150) {
	o <- rev(occ_mutfix[t,])
	R_mutfix[t] <- length(o[!is.na(o)])
}

R <- numeric()
for (t in 1:150) {
	o <- rev(occ[t,])
	R[t] <- length(o[!is.na(o)])
}

R_lowmutfix <- numeric()
for (t in 1:150) {
	o <- rev(occ_mutlowfix[t,])
	R_lowmutfix[t] <- length(o[!is.na(o)])
}

R_dmutfix <- numeric()
for (t in 1:150) {
  o <- rev(occ_d02_mutfix[t,])
  R_dmutfix[t] <- length(o[!is.na(o)])
}

R_dlowmutfix <- numeric()
for (t in 1:150) {
  o <- rev(occ_d02_lowmutfix[t,])
  R_dlowmutfix[t] <- length(o[!is.na(o)])
}

R_dfix <- numeric()
for (t in 1:150) {
  o <- rev(occ_d02[t,])
  R_dfix[t] <- length(o[!is.na(o)])
}

x11(width=3.75,height=6.67)

wd <- 20

ts <- c(5,10,15,50)
cols <- gray(seq(0.8,0,len=length(ts)))
cntbreak <- 20

par(mfrow=c(2,1),bty="l",mar=c(2,6,2,1.5),lwd=2,cex.axis=1.25,cex.lab=1.5,oma=c(4,0,0,0))


plot(1,1,type="n",xlab="time (generations)",ylab="range border position",xlim=c(0,5000),ylim=c(75,250),main="         evolving dispersal")

lines(runmean(R[101:150],00)~seq(1,5000,len=50),lwd=2)
lines(runmean(R_mutfix[101:150],00)~seq(1,5000,len=50),lwd=2,lty="dashed")
lines(runmean(R_lowmutfix[101:150],00)~seq(1,5000,len=50),lwd=2,lty="dotted")

#legend("bottomright",bty="n",lwd=2,lty=c("solid","dashed","dotted"),cex=1.25,c("control",expression(m==10^-4),expression(m==10^-5)))
text(200,par("usr")[[4]]*1.05,"A",xpd=T,cex=1.75)


plot(1,1,type="n",xlab="time (generations)",ylab="range border position",xlim=c(0,5000),ylim=c(75,250),main="         fixed dispersal")

lines(runmean(R_dfix[101:150],00)~seq(1,5000,len=50),lwd=2)
lines(runmean(R_dmutfix[101:150],00)~seq(1,5000,len=50),lwd=2,lty="dashed")
lines(runmean(R_dlowmutfix[101:150],00)~seq(1,5000,len=50),lwd=2,lty="dotted")

legend("bottomright",bty="n",lwd=2,lty=c("solid","dashed","dotted"),cex=1.15,c("m = evol.",expression(m==10^-4),expression(m==10^-5)))
text(200,par("usr")[[4]]*1.05,"B",xpd=T,cex=1.75)

title(xlab="time (generations)", outer=T, line=1.5, adj=0.6)

dev.copy2eps(file="figure_2.eps",title="Cobben & Kubisch - Figure 2")
