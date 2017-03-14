repli <- 200
ts <- 150

dens <- matrix(0,nrow=ts-100,ncol=250)
cnts <- matrix(0,nrow=ts-100,ncol=250)
disp <- matrix(0,nrow=ts-100,ncol=250)
evol <- matrix(0,nrow=ts-100,ncol=250)
div  <- matrix(0,nrow=ts-100,ncol=250)
neut_div  <- matrix(0,nrow=ts-100,ncol=250)
disp_div  <- matrix(0,nrow=ts-100,ncol=250)

geom_mat_mean <- function(mat) {
  mat[mat==0] <- NA
  xdim <- dim(mat)[[2]]
  res <- numeric(xdim)
  for (x in 1:xdim) {
    y <- mat[,x]; y <- y[which(!is.na(y))]
    m <- prod(y,na.rm=T)^(1/length(y))
    res[x] <- m
  }
  return(res)
}


for (t in 101:ts) {
  print(t)
  denst <- matrix(NA,nrow=repli,ncol=250)
  dispt <- matrix(NA,nrow=repli,ncol=250)
  evolt <- matrix(NA,nrow=repli,ncol=250)
  divt <- matrix(NA,nrow=repli,ncol=250)
  neut_divt <- matrix(NA,nrow=repli,ncol=250)
  disp_divt <- matrix(NA,nrow=repli,ncol=250)
  cntst <- matrix(NA,nrow=repli,ncol=250)

  for (r in 1:repli) {
    print(r)
    densr <- as.matrix(read.table(paste("results/densities_",r,".txt",sep="")))
    occr <- as.matrix(read.table(paste("results/occupancy_",r,".txt",sep="")))
    dispr <- as.matrix(read.table(paste("results/dispersal_",r,".txt",sep="")))
    evolr <- as.matrix(read.table(paste("results/evolvability_",r,".txt",sep="")))
    divr  <- as.matrix(read.table(paste("results/diversity_",r,".txt",sep="")))
    neut_divr  <- as.matrix(read.table(paste("results/neutral_diversity_",r,".txt",sep="")))
    disp_divr <- as.matrix(read.table(paste("results/disp_diversity_",r,".txt",sep="")))
    
    denst[r,] <- densr[t,]
    dispt[r,] <- dispr[t,]
    evolt[r,] <- evolr[t,]
    divt[r,] <- divr[t,]
    neut_divt[r,] <- neut_divr[t,]
    disp_divt[r,] <- disp_divr[t,]
    cntst[r,] <- occr[t,]>0
  }
  
  dens[t-100,] <- apply(denst,2,mean,na.rm=T)
  disp[t-100,] <- apply(dispt,2,mean,na.rm=T)
  evol[t-100,] <- geom_mat_mean(evolt)
  div[t-100,] <- apply(divt,2,mean,na.rm=T)
  neut_div[t-100,] <- apply(neut_divt,2,mean,na.rm=T)
  disp_div[t-100,] <- apply(disp_divt,2,mean,na.rm=T)
  cnts[t-100,] <- apply(cntst,2,sum)
}

write.table(dens,file="processed_results/dens.txt",col.names=FALSE,row.names=FALSE)
write.table(disp,file="processed_results/disp.txt",col.names=FALSE,row.names=FALSE)
write.table(evol,file="processed_results/evol.txt",col.names=FALSE,row.names=FALSE)
write.table(div,file="processed_results/div.txt",col.names=FALSE,row.names=FALSE)
write.table(neut_div,file="processed_results/neut_div.txt",col.names=FALSE,row.names=FALSE)
write.table(cnts,file="processed_results/cnts.txt",col.names=FALSE,row.names=FALSE)
write.table(disp_div,file="processed_results/disp_div.txt",col.names=F,row.names=F)

























