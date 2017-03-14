write("simulation package initiated...", file="progress.sim", append=F)

for (r in 1:200) {
	write(paste(" -> replication:",r,sep=" "),file="progress.sim",append=T)

	system("./range_exp")

	res_file  <- "data/adaptation.txt"
  	new_loc   <- paste("results/adaptation_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

  	res_file  <- "data/densities.txt"
  	new_loc   <- paste("results/densities_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

  	res_file  <- "data/dispersal.txt"
  	new_loc   <- paste("results/dispersal_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

  	res_file  <- "data/diversity.txt"
  	new_loc   <- paste("results/diversity_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

  	res_file  <- "data/neutral_diversity.txt"
  	new_loc   <- paste("results/neutral_diversity_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)
  	
  	res_file  <- "data/mut_diversity.txt"
  	new_loc   <- paste("results/mut_diversity_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)
  	
  	res_file  <- "data/disp_diversity.txt"
  	new_loc   <- paste("results/disp_diversity_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)
  	
  	res_file  <- "data/evolvability.txt"
  	new_loc   <- paste("results/evolvability_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

  	res_file  <- "data/occupancy.txt"
  	new_loc   <- paste("results/occupancy_",r,".txt",sep="")
  	file.copy(res_file,new_loc,overwrite=T)

}

write("done.",file="progress.sim",append=T)
