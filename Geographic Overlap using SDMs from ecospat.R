for (sp in seq_along(species_list))
{
  ### Subset the DATABASE per species
  species_records <- subset(species_DB, Species==species_list[sp])
  ### Sort the time intervals for the selected species
  Species_time_bins <- sort(unique(species_records$Interval))
  spec_nam <- as.vector(unique(species_records$Species))
  vec <- as.vector(NULL)
  for (bi in seq_along(Species_time_bins))
  { 
    Occurrences_per_bin <- subset(species_records, Interval==Species_time_bins[bi])	 ## a subset data.frame for each time bin ###
    vec <- c(vec, nrow(Occurrences_per_bin)) ### Construct a vector with the number of Occurrences for each time interval
  }
  list_len <- length(vec[vec >= 2*n_var])
  sdm_bins <- Species_time_bins[vec >= 2*n_var] ### Select the time intervals with enough occurrences to run the analysis
  used_occ <- vec[vec >= 2*n_var]
  if(list_len > 0)
  {
  #### have to be specified based on the unique bins of the species
  for(bi in seq_along(sdm_bins)) 
    {
      clim_per_bin <- subset(clim.all, Interval==sdm_bins[bi])
      temp_hum_per_bin <- subset(Human_DB, Interval==sdm_bins[bi])
      #print(nrow(temp_hum_per_bin))
      temp_sp_per_bin <- subset(species_records, Interval==sdm_bins[bi]) 
      #print(nrow(temp_sp_per_bin))
      temp_row_species <- sample.sp.globvar(dfsp=clim_per_bin,colspxy=1:2,colspkept=NULL,dfvar=temp_sp_per_bin,colvarxy=10:9,colvar=13,resolution=0.1) ### use the grid cell coordinates (specificaly: x=cell.Latitude, y=cell.Longitude) 
      temp_row_hums <- sample.sp.globvar(dfsp=clim_per_bin,colspxy=1:2,colspkept=NULL,dfvar=temp_hum_per_bin,colvarxy=10:9,colvar=13,resolution=0.1)
      assign(paste("row.pa.sp.",sdm_bins[bi], sep=""), temp_row_species) ### rows of presence absence/background values for the species
      assign(paste("row.pa.hums.", sdm_bins[bi], sep=""), temp_row_hums)  ### rows of presence absence/background values for humans
    }
	sdm_bins2 <- as.vector(NULL)
	for(ro in seq_along(ls(pattern = "row.pa.sp.")))
    {
      clim_per_bin2 <- subset(clim.all, Interval==sdm_bins[ro])
      temp_sp_pa <- data.frame((!is.na(get(ls(pattern="row.pa.sp")[ro])))*1);names(temp_sp_pa)<-c("pa1") ## convert to presence absence data frame
      temp_sp_pa_b <- cbind(clim_per_bin2, temp_sp_pa)
      print(nrow(temp_sp_pa_b[temp_sp_pa_b$pa1>0,]))  ### the number of "unique" occurrecnes (the numnber of grid cell where we have presences)
      if(nrow(temp_sp_pa_b[temp_sp_pa_b$pa1>0,]) >= 6)
      {
        assign(paste("species.pa.", sdm_bins[ro], "b", sep=""), temp_sp_pa_b)
        temp_hum_pa <- data.frame((!is.na(get(ls(pattern="row.pa.hums.")[ro])))*1);names(temp_hum_pa)<-c("pa1") ## convert to presence absence data frame
        temp_hum_pa_b <- cbind(clim_per_bin2, temp_hum_pa)
        print(nrow(temp_hum_pa_b[temp_hum_pa_b$pa1>0,])) ### the number of "unique" occurrecnes (the numnber of grid cell where we have presences)
        assign(paste("human.pa.", sdm_bins[ro], "b", sep=""), temp_hum_pa_b)
        sdm_bins2 <- c(sdm_bins2, sdm_bins[ro])
      }
    }
    sdm_bins <- sdm_bins2 ## We keep only those time bins for which we have minimum 6 unique occurrences
    PROJ = F
    names(clim.all)
    Xvar<-match(vars, names(clim.all)) ### We have to specify which of all the variables we are going to use again
    nvar<-length(Xvar)
    iterations<-500
    R=100
	if(PROJ == F)
    {
      for (ro in seq_along(ls(pattern = "human.pa")))
      {
        human_data <- get(ls(pattern = "human.pa")[ro])
        species_data <- get(ls(pattern = "species.pa")[ro])
        par(mfrow=c(1,2), mar=c(2,2,2,2))
		## humans Maxent ###
        human_xyz <- human_data[,c(2,1,Xvar[1])];names(human_xyz)<-names(human_data[,c(1,2,Xvar[1])])
        human_rasts <- rasterFromXYZ(human_xyz)
        for(r in Xvar[2:nvar]) ## will depend based on the number of variables
        {
          human_xyz <- human_data[,c(2,1,r)]
          human_nrast <- rasterFromXYZ(human_xyz)
          human_rasts <- stack(human_rasts, human_nrast)
        }
        human_occur <- human_data[which(human_data$pa1==1),2:1]
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='') #check if maxent.jar is located in the right folder
        human_me <- maxent(human_rasts, human_occur)
        human_r<-predict(human_me, human_rasts)
        
        ## species Maxent ##
        species_xyz <- species_data[,c(2,1,Xvar[1])];names(species_xyz)<-names(species_data[,c(1,2,Xvar[1])])
        species_rasts <- rasterFromXYZ(species_xyz)
        for(r in Xvar[2:nvar])
		{
          species_xyz <- species_data[,c(2,1,r)]
          species_nrast <- rasterFromXYZ(species_xyz)
          species_rasts <- stack(species_rasts, species_nrast)
        }
        species_occur <- species_data[which(species_data$pa1==1),2:1]
        jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        species_me <- maxent(species_rasts, species_occur)
        species_r<-predict(species_me, species_rasts)
        temp_clim <- subset(clim.all, Interval==sdm_bins[ro])
        temp_clim2 <- temp_clim[,c(2,1,3:10)]; names(temp_clim2)<-names(temp_clim[,c(1:10)])
        temp_species <- species_data[species_data$pa1>0,]
        temp_species2 <- temp_species[,c(2,1,3:10)]; names(temp_species2)<-names(temp_species[,c(1,2,3:10)])
        temp_humans <- human_data[human_data$pa1>0,]
        temp_humans2 <- temp_humans[,c(2,1,3:10)]; names(temp_humans2)<-names(temp_humans[,c(1,2,3:10)])
        assign(paste("clim", sdm_bins[ro], "b", sep=""), temp_clim2)
        assign(paste("hum.occ.", sdm_bins[ro], "kb", sep=""), temp_humans2)
        assign(paste("spec.occ.", sdm_bins[ro], "kb", sep=""), temp_species2)
        ## HUMANS MODELING ###
        clim.allb <- temp_clim2
        #clim.allb<-rbind(temp_clim2,temp_clim2) ### here we are using both ranges to calibrate the model
        scores.clim.all.MAXENT.humans <- data.frame(predict(human_me, clim.allb[,Xvar], progress="text"))
        temp_scores_clim_hum <- data.frame(predict(human_me, temp_clim2[,Xvar], progress="text"))
        assign(paste("hum.scores.clim.", sdm_bins[ro], "k.MAXENT", sep=""), temp_scores_clim_hum)
        temp_hum_occ <- data.frame(predict(human_me, get(paste("hum.occ.", sdm_bins[ro],"kb", sep=""))[,Xvar], progress="text"))
        assign(paste("scores.humans.", sdm_bins[ro], "k.MAXENT", sep=""), temp_hum_occ)
        zz_hum <- grid.clim(scores.clim.all.MAXENT.humans, get(ls(pattern="hum.scores.clim.")[ro]), get(ls(pattern="scores.humans.")[ro]),R)
        assign(paste("humans_z", sdm_bins[ro], sep=""), zz_hum)
        ###  MAMMALS MODELLING ###
        scores.clim.all.MAXENT.species <- data.frame(predict(species_me, clim.allb[,Xvar], progress="text"))
        temp_scores_clim_spec <- data.frame(predict(species_me, temp_clim2[,Xvar], progress="text"))
        assign(paste("spec.scores.clim.", sdm_bins[ro], "k.MAXENT", sep=""), temp_scores_clim_spec)
        temp_species_occ <- data.frame(predict(species_me, get(paste("spec.occ.", sdm_bins[ro],"kb", sep=""))[,Xvar], progress="text"))
        assign(paste("scores.species.", sdm_bins[ro], "k.MAXENT", sep=""), temp_species_occ)
        zz_spec <- grid.clim(scores.clim.all.MAXENT.species, get(ls(pattern="spec.scores.clim.")[ro]),get(ls(pattern="scores.species.")[ro]),R)
        assign(paste("species_z", sdm_bins[ro], sep=""), zz_spec)
      }
    }
	setwd("C:/Users/Dest/Desktop/Project 3 PHD/ANALYSES/Kostis_data/All sequences/Fossil data for SDMs/Geographic Overlap New humans data")
	Overlap <- as.vector(NULL)
    len <- length(sdm_bins)
    for (z in seq_along(sdm_bins))
    {
      a <-niche.equivalency.test(get(rev(ls(pattern="humans_z"))[z]), get(rev(ls(pattern="species_z"))[z]), rep=500)
      assign(paste("a_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep=""), a) ### speciffy soemwhere for which species is this
      b <-niche.similarity.test(get(rev(ls(pattern="humans_z"))[z]), get(rev(ls(pattern="species_z"))[z]), rep=500)
      assign(paste("b_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep=""), b) ### speciffy soemwhere for which species is this
      b2 <-niche.similarity.test(get(rev(ls(pattern="species_z"))[z]), get(rev(ls(pattern="humans_z"))[z]), rep=500)
      assign(paste("b2_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep=""), b2)
      layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
      plot.niche(get(rev(ls(pattern="humans_z"))[z]), title= paste0("Maxent humans- ", rev(sdm_bins)[z],"k_niche", sep=""), name.axis1="Probability of occurence",name.axis2="PC2")
      box()
      plot.niche(get(rev(ls(pattern="species_z"))[z]), title= paste0("Maxent - ",species_list[sp], "_", rev(sdm_bins)[z], "k_niche", sep=""), name.axis1="Probability of occurence",name.axis2="PC2")
      box()
      plot(species_me)
      plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(niche.overlap(get(rev(ls(pattern="humans_z"))[z]), get(rev(ls(pattern="species_z"))[z]),cor=T)[1]),3)))
      Over <- round(as.numeric(niche.overlap(get(rev(ls(pattern="humans_z"))[z]), get(rev(ls(pattern="species_z"))[z]),cor=T)[1]),3)
      Overlap <- c(Overlap, Over)
      plot.overlap.test(get(paste("a_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep="")),"D","Equivalency")
      plot.overlap.test(get(paste("b_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep="")),"D","Similarity 2->1")
      plot.overlap.test(get(paste("b2_", species_list[sp], "_", rev(sdm_bins)[z], "k", sep="")),"D","Similarity 1->2")
      #dev.off()
    }
	Overlap2 <- cbind.data.frame(Time=rev(sdm_bins), Overlap) ### check the 
    assign(paste("Geographic_overlap_", species_list[sp], "_humans", sep=""), Overlap2)
    write.table(Overlap2, file=paste("Geographic_overlap_", species_list[sp], "_humans.txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    rm(list=c(ls(pattern = "row.pa."), ls(pattern="human.pa."), ls(pattern = "species.pa."), ls(pattern = "hum.occ."), ls(pattern = "spec.occ."),
              ls(pattern = "hum.scores.clim."), ls(pattern = "scores.humans."), ls(pattern = "spec.scores.clim."), ls(pattern = "scores.species."), 
              ls(pattern = "humans_z"), ls(pattern = "species_z"), ls(pattern = "clim_")))
	}
	setwd("C:/Users/kgiab/Desktop/Project/Project 3 PHD/ANALYSES/Kostis_data/All sequences/Fossil data for SDMs/Geographic Overlap New humans data")
}
