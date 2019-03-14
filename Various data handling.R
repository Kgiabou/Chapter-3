require(biomod2)
require(ggplot2)
require(gridExtra)
require(SDMTools)
require(raster)
require(rasterVis)
require(ade4)

bins <-c(seq(8,22,1), seq(24,46,2))
par(mfrow=c(3,1), mar=c(2,2,2,2), oma=c(2,2,2,2))
for(i in seq_along(bins))
{
  Climate <- read.delim(dir(pattern = paste("clim_", bins[i], "kyr", sep = "")),h=T,sep="\t", stringsAsFactors = FALSE, quote="")
  colnames(Climate) <- c("Latitude", "Longitude", "Prec_Sp", "Prec_Su", "Prec_Au", "Prec_Wi", "Temp_Sp", "Temp_Su", "Temp_Au", "Temp_Wi")
  Interval<- bins[i]
  ### Convert climatic dataframe to raster #####
  Climate3 <- na.exclude(Climate)
  Climatexyz <- Climate3[,c(2,1,3)]; names(Climatexyz) <- names(Climatexyz[,c(1,2,3)])
  clim_rasts <- rasterFromXYZ(Climatexyz)
  plot(clim_rasts, main=paste0(bins[i],"kyr_",names(clim_rasts)[[1]], sep=""), cex.main = 0.9, xlim=c(-180,180), ylim=c(30,90))
  for(r in c(4:10)) ## also will depend based on the variables
  {
    clim_xyz <- Climate3[,c(2,1,r)]
    clim_nrast <- rasterFromXYZ(clim_xyz)
    clim_rasts <- stack(clim_rasts, clim_nrast)
    crs(clim_rasts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  }
  assign(paste("clim_", Interval,"k", sep=""), clim_rasts)
}

### PCA to assess the collinearity and importance of climatic variables ###
setwd("C:/Users/kgiab/Desktop/Project 3 PHD/ANALYSES/Kostis_data/All sequences/Fossil data for SDMs/Geographic Overlap New humans data")


library(ade4)
library(factoextra)

humans_final3 <- read.delim("Humans_project3_FINAL3.txt", h=T,sep = "\t", stringsAsFactors = FALSE, quote = "") 
points_humans <- humans_final3[,c("cell_Longitude", "cell_Latitude", "Interval")]
bins3 <- as.vector(unique(humans_final3$Interval))
for (i in seq_along(bins3))
{
 points_humans <- subset(humans_final3, Interval==bins3[i])
 points_h <- points_humans[,c("cell_Longitude", "cell_Latitude")]
 cl_ras <- get(ls(pattern = paste("clim_", bins3[i], "k", sep="")))
 ext <- extract(cl_ras, points_h, cellnumbers=TRUE, df=TRUE, nl=8, na.rm=FALSE)
 hum_id <- cellFromXY(subset(cl_ras,1), points_h)
 cl_df <- na.omit(as.data.frame(cl_ras))
 pca_cl <- dudi.pca(cl_df, scannf = F, nf = 5) ### pca on all variables for contirbution
 fviz_pca_var(pca_cl, col.var = "contrib", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) ### plotting variable loadings
 res.var <- get_pca_var(pca_cl) ### results of pca 
 s.class(pca_cl$li[, 1:2], fac= factor(rownames(cl_df)%in% hum_id, levels = c("FALSE", "TRUE" ),labels = c("background", "Homo_Sapiens")), col=c("red", "blue"),
 csta = 0, cellipse = 2, cpoint = .6, pch = 19, sub=paste0("Asia_", bins3[i], "ka", sep=""))
 assign(paste0("Extracted_humans_", bins3[i], "k", sep=""), data.frame(points_humans, ext[,c(2:10)]))
 }
  alll <- ls(pattern="Extracted_humans_")[c(26:27,1:25)]
  alls <- lapply(alll, get)
  Human_all <- do.call(rbind, alls)
  
  Human_all2 <- Human_all[which(is.na(Human_all$Prec_Sp)),]
  rows <- vector(as.numeric(rownames(Human_all2)))
  final_humans <- Human_all[-rows,]
  
  write.table(final_humans, file = "Humans_project3_SDMs.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
 ### the less correlated variables(Prec_Su less so).
  vars < - c("Temp_Wi", "Prec_Au", "Temp_Su", "Prec_Su")  

##############################################################
# Extract climatic variables from georeferenced occurrences #
##############################################################


spec_final <- read.delim("Species_project3_FINAL.txt", h=T,sep = "\t", stringsAsFactors = FALSE, quote = "") 
species_list <- as.vector(unique(spec_final$Species))

for (sp in seq_along(species_list))
{
sps <- subset(spec_final, Species==species_list[sp])
points_species <- spec_final[,c("cell_Longitude", "cell_Latitude", "Interval")]
bins3 <- sort(as.vector(unique(sps$Interval)), decreasing=FALSE)
species_data <- as.data.frame(matrix(NA, ncol=ncol(spec_final) + 9))
colnames(species_data) <- c(colnames(spec_final), "cells", "Prec_Sp", "Prec_Su", "Prec_Au", "Prec_Wi", "Temp_Sp", "Temp_Su", "Temp_Au", "Temp_Wi")
for (i in seq_along(bins3)) ### Prec_Au, Temp_Su, Temp_Au ###
{
 points_sps <- subset(sps, Interval==bins3[i])
 points_s <- points_sps[,c("cell_Longitude", "cell_Latitude")]
 cl_ras <- get(ls(pattern = paste("clim_", bins3[i], "k", sep="")))
 ext <- extract(cl_ras, points_s, cellnumbers=TRUE, df=TRUE, nl=8, na.rm=FALSE)
 sub_sp <- data.frame(points_sps, ext[,c(2:10)])
 #spec_id <- cellFromXY(subset(cl_ras,1), points_s)
 cl_df <- na.omit(as.data.frame(cl_ras))
 species_data <- rbind(species_data, sub_sp)
}
assign(paste0("Extracted_", species_list[sp], sep=""), species_data[-c(1),])
}

  Sp_alll <- ls(pattern="Extracted_")
  Sp_alls <- lapply(Sp_alll, get)
  Species_all <- do.call(rbind, Sp_alls)
  Species_all2 <- Species_all[which(is.na(Species_all$Prec_Sp)),]
  rows <- as.vector(as.numeric(rownames(Species_all2)))
  final_sp <- Species_all[-rows,]
 
