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
  #Climate2 <- Climate[which(Climate$Longitude < -168 | Climate$Longitude >= 60),]
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
    crs(clim_rasts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" ### find in one of our scripts the coordinate system we want to use
  }
  assign(paste("clim_", Interval,"k", sep=""), clim_rasts)
}
