#' geography.covariates
#'
#' geography covariates
#' @param dem.grid,
#' @param Nx_dem,  
#' @param Ny_dem,
#' @param cellsize_dem,
#' @param coverage.grid,
#' @param cellsize_coverage,
#' @param UTM,
#' @param border, 
#' @param coordinates_coast, 
#' @param coordinates_border,
#' @param FIGURE=TRUE,
#' @param path
#' @keywords geography
#' @export
#' @examples
#' geography.covariates()
geography.covariates<-function(dem.grid,
                        Nx_dem,  
						Ny_dem,
						cellsize_dem,
					    coverage.grid,
						cellsize_coverage,
						UTM,
                        border, 
						coordinates_coast, 
						coordinates_border,
						FIGURE=TRUE,
						path){
						
####################
# INITIALIZING
####################
   xMax_dem<-max(dem.grid[,1],na.rm=TRUE)
   xMin_dem<-min(dem.grid[,1],na.rm=TRUE)
   yMax_dem<-max(dem.grid[,2],na.rm=TRUE)
   yMin_dem<-min(dem.grid[,2],na.rm=TRUE)
   
   #focus on an assumed square coverage
   xMax_coverage<-max(coverage.grid[,1],na.rm=TRUE)
   xMin_coverage<-min(coverage.grid[,1],na.rm=TRUE)
   yMax_coverage<-max(coverage.grid[,2],na.rm=TRUE)
   yMin_coverage<-min(coverage.grid[,2],na.rm=TRUE)
   Nx_coverage<-(xMax_coverage-xMin_coverage)/cellsize_coverage
   Ny_coverage<-(yMax_coverage-yMin_coverage)/cellsize_coverage
   
   #square shape indice on the dem grid
   indice_square <- which( (dem.grid[,1]<=(xMax_coverage+cellsize_coverage)) &
                           (dem.grid[,1]>=xMin_coverage) &
						   (dem.grid[,2]<=(yMax_coverage+cellsize_coverage)) &
						   (dem.grid[,2]>=yMin_coverage ) )

   # if not a square grid						    
   if ( (Nx_coverage*Ny_coverage)>nrow(coverage.grid)) {
      #coverage indice on the dem grid
      indice_xy<-rep(0,nrow(coverage.grid))
      for (i in 1:nrow(coverage.grid)) {
         indice_x<-round((coverage.grid[i,1]-xMin_dem)/cellsize_dem)+1
         indice_y<-round((yMax_dem-coverage.grid[i,2])/cellsize_dem)+1
	     temp<-indice.mat2vec(indice_x,indice_y,Nx_dem,Ny_dem)
	     indice_xy[i]<-temp$indice_xy
       }
	   indice_coverage<-indice_xy
	   
   } else indice_coverage<-indice_square   
					

					   
# geography.covariate x,y,distance to the coast
geography.covariate<-matrix(NA,nc=3,nr=length(indice_coverage))
geography.covariate[,1]<-dem.grid[indice_coverage,1] #x
geography.covariate[,2]<-dem.grid[indice_coverage,2] #y
#dem.covariate[,3]                             #distance from the coast


##########################
# Distance to the coast  #
##########################
for (i_pixel in 1: length(indice_coverage)) {
  
   tmp<-(dem.grid[indice_coverage[i_pixel],1]-coordinates_coast[,1])^2
   tmp2<-(dem.grid[indice_coverage[i_pixel],2]-coordinates_coast[,2])^2
   tmp3<-sqrt(tmp+tmp2)
   geography.covariate[i_pixel,3]<-tmp3[which(tmp3==min(tmp3,na.rm=TRUE))[1]]
}
#putting pixel located in the ocean to NA
geography.covariate[which(is.na(dem.grid[indice_coverage,3])),3]<-NA


#######################
# PLOT
#######################
if (FIGURE) { 

dir.create(paste(path, "figures",sep=""), showWarnings = FALSE)
xMin<-min(dem.covariate[,1],na.rm=TRUE)
xMax<-max(dem.covariate[,1],na.rm=TRUE)
yMin<-min(dem.covariate[,2],na.rm=TRUE)
yMax<-max(dem.covariate[,2],na.rm=TRUE)

# Distance to the coast
Min_diff<-min(geography.covariate[,3],na.rm=TRUE)
Max_diff<-max(geography.covariate[,3],na.rm=TRUE)
colour_indice<-round(100*(geography.covariate[,3]-Min_diff)/(Max_diff-Min_diff))/100
colour_indice[is.na(colour_indice)]<-0
colour_indice<-1-colour_indice
colour_indice<-10^colour_indice/10
 #colour_indice<-log(colour_indice)
colour_indice<-grey(colour_indice)
#colcode<-heat.colors(colour_indice)
tiff(filename = paste(path,"figures/distance to the coast.tif",sep=""), width = width_mm, height = height_mm,
     units = "mm", pointsize = 12,
     compression = "lzw",
     bg = "transparent", res = 300,
     restoreConsole = TRUE) #col=colcode[dem.grid[,3]],
par(fig=c(0,0.8,0,1))
plot(geography.covariate[,1],geography.covariate[,2],col=colour_indice,pch=15,cex=0.4,xlim=c(xMin,xMax),ylim=c(yMin,yMax),axes=FALSE,xlab="", ylab="",asp=1)
par(new=TRUE)
plot(border_norway,border="blue",lwd=1,axes=FALSE,xlim=c(xMin,xMax),ylim=c(yMin,yMax),asp=1)
box() 
axis(2, ylim=c(yMin,yMax))
mtext("y UTM33 (meters)",side=2,line=2.5)	 
axis(1, ylim=c(xMin,xMax))
mtext("x UTM33 (meters)",side=1,line=2.5)
#legend
max_legend<-Max_diff
min_legend<-Min_diff
par(fig=c(0.7,1,0,1),new=TRUE)
i1 <- seq(min_legend,(max_legend-abs(min_legend-max_legend)/length(colour_indice)),abs(min_legend-max_legend)/length(colour_indice))
i2 <- seq((min_legend+abs(min_legend-max_legend)/length(colour_indice)),max_legend,abs(min_legend-max_legend)/length(colour_indice))
min_length<-min(length(i1),length(i2),na.rm=TRUE)
plot(i1[1:min_length],i2[1:min_length],type = "n", xlim=c(min_legend,max_legend),ylim=c(min_legend,max_legend),axes=FALSE,xlab="", ylab="")
rect(min(i1), i1,max(i1), i2,col=grey(seq(1,0,-1/length(i1))),border="transparent")
axis(4, ylim=c(0,max(i1,na.rm=TRUE)))
#title("mm",line=2.5) #,cex.main = 0.5)
dev.off()


}

results<-list(geography.covariate=geography.covariate)
return(results)

}
	
