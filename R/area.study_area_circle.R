#' study_area_circle
#'
#' study_area_circle
#' @param xCenter
#' @param yCenter
#' @param radius
#' @param name
#' @param background.grid
#' @param cellsize_background
#' @param Nx_background
#' @param Ny_background
#' @keywords study area
#' @export
#' @examples
#' study_area_circle()
study_area_circle<-function(xCenter,
                            yCenter,
							radius,
                            name,
					        background.grid,
					        cellsize_background,
							Nx_background,
							Ny_background
                            ) {
 
   #location of center of circle on the background grid
   indice_x_center<-round((xCenter-min(background.grid[,1],na.rm=TRUE))/cellsize_background)
   indice_y_center<-round((max(background.grid[,2],na.rm=TRUE)-yCenter)/cellsize_background)
   temp<-indice.mat2vec(indice_x_center,indice_y_center,Nx_background,Ny_background)
   indice_xy<-temp$indice_xy
   rm(temp)
   
   indice<-indice.within_a_circle(indice_i=indice_xy, #center of the circle
                          Nx=Nx_background,
						  Ny=Ny_background,
						  delta=(radius/cellsize_background),
						  delta_res=1) 
	
   area_name<-paste("study_area",name,sep="")
   xcoord<-background.grid[indice$indice_deltaxy,1]
   ycoord<-background.grid[indice$indice_deltaxy,2]
   xMin<-min(xcoord,na.rm=TRUE)
   xMax<-max(xcoord,na.rm=TRUE)
   yMin<-min(ycoord,na.rm=TRUE)
   yMax<-max(ycoord,na.rm=TRUE)
   Nx<-(xMax-xMin)/cellsize_background+1
   Ny<-(yMax-yMin)/cellsize_background+1
   
   #For later study, knowing the indice of each pixel within a rectangular grid is important. So computing indices in a rectangular grid
   #indice_xy<-rep(0,length(xcoord))
   #for (i in 1:length(xcoord)) {
   #   indice_x<-round((xcoord[i]-xMin)/cellsize_background)+1
   #   indice_y<-round((yMax-ycoord[i])/cellsize_background)+1
   #   temp<-indice.mat2vec(indice_x,indice_y,Nx,Ny)
   #   indice_xy[i]<-temp$indice_xy
   #}

   result<-list(area_name=area_name,
                grid=cbind(xcoord,ycoord),
				xMin=xMin,
				xMax=xMax,
				yMin=yMin,
				yMax=yMax,
				Nx=Nx,
				Ny=Ny,
				indice=indice)
						 
   return(result)						 
}
