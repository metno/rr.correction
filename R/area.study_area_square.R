#' study_area_square
#'
#' study_area_square
#' @param xMin
#' @param xMax
#' @param yMin
#' @param yMax
#' @param name
#' @param background.grid
#' @param cellsize_background
#' @keywords study area
#' @export
#' @examples
#' study_area_square()
study_area_square<-function(xMin,
                            xMax,
                            yMin,
                            yMax,
                            name,
					        background.grid,
					        cellsize_background
                            ) {
					
   indice<-which((background.grid[,1]>=xMin) & (background.grid[,1]<=xMax) &
	                                   (background.grid[,2]>=yMin) & (background.grid[,2]<=yMax)
                                       )
   area_name<-paste("study_area",name,sep="")
   xcoord<-background.grid[indice,1]
   ycoord<-background.grid[indice,2]
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
