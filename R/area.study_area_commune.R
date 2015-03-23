#' study_area_commune
#'
#' study_area_commune (Only in Norway so far)
#' @param shp shapefile
#' @param name
#' @param background.grid
#' @param cellsize_background
#' @keywords study area
#' @export
#' @examples
#' study_area_commune()
study_area_commune<-function(shp,
                    name,
					background.grid,
					cellsize_background
                    ) {
					
   indice<-which(shp@data$NAVN==name)
   coordinates<-shp@polygons[[indice]]@Polygons[[1]]@coords
   surface<-shp@data$AREA[indice]/1e6 #m^2 ->km^2
					
   temp_background<-as.data.frame(background.grid[,1:2])
   coordinates(temp_background)<- ~x+y
					
   #process
   temp<-as.data.frame(coordinates)
   temp2<-Polygons(list(Polygon(coordinates)),name)
   temp3<-SpatialPolygons(list(temp2))
   background<-overlay(temp3, temp_background)
   indice<-which(!is.na(background))
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
   # }

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
