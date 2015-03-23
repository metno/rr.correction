#################################################
#################################################
##                                             ##
##  RR correction                              ##
##  The whole Norwegian territory (continent)  ##
##  Data from 1970-2011                        ##
##                                             ##
#################################################
################################################# 

rm(list=ls())

###########################################################################
###########################################################################
##                                                                       ##
##   LIBRARIES                                                           ##
##                                                                       ##
###########################################################################
###########################################################################
require(rgdal) 
require(gstat) 

library(rr.correction)


###########################################
###########################################
##  LOAD DATA:                           ##
##         - GIS                         ##
##         - DEM                         ##  
##         - figure                      ##
##         - QPE radar HURUM             ##
###########################################
###########################################

data(border_communes)
data(border_norway)
data(border_norway)
data(coordinates_norway)
data(coordinates_norway_coast)
data(dem.grid)
data(dem.parameter)
data(figure.parameter)


###########################################################################
###########################################################################
##                                                                       ##
##   PARAMETERS                                                          ##
##                                                                       ## 
###########################################################################
###########################################################################

#PATH
source_path<-"/path/"                       
path_result<-"/path/results/"
dir.create(path_result, showWarnings = FALSE)

# AREA #

case_study<-c("0","Norway")     # circular area over HURUM
name<-case_study[2]

							  
############################################################################
############################################################################
## AREAS OF THE CASES STUDIES                                             ##
##                                                                        ##
############################################################################
############################################################################

area_case_study<-study_area_norway(shp=border_norway,
                              name=case_study[2],
				              background.grid=dem.grid[,1:2],
				              cellsize_background=dem.parameter$cellsize 


							  
#################################################
#################################################
##                                             ##
##  Correction                                 ## 
##                                             ##
#################################################
################################################# 

data(TAMRR)
TAMRR<-gages
data(TAMRR.covariate)
TAMRR.covariate<-gages.covariate

data(RRNSDNUD)
RRNSDNUD<-gages
data(RRNSDNUD.covariate)
RRNSDNUD.covariate<-gages.covariate

RR_corrected<-rr.correction(RR_values=RRNSDNUD$values,
                                    RR_metadata=RR$metadata,
									RR_covariate=RR.covariate,
					                TAMRR_values=TAMRR$values,
									TAMRR_metadata=TAMRR$metadata,
									TAMRR_covariate=TAMRR.covariate,
									dateTS=RR$dateTS)

		   
#dir.create(paste(path_result, "data",sep=""), showWarnings = FALSE)
#save(gages,file= paste(path_result,"data/",variable_p,type_p,".RData",sep=""))

