#' tam.interpolation.point.senorge
#'
#' TAM SeNorge Interpolation method v1.1
#'  Tveito, O. E., E. Førland, R. heino, I. Hansen-Bauer, H. Alexandersson, B. Dahlström, A.     
#'           Drebs, C. Kern-Hansen, T. Jónsson, E. Vaarby Laursen and Y. Westman, 2000: “Nordic  
#'           temperature maps”, DNMI Report no. 09/00                                              
#'  Tveito, O. E., I. Bjørdal, A. O. Skjelvåg and B. Aune, 2005: ”A GIS-based agro-ecological    
#'            decision system based on gridded climatology, Met. Apps., Vol. 12: 1, p. 57-68     
#'                                                                                               
#'  Synthesis note:                                                                              
#'  Mohr M.  New Routines for Gridding of Temperature and Precipitation Observations for         
#'            “seNorge.no”, MET report No. 08/2008                                               
#'                                                                                               
#' @param TAM_metadata
#' @param TAM_values
#' @param TAM_covariate
#' @param target_metadata 
#' @param target_covariate
#' @param month
#' @param CV
#' @keywords TAM
#' @export
#' @examples
#' tam.interpolation.point.senorge()   
tam.interpolation.point.senorge<-function(TAM_metadata,
                                  TAM_values,
								  TAM_covariate,
                                  target_metadata,
								  target_covariate,
								  month,
								  CV) { #Cross-validation   
								  
   ####################
   # CONSTANTS
   ####################
   
   #Detrending constants
   k <- c(19.4879,14.7045,19.8752,26.0346,35.1963,39.4180,35.9240,36.8990,34.4318, 29.1849,25.9764,18.0301)
   v1 <-c(-.0012,-.0019, -.0046,-.0061,-.0063,-.0063,-.0061, -.0057,-.0055,-.0046, -.0032, -.0016)
   v2 <-c(-.0051,-.0043, -.0021, -.0006, 0.0007,0.0011, 0.0009,0.0000,-.0011,-.0016, -.0031, -.0041)
   v3 <-c(-.0083,-.0062,-.0033,-.0008, 0.0000,0.0021, 0.0016, 0.0005,0.0000,-.0017, -.0054, -.0089)
   v4 <-c(-.2694,-.1918,-.2579,-.3288, -.4166,-.4460,-.3700, -.3763,-.3764, -.3349,-.3491,-.2459)
   v5 <-c(-.4395,-.4282, -.3044, -.1505,-.0363,0.091, 0.1290,0.0370,-.0543,-.1581, -.2194,-.3470)
   
   #Fennoms from DEM used as regressors when de-trending TAM are the ones described in the MET report Mohr,2008. 

   # Variogram models 
   sill<-c(7,5,1.3,.33,.7,1,.6,.19,.3,1.3,3.5,7)
   nugget<-c(0,0,0,0,0,0,0,0,0,0,0,0)
   model<-c("Exp","Exp","Exp","Exp","Exp","Exp","Sph","Exp","Exp","Exp","Exp","Exp")
   range<-c(250,250,200,100,75,150,500,100,75,175,200,250)
  
   #########################
   # TAM DATASET
   #########################
   nr_meta<-dim(TAM_metadata)[1]
  
   TAM_metadata2<-data.frame(matrix(nrow=nr_meta,ncol=10))
   TAM_metadata2[ , 1:4 ]<-TAM_metadata
   TAM_metadata2[,5]<-TAM_covariate[,7] #$meanDEM
   TAM_metadata2[,6]<-TAM_covariate[,6] #$minDEM
   TAM_metadata2[,7]<-TAM_covariate[,3] #$Lat
   TAM_metadata2[,8]<-TAM_covariate[,4] #$Lon
   TAM_metadata2[,9]<-TAM_values
   colnames(  TAM_metadata2)<-c("Stnr","HoH","X","Y","meanDEM","minDEM","LAT","LON","TAM","TM0")
   
   
   #########################
   # Target DATASET
   #########################
   nr_meta<-dim(target_metadata)[1]
  
   target_metadata2<-data.frame(matrix(nrow=nr_meta,ncol=10))
   target_metadata2[ , 1:4 ]<-target_metadata[,1:4]
   target_metadata2[,5]<-target_covariate[,7] #$meanDEM
   target_metadata2[,6]<-target_covariate[,6] #$minDEM
   target_metadata2[,7]<-target_covariate[,3] #$Lat
   target_metadata2[,8]<-target_covariate[,4] #$Lon
   colnames(  target_metadata2)<-c("Stnr","HoH","X","Y","meanDEM","minDEM","LAT","LON","TAM","TM0")
   
   #########################
   # TAKING OUT THE TREND (elevation, lat, lon.,,,)
   #########################
  
   TAM_metadata2$TM0<- TAM_metadata2$TAM - k[month] - (v1[month]*TAM_metadata2$HoH) -
                                                            (v2[month]*TAM_metadata2$meanDEM) -
								                            (v3[month]*TAM_metadata2$minDEM) -
								                            (v4[month]*TAM_metadata2$LAT) -
								                            (v5[month]*TAM_metadata2$LON) 

   
   TAM_metadata2<-na.omit( TAM_metadata2)
	
   ##########################
   # KRIGING
   ##########################  
   #model    
   variog_model <-vgm(sill[month],model[month],range[month],nugget=0)
    
   # krige(formula, locations, data, newdata, model,
   TAM_metadata2$X<-TAM_metadata2$X/1000
   TAM_metadata2$Y<-TAM_metadata2$Y/1000
   coordinates(TAM_metadata2) = ~X+Y
   
   target_metadata2<-target_metadata2
   target_metadata2$X<-target_metadata2$X/1000
   target_metadata2$Y<-target_metadata2$Y/1000
   coordinates(target_metadata2) = ~X+Y

   #Kriging and Cross-Validation
   if (!CV) {
      TAM_cv<- NULL
	  TAM_k<-krige(TM0~1,TAM_metadata2,target_metadata2,model=variog_model,nmax=6)
   } else {
      TAM_cv<-krige.cv(TM0~1,TAM_metadata2,target_metadata2,model=variog_model)#,nmax=6) #nfold=10
	  TAM_k<-NULL
   }  

   #########################
   # PUTTING BACK THE TREND (elevation, lat, lon.,,,)
   ######################### 
   if (!CV) {
      target_metadata2$TAM <-  TAM_k$var1.pred + k[month] + (v1[month]*target_metadata2$HoH) +
                                                            (v2[month]*target_metadata2$meanDEM) +
								                            (v3[month]*target_metadata2$minDEM) +
								                            (v4[month]*target_metadata2$LAT) +
								                            (v5[month]*target_metadata2$LON)
   }
   
   ###############
   # RESULTS     
   ###############
   results<-list(krig=list(var1.pred=target_metadata2$TAM,
                           var1.var=TAM_k$var1.var),   
                 cv=TAM_cv)
				 
   return(results)				 

}
