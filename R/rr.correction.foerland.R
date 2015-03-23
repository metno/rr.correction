#' rr.correction.foerland
#'
#' daily RR correction. Based on DNMI KLIMA, Manual for operational correction of Nordic precipitation data, 24/96
#' P_cor<-exposition_coef*(P_mes+mean daily evaporation loss+wetting loss). All the constants come from the report
#' @param RR_values
#' @param RR_metadata
#' @param RR_covariate
#' @param TAMRR_values
#' @param TAMRR_metadata
#' @param TAMRR_covariate
#' @param month
#' @param CV  cross validation default FALSE
#' @keywords RR correction
#' @export
#' @examples
#' rr.correction.foerland()
rr.correction.foerland<-function(RR_values,
                     RR_metadata,
					 RR_covariate,
					 TAMRR_values,
					 TAMRR_metadata,
					 TAMRR_covariate,
					 month=month,
					 CV=FALSE
                    ) {
					
				
					
   ###################
   # INTITIALIZING
   ###################
   RR_values_cor<-rep(NA,length(RR_values))
   
   ####################
   # CONSTANTS
   ####################   
   #mean daily evaporation loss (mm/day) for Norway (DNMI KLIMA, Manual for operational correction of Nordic precipitation data, 24/96 p38, table 5.2)
   mean_daily_evaporation_loss<-vector(length=12) #12 values for 12 months)
   mean_daily_evaporation_loss[1]<-0.02 
   mean_daily_evaporation_loss[2]<-0.02
   mean_daily_evaporation_loss[3]<-0.03
   mean_daily_evaporation_loss[4]<-0.16
   mean_daily_evaporation_loss[5]<-0.04
   mean_daily_evaporation_loss[6]<-0.06
   mean_daily_evaporation_loss[7]<-0.06
   mean_daily_evaporation_loss[8]<-0.05
   mean_daily_evaporation_loss[9]<-0.03
   mean_daily_evaporation_loss[10]<-0.02
   mean_daily_evaporation_loss[11]<-0.02
   mean_daily_evaporation_loss[12]<-0.02

   #wetting loss (mm/case) (p38, table 5.2)
   #Four rows : rain, drizzle,snow,mixed
   #two columns: Norwegian gages (manual gage), weighting and tipping bucket
   # !! DRIZZLE ISN'T TAKEN INTO ACCOUNT!!
   # snow: T<0
   # mixed : 0<=T<2
   # precipitation: T>2
   wetting_loss<-matrix(nc=4,nr=2) 
   wetting_loss[1,]<-c(0.15,0.14,0.05,0.13)  #Norwegian gages, manual gage  
   wetting_loss[2,]<-c(0.15,0.15,0.1,0.15)   #weighting and tipping bucket

   ######################################
   # Because not all                    #
   # RR stations have thermometer       #
   # AND                                #
   # Because not all                    #
   # TAMRR stations get a measure       #   
   ######################################
   #RR stations with thermometer
   temp<-rbind(TAMRR_metadata[,c("X","Y")],RR_metadata[,c("X","Y")])
   indice_notarget<-which(duplicated(temp)) -dim(TAMRR_metadata)[1]
    
   indice_target<-seq(1,dim(RR_metadata)[1])	
   indice_target2<-seq(1,dim(RR_metadata)[1])
   if (length(indice_notarget)>0) {
      #RR stations without thermometer
      indice_target<-indice_target[-indice_notarget]
      
      #TAMRR station with rain gage but with NA TAMRR value
      temp<-rbind(RR_metadata[indice_notarget,c("X","Y")],TAMRR_metadata[,c("X","Y")])
      indice_tamrr<-which(duplicated(temp)) -length(indice_notarget)
      if (length(indice_tamrr)>0) {
	     indice_na<-which(is.na(TAMRR_values[indice_tamrr]))
		 indice_nona<-which(!is.na(TAMRR_values[indice_tamrr]))
	  }	 
      if (length(indice_na)>0) indice_target2<-indice_target2[indice_notarget[indice_na]]
	  
	  indice_target3<-c(indice_target,indice_target2)	  
	  indice_tamrr<-indice_tamrr[indice_nona]
	  
   } else {
      indice_target3<-indice_target
	  indice_tamrr<-c()
   }
   
   ###########################
   # TEMPERATURE INTERPOLATION
   ###########################
   # Because temperature is not always on precip gages
   # Temperature is interpolated on these locations
   
   TAM_interpolated<-tam.interpolation.point(TAM_metadata=TAMRR_metadata,
                                                TAM_values=TAMRR_values,
												TAM_covariate=TAMRR_covariate,
                                                target_metadata=RR_metadata[indice_target3,],
												target_covariate=RR_covariate[indice_target3,],
								                month=month,
												CV=CV,
												method="SeNorge")
	
   ###################################
   # TAMRR field used for correction #
   ###################################
   TAMRR_rr<-rep(NA,dim(RR_metadata)[1])
   indice_rr<-seq(1,dim(RR_metadata)[1])
   indice_rr<-indice_rr[-indice_target3]
   if (length(indice_tamrr)>0) TAMRR_rr[indice_rr]<-TAMRR_values[indice_tamrr]
   TAMRR_rr[indice_target3]<-TAM_interpolated$krig$var1.pred
   
   
   #################
   # CORRECTION    #
   #################   
   indice_snow<-which(TAMRR_rr<=0)	                 # snow										
   indice_mix<-which( (TAMRR_rr>0) & (TAMRR_rr<=2) ) # mix
   indice_precipitation<-which(TAMRR_rr>2)           # precipitation
   indice_dry<-which(RR_values==0)                   # dry, not to be forgotten. Not written in Førland report
	  
   #snow correction
   RR_values_cor[indice_snow]<-RR_metadata$k_sol[indice_snow]*(RR_values[indice_snow]+mean_daily_evaporation_loss[month]+wetting_loss[1,3]) 
	  
   #mix correction
   RR_values_cor[indice_mix]<-RR_metadata$k_liq[indice_mix]*(RR_values[indice_mix]+mean_daily_evaporation_loss[month]+wetting_loss[1,4]) 

   #precipitation correction
   RR_values_cor[indice_precipitation]<-RR_metadata$k_liq[indice_precipitation]*(RR_values[indice_precipitation]+mean_daily_evaporation_loss[month]+wetting_loss[1,1]) 
 
   #No exposition coef
   RR_values_cor[is.na(RR_metadata$exp_class)]<-RR_values[is.na(RR_metadata$exp_class)]
   
   #Dry 
   RR_values_cor[indice_dry]<-0

   #if CROSS VALIDATION CV
   #TAM_interpolated$cv$residual
   
   ################
   # RESULTS
   ################
   results<-list(RR_values_cor=RR_values_cor,
				 TAM_values_int=TAMRR_rr
				 )
   
   return(results)
	  
}
