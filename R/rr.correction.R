#' rr.correction
#'
#' daily RR correction. Based on DNMI KLIMA, Manual for operational correction of Nordic precipitation data, 24/96
#' P_cor<-exposition_coef*(P_mes+mean daily evaporation loss+wetting loss). All the constants come from the report
#' @param RR_values
#' @param RR_metadata
#' @param RR_covariate
#' @param TAMRR_values
#' @param TAMRR_metadata
#' @param TAMRR_covariate
#' @param dateTS
#' @keywords RR correction
#' @export
#' @examples
#' rr.correction()
rr.correction<-function(RR_values,
                        RR_metadata,
						RR_covariate,
					    TAMRR_values,
						TAMRR_metadata,
						TAMRR_covariate,
					    dateTS) {


   ####################################
   # INITIALIZING
   ####################################
   RR_cor<-matrix(NA,nr=nrow(RR_values),nc=ncol(RR_values))
   TAM_int<-matrix(NA,nr=nrow(RR_values),nc=ncol(RR_values))
 
   ######################################################
   # Mean daily evaporation loss according to the month
   ######################################################
   dataTS<-as.numeric(sub("-","",sub("-","",sub("UTC","",dateTS)))) #strptime(timeserie, "%Y-%m-%d" , tz = "UTC")
   year<-round(dataTS/1e4)
   month<-round((dataTS-year*1e4)/100)
   day<-round((dataTS-year*1e4-month*100))
   nbtimestep<-length(dateTS)
   

   ##########################################
   # CORRECTION
   ##########################################
   exposition_coef<-RR_metadata[,5:7]

   flag<-rep(NA,nbtimestep)
   for (i_t in 1: nbtimestep) {   
      temp<-try(rr.correction.foerland(RR_values=RR_values[i_t,],
                           RR_metadata=RR_metadata,
						   RR_covariate=RR_covariate,
					       TAMRR_values=TAMRR_values[i_t,],
						   TAMRR_metadata=TAMRR_metadata,
						   TAMRR_covariate=TAMRR_covariate,
						   month=month[i_t],
						   CV=FALSE
						   ),TRUE)
						   
      				   
      RR_cor[i_t,]<-temp$RR_values_cor 
      TAM_int[i_t,]<-temp$TAM_values_int
	  
	  flag[i_t]<-is.character(temp$RR_values_cor)
   }
   
   RR_cor[which(flag),]<-RR_values[which(flag),]
   TAM_int[which(flag),]<-TAM_int[which(flag),]
   date_No_correction<-dateTS[which(flag)]

   ##############################
   # RESULTS
   ##############################
   results<-list(cor=RR_cor,
                 TAM_int=TAM_int,
                 date_No_correction=date_No_correction)
   
   return(results)
   

}