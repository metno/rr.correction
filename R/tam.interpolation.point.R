#' tam.interpolation.point
#'
#' tam.interpolation.point                                                                       
#' @param TAM_metadata
#' @param TAM_values
#' @param TAM_covariate
#' @param target_metadata 
#' @param target_covariate
#' @param month
#' @param CV
#' @param method
#' @keywords TAM
#' @export
#' @examples
#' tam.interpolation.point() 
tam.interpolation.point<-function(TAM_metadata,
                                  TAM_values,
								  TAM_covariate,
                                  target_metadata,
								  target_covariate,
								  month,
								  CV, #Cross-validation
								  method="SeNorge") {  #v1.0
								  
								  
			
   if (method=="SeNorge") {
      tam.interpolation.point.senorge(TAM_metadata=TAM_metadata,
	                                 TAM_values=TAM_values,
									 TAM_covariate=TAM_covariate,
									 target_metadata=target_metadata,
									 target_covariate=target_covariate,
									 month=month,
									 CV=CV)
   } else { } # to be written
}