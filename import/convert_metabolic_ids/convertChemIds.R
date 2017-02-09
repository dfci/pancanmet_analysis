#' Small Molecule Chemical ID Converter
#'
#' Converts an ID from one ID type to another using the 
#' Chemical Translation Service (http://cts.fiehnlab.ucdavis.edu/)
#'
#' @param id small molecule ID
#' @param fromIdType ID type being converted "from" from this list: 
#'   http://cts.fiehnlab.ucdavis.edu/service/convert/fromValues
#' @param toIdType ID type being converted "to" from this list: 
#'   http://cts.fiehnlab.ucdavis.edu/service/convert/toValues
#' @return a named vector of IDs 
#' @author Augustin Luna (lunaa@cbio.mskcc.org)
#' @example
#'   convertChemIds("6305", "PubChem CID", "ChEBI")
convertChemIds <- function(id, fromIdType, toIdType, debug=TRUE) {
	require(RCurl)
	require(stringr) 
	require(rjson) 
	
	toId <- NULL
	
	if(debug) {
		cat("ID: ", id, "\n")
	}

	#Example: http://cts.fiehnlab.ucdavis.edu/service/convert/PubChem%20CID/ChEBI/6305
	url <- paste("http://cts.fiehnlab.ucdavis.edu/service/convert/", fromIdType, "/", toIdType, "/", id, sep="")
	url <- URLencode(url)

	if(debug) {		
		cat("URL: ", url, "\n")
	}	

	results <- getURL(url)
	results <- fromJSON(results)
	
	toId <- results[[1]]$result
	
	if(is.null(toId) | length(toId) == 0) {
		return(NA)
	} else {
		return(toId) 
	}
}

