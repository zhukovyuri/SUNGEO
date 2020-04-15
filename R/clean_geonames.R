#' Clean GeoNames
#'
#' Function to clean GeoNames gazetteer files.
#'
#' @param gn Raw GeoNames data file. \code{data.frame} or \code{data.table} object.
#' @param ppl Keep only populated place names? Logical.
#' @return Returns a \code{data.table} with properly-formatted GeoNames information.
#' @export
#' @import tidyverse data.table
#' @importFrom data.table last first between
#' @details Function used internally by \code{geocode_gn}

clean_geonames <- function(gn,ppl=FALSE){

  # Tab-delimited
  geonames_t <- apply(gn,1,function(x){paste(x,collapse="\t")})
  geonames_t <- as.character(geonames_t)

  # Fix incomplete rows
  nskips <- grep("\n",geonames_t,fixed=T)
  while(length(nskips)>0){print(length(nskips))
    wch <- nskips[1]
    subtext <- geonames_t[wch]
    subtext <- unlist(strsplit(subtext,"\n",fixed=T))
    geonames_t <- c(geonames_t[1:(wch-1)],subtext,geonames_t[(wch+1):length(geonames_t)])
    nskips <- grep("\n",geonames_t,fixed=T)
  }
  geonames_t <- strsplit(geonames_t,split="\t")
  colz <- c()
  for(i in 1:length(geonames_t)){
    colz[i] <- length(geonames_t[[i]])
  }
  geonames_new <- list()
  for(i in which(colz==19)){
    geonames_new[[i]] <- t(geonames_t[[i]])
  }
  geonames_t <- do.call(rbind,geonames_new) %>% as.data.frame()

  # Add column names
  names(geonames_t) <- c("geonameid","name","asciiname","alternatenames","latitude","longitude","feature_class","feature_code","country_code","cc2","admin1_code","admin2_code","admin3_code","admin4_code","population","elevation","dem","timezone","modification_date")
  classez <- c("numeric","character","character","character","numeric","numeric","character","character","character","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character")
  geonames_t <- suppressWarnings(
    geonames_t %>% dplyr::mutate_if(classez%in%"numeric",as.character)  %>% dplyr::mutate_if(classez%in%"numeric",as.numeric)
    )
  geonames_t <- suppressWarnings(
    geonames_t %>% dplyr::mutate_if(classez%in%"character",as.character)
  )

  # Output
  geonames <- geonames_t %>% as.data.table()
  if(ppl){geonames <- geonames[grepl("PPL",feature_code),]}
  return(geonames)
}
