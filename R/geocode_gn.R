#' Geocode addresses with GeoNames
#'
#' Function to find geographic coordinates of addresses and place names, using GeoNames gazetteer data.
#'
#' @param query Address(es) or place name(s) to be geocoded. Character string or vector.
#' @param country_name Name of country. Must supply either \code{country_name} or \code{country_iso3}. Character string.
#' @param country_iso3 ISO 3166-1 alpha-3 country code. Must supply either \code{country_name} or \code{country_iso3}. Character string.
#' @param str_meth String distance metric (\code{"osa"}, \code{"lv"}, \code{"dl"}, \code{"hamming"}, \code{"lcs"}, \code{"qgram"}, \code{"cosine"}, \code{"jaccard"}, \code{"jw"}, \code{"soundex"}). Default is \code{"cosine"}. See \link[stringdist]{stringdist-metrics}. Character string.
#' @param match_all Attempt to match all queries? Default is \code{FALSE}. Logical.
#' @param details Should detailed results be returned? Default is \code{FALSE}. Logical.
#' @param ncores Number of cores to use for parallel socket cluster. Default is 1. Numeric.
#' @param verbose Print status messages and progress bar? Default is \code{FALSE}. Logical.
#' @return A \code{data.table} object. If \code{details=FALSE}, contains fields
#' \itemize{
##'  \item{"query". }{User-supplied address query(ies). Character string.}
##'  \item{"geonameid". }{GeoNames ID. Character string.}
##'  \item{"address". }{GeoNames place name. Character string.}
##'  \item{"longitude". }{Horizontal coordinate. Numeric.}
##'  \item{"latitude". }{Vertical coordinate. Numeric.}
##'  }
##' If \code{details=TRUE}, contains additional fields
#' \itemize{
##'  \item{"country_iso3". }{ISO 3166-1 alpha-3 country code. Character string.}
##'  \item{"country_name". }{Country name Character string.}
##'  \item{"strdist". }{String distance. Numeric.}
##'  \item{"str_meth". }{String distance metric. Character string.}
##'  \item{"match_type". }{Type of match. Character string.}
##'  }
#' @export
#' @import tidyverse data.table countrycode stringdist parallel
#' @importFrom data.table last first between
#' @importFrom rvest html_session
#' @details This function geocodes one or multiple addresses by computing string distances between the user-supplied address query(ies) and place names (including alternate, historical and non-English spellings) in the GeoNames' free place name gazetteer data (\url{http://www.geonames.org/}). Returns geographic coordinates and other information for the matched place name. Three types of string matches, along with
#' \itemize{
##'  \item{"perfect". }{Perfect string match (e.g. "Frankfurt am Mein" -- "Frankfurt am Mein").}
##'  \item{"partial". }{Partial string match (e.g. "Frankfurt" -- "Frankfurt am Mein").}
##'  \item{"fuzzy". }{Approximate string match (e.g. "Frnkfurt" -- "Frankfurt am Mein").}
##'  }
#' Note that -- while capabale of geocoding a single address at a time -- this function is optimized for batch geocoding of multiple addresses at once. For single addresses, we recommend using \code{\link[SUNGEO]{geocode_osm}}.
#'
#' @examples
#' # Geocode an address
#' \dontrun{
#' geocode_gn("Chisinau", country_name = "Moldova")
#' }
#' # Return detailed results
#' \dontrun{
#' geocode_gn("Chisinau", country_name = "Moldova", details = TRUE)
#' }
#' # Return detailed results for multiple addresses, with progress reports
#' \dontrun{
#' geocode_gn(query = c("Chisinau","Buiucani, Chisinau","Chisinau centru"),
#'            country_name = "Moldova", details = TRUE, verbose = TRUE)
#' }

geocode_gn <- function(query,
                       country_name=NULL,
                       country_iso3=NULL,
                       str_meth="cosine",
                       match_all=FALSE,
                       details=FALSE,
                       ncores=1,
                       verbose=FALSE){

  # Country codes
  if(length(country_name)==0&length(country_iso3)==0){stop("Please supply country code or name.")}
  if(length(country_name)>0&length(country_iso3)==0){
    country_iso3 <- countrycode::countrycode(country_name,"country.name","iso3c")
    country_iso2 <- countrycode::countrycode(country_name,"country.name","iso2c")
  }
  if(length(country_name)==0&length(country_iso3)>0){
    country_name <- countrycode::countrycode(country_iso3,"iso3c","country.name")
    country_iso2 <- countrycode::countrycode(country_iso3,"iso3c","iso2c")
  }
  if(length(country_name)>0&length(country_iso3)>0){
    country_iso2 <- countrycode::countrycode(country_iso3,"iso3c","iso2c")
  }

  # Download and parse GeoNames (with error handling)
  downloadFail <- FALSE
  tryCatch({
    if(verbose){paste0("Fetching gazetteer data for ",country_name,"...") %>% print()}
    # Download and unzip
    tempf <- tempfile()
    download.file(paste0("http://download.geonames.org/export/dump/",country_iso2,".zip"),destfile=tempf,cacheOK = TRUE,quiet = (!verbose), method = ifelse(Sys.info()['sysname']!="Windows","auto","wb") )
    con_gn <- unzip(zipfile = tempf,files = paste0(country_iso2,".txt"))
    gn <- read.delim(con_gn[grep(paste0(country_iso2,".txt"),con_gn)],sep="\t",header=F,quote = "", stringsAsFactors = FALSE,encoding = "utf-8") %>% as.data.table()
    gn <- clean_geonames(gn,ppl=TRUE)
    gn[,alternatenames := gsub(",","|",alternatenames)]
    gn[,TERM := paste0(name,"|",asciiname,"|",alternatenames,collapse="|") %>% gsub("\\|$","",.), by = geonameid]
    con_gn %>% textConnection() %>% close() %>% on.exit()
    unlink(tempf)
    gc()
    if(file.exists(con_gn)){file.remove(con_gn)}
  }, warning = function(w) {
    downloadFail <<- TRUE
  }, error = function(e) {
    message("Cannot access GeoNames server. Please check your internet connection and try again.")
    downloadFail <<- TRUE
  }, finally = {
  })
  if(downloadFail){
    cat("Cannot access GeoNames server. Please check your internet connection and try again.")
    return()
  } else {

    # Progress report
    if(verbose){paste0("Starting string distance matching...") %>% print()}

    # Reverse string distance matrix (parallel)
    if(match_all){
      cl <- parallel::makePSOCKcluster(ncores)
      parallel::setDefaultCluster(cl)
      parallel::clusterExport(NULL,c("query","str_meth","country_iso3","gn"),envir = environment())
      parallel::clusterEvalQ(NULL, expr=library(stringdist))
      parallel::clusterEvalQ(NULL, expr=library(dplyr))
      gn_list <- parallel::parLapply(NULL,1:nrow(gn),function(j){
        # cat(paste0(country_iso3,": ",j,"/",nrow(gn),"\r"))
        stringdist::stringdistmatrix(query %>% tolower(),strsplit(gn$TERM[j],"\\|")[[1]] %>% gsub("[^[:alnum:] ]", "",.) %>% tolower, method=str_meth, nthread=1) %>% apply(1,min)
      })
      parallel::stopCluster(cl)
      gc()
      # Drop null elements
      gn_mat_ <- data.frame(query=query,stringsAsFactors = FALSE) %>% cbind(gn[gn_list %>% bind_cols() %>% apply(1,which.min),.(geonameid,feature_code,asciiname,longitude,latitude)]) %>% as.data.table()
      gn_mat_[,strdist := gn_list %>% bind_cols() %>% apply(1,min)]
      gn_mat_[,str_met := str_meth]
    }
    gc()

    # Loop over ppl's
    cl <- parallel::makePSOCKcluster(ncores)
    parallel::setDefaultCluster(cl)
    parallel::clusterExport(NULL, c("query","str_meth","country_iso3","country_iso2","country_name","gn"),envir = environment())
    parallel::clusterEvalQ(NULL, library(stringdist))
    parallel::clusterEvalQ(NULL, library(dplyr))
    parallel::clusterEvalQ(NULL, library(data.table))
    gn_mat <- parallel::parLapply(NULL,1:nrow(gn),function(j){
      cat(paste0(country_iso3,": ",j,"/",nrow(gn),"\r"))
      # Create holder matrix
      gn_j <- data.frame(country_iso3=country_iso3,country_name=country_name,geonameid=gn$geonameid[j],gn_type=gn$feature_code[j],address=gn$asciiname[j],longitude=gn$longitude[j],latitude=gn$latitude[j],gn_perfect=NA_character_,gn_pt=NA_character_,gn_fz=NA_character_,strdist=NA_real_,str_meth=str_meth,stringsAsFactors = FALSE) %>% as.data.table()
      tryCatch({
        # Partial and fuzzy matches
        match_perfect <- which((query %>% tolower()) %in% (strsplit(gn$TERM[j],"\\|")[[1]] %>% tolower()))
        match_pt <- grep(gn$TERM[j] %>% tolower(),query %>% tolower())
        match_fz <- stringdist::stringdist(gn$TERM[j] %>% tolower(),query %>% tolower(),method=str_meth,nthread=1)
        # Take closest fuzzy match & exact match (if any)
        gn_j[,strdist := match_fz[match_fz %>% which.min()]]
        gn_j[,gn_fz := query[match_fz %>% which.min()]]
        if(length(match_perfect)>0){
          gn_j[,gn_perfect := query[match_perfect[1]]]
        }
        if(length(match_pt)>0){
          gn_j[,gn_pt := query[match_pt[1]]]
        }
      }, error=function(e){})
      gn_j
    }) %>% dplyr::bind_rows()
    parallel::stopCluster(cl)
    gc()

    # Progress report
    if(verbose){paste0("Extracting best matches...") %>% print()}

    # Consolidate
    gn_mat$ix <- 1:nrow(gn_mat)
    if(match_all){gn_mat_$ix <- 1:nrow(gn_mat_)}

    # Extract best matches
    gn_perfect <- gn_mat[!is.na(gn_perfect),][match(query,gn_perfect),] %>% as.data.frame()
    gn_pt <- gn_mat[!is.na(gn_pt),][match(query,gn_pt),] %>% as.data.frame()
    gn_fz <- gn_mat[ix%in%gn_mat[,ix[which.min(strdist)],by=gn_fz][,V1],][match(query,gn_fz),] %>% as.data.frame()

    # Replace missing exact matches with best fuzzy matches
    gn_best <- gn_perfect
    gn_best[is.na(gn_perfect$gn_pt),] <- gn_pt[is.na(gn_perfect$gn_pt),]
    gn_best[is.na(gn_pt$gn_pt)&is.na(gn_pt$gn_perfect),] <- gn_fz[is.na(gn_pt$gn_pt)&is.na(gn_pt$gn_perfect),]
    gn_best <- gn_best %>% as.data.table()
    gn_best[,query := query]
    gn_best[!is.na(gn_perfect),match_type:="perfect"]
    gn_best[!is.na(gn_pt)&is.na(gn_perfect),match_type:="partial"]
    gn_best[!is.na(gn_fz)&is.na(gn_pt)&is.na(gn_perfect),match_type:="fuzzy"]
    gn_best[, c("gn_pt","gn_fz","gn_perfect","ix"):=NULL]
    gn_best[match_type=="perfect",strdist := 0]
    gn_best[match_type=="perfect",str_meth := "identity"]
    gn_best <- gn_best %>% dplyr::select(query,geonameid,country_iso3,country_name,address,everything())

    # Grab missing from reverse matrix
    if(match_all){
      if(gn_mat_[query%in%gn_best[is.na(geonameid),query] ,] %>% nrow() > 0){
        gn_best[is.na(geonameid),query] == gn_mat_$query
        gn_best[is.na(geonameid),strdist := gn_mat_$strdist]
        gn_best[is.na(geonameid),str_meth := gn_mat_$str_meth]
        gn_best[is.na(geonameid),address := gn_mat_$asciiname]
        gn_best[is.na(geonameid),longitude := gn_mat_$longitude]
        gn_best[is.na(geonameid),latitude := gn_mat_$latitude]
        gn_best[is.na(geonameid),gn_type := gn_mat_$feature_code]
        gn_best[is.na(geonameid),match_type := "fuzzy"]
        gn_best[is.na(geonameid),geonameid := gn_mat_$geonameid]
      }}

    # Output
    if(!details){gn_best <- gn_best %>% select(query,geonameid,address,longitude,latitude)}
    return(gn_best)
  }

}
