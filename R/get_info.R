#' Information on available SUNGEO data files
#'
#' This function reports the availability of data files on the SUNGEO server, searchable by country and topic.
#'
#' @param country_names Country name(s). Character string (single country) or vector of character strings (multiple countries).
#' @param country_iso3s Country code (ISO 3166-1 alpha-3). Character string (single country) or vector of character strings (multiple countries).
#' @param topics Data topics. See \code{get_info()} for full list. Character string (single topic) or vector of character strings (multiple topics).
#' @return list object, with three slots: 'summary', 'topics', and 'geoset'.
#' @importFrom data.table data.table copy
#' @importFrom stringr str_split
#' @importFrom utils data
#' @examples
#' # Get list of all available data
#' \dontrun{
#' out_1 <- get_info()
#' out_1["summary"]
#' out_1["topics"]
#' out_1["geosets"]
#' }
#'
#' # Get list of available data for a single country
#' \dontrun{
#' out_2 <- get_info(country_names="Afghanistan")
#' out_2
#' }
#'
#' # Get list of available data for a single topic
#' \dontrun{
#' out_3 <- get_info(topics="Elections:LowerHouse:CLEA")
#' out_3
#' }
#'
#' # Get list of available data for a multiple countries and topics
#' \dontrun{
#' out_4 <- get_info(
#'                  country_names=c("Afghanistan","Zambia"),
#'                  topics=c("Elections:LowerHouse:CLEA","Events:PoliticalViolence:GED"))
#' out_4
#' }
#' @seealso \code{\link{get_data}}
#' @export
get_info <- function(
    country_names=NULL,
    country_iso3s=NULL,
    topics=NULL
){

  # Get object
  #utils::data("available_data")
  country_iso3 <- desc <- geoset_years <- geosets <- space_units <- year_range <- country_name_alt <- NULL

  # Country codes
  if(length(country_iso3s) == 0 & length(country_names) > 0){
    country_iso3s <- unique(SUNGEO::cc_dict[match(country_names,SUNGEO::cc_dict[,country_name_alt]),country_iso3])
  }

  # Default behavior
  if(is.null(country_iso3s)&is.null(topics)){
    avail_subset <- data.table::copy(SUNGEO::available_data)
    # Summarize
    country_by_topic <- sapply(sapply(avail_subset,function(x){x[,country_iso3]}),function(x){paste0(x,collapse=", ")})
    units_by_topic <- sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){paste0(x,collapse=", ")})
    years_by_topic <- sapply(sapply(sapply(avail_subset[!grepl("^Geoset",names(avail_subset))],function(x){x[,year_range]}),unique),function(x){paste0(x,collapse=", ")})
    years_by_geoset <- sapply(sapply(avail_subset[grepl("^Geoset",names(avail_subset))],function(x){x[,geoset_years]}),unique)
    years <- unlist(stringr::str_split(years_by_topic,", |-"))
    years <- paste0(range(years[is.finite(as.numeric(years))]),collapse="-")
    n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
    n_units <- length(unique(unlist(sapply(stringr::str_split(units_by_topic,", "),length))))
    n_geosets <- length(grep("^Geoset",names(avail_subset)))
    avail_geosets <- names(avail_subset)[grep("^Geoset",names(avail_subset))]
    avail_geosets <- data.table::data.table(geosets=avail_geosets)[,country_iso3:=country_by_topic[match(geosets,names(country_by_topic))]][,desc:=paste0(geosets," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(geosets," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    n_topics <- sum(!grepl("^Geoset",names(avail_subset)))
    avail_topics <- names(avail_subset)[!grepl("^Geoset",names(avail_subset))]
    avail_topics <- data.table::data.table(topics=avail_topics)[,country_iso3:=country_by_topic[match(topics,names(country_by_topic))]][,desc:=paste0(topics," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(topics," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    # Output
    summary_message <- c("-----------------------------------------","Summary of data available through SUNGEO","-----------------------------------------","",paste0("Found data for ",n_countries," countries, ",n_topics," topics"),paste0(n_units," spatial units, ",n_geosets," geographic boundary sets"),"","-------","Topics","-------",avail_topics[,(desc)],"","-------","Geosets","-------",avail_geosets[,(desc)],"","-------","See 'topics' and 'geosets' slots for details on spatial units and years.","Please narrow your search criteria for more informative results.","-------"
    )
    print(summary_message)
    out <- list(summary=summary_message,topics=avail_subset[!grepl("^Geoset",names(avail_subset))],geosets=avail_subset[grepl("^Geoset",names(avail_subset))])
  }

  # By topic
  if(is.null(country_iso3s)&!is.null(topics)){
    # Take subset
    avail_subset <- SUNGEO::available_data[topics]
    # Summarize
    if(length(topics)==1){
      country_by_topic <- paste0(sapply(avail_subset,function(x){x[,country_iso3]}),collapse=", ")
      names(country_by_topic) <- topics
      units_by_topic <- paste0(sort(unique(unlist(sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){stringr::str_split(x,", ")})))),collapse=", ")
      names(units_by_topic) <- topics
      years_by_topic <- paste0(range(unlist(sapply(sapply(avail_subset,function(x){x[,year_range]}),function(x){stringr::str_split(x,", |-")}))),collapse="-")
      names(years_by_topic) <- topics
      years <- paste0(range(unlist(stringr::str_split(years_by_topic,", |-"))),collapse="-")
      n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
      n_units <- length(unique(unlist(unique(stringr::str_split(units_by_topic,", ")))))
      n_geosets <- length(unique(unlist(sapply(sapply(sapply(avail_subset,function(x){x[,geosets]}),function(x){stringr::str_split(as.character(x),"-|, ")}),unlist))))
    }
    if(length(topics)>1){
      country_by_topic <- sapply(sapply(avail_subset,function(x){x[,country_iso3]}),function(x){paste0(x,collapse=", ")})
      names(country_by_topic) <- topics
      units_by_topic <- sapply(sapply(sapply(sapply(sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){stringr::str_split(x,", ")}),unlist),unique),sort),function(x){paste0(x,collapse=", ")})
      years_by_topic <- apply(sapply(sapply(sapply(avail_subset,function(x){x[,year_range]}),function(x){stringr::str_split(as.character(x),"-|, ")}),range),2,function(x){paste0(unique(x),collapse="-")})
      years <- paste0(range(unlist(stringr::str_split(years_by_topic,", |-"))),collapse="-")
      n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
      n_units <- length(unique(unlist(unique(stringr::str_split(units_by_topic,", ")))))
      n_geosets <- length(unique(unlist(sapply(sapply(sapply(avail_subset,function(x){x[,geosets]}),function(x){stringr::str_split(as.character(x),"-|, ")}),unlist))))
    }
    avail_countries <- unique(unlist(stringr::str_split(country_by_topic,", ")))
    avail_geosets <- paste0("Geoset:",unique(unlist(stringr::str_split(unlist(sapply(avail_subset,function(x){x[,geosets]})),", "))))
    avail_countries <- sapply(SUNGEO::available_data[avail_geosets],function(x){avail_countries[avail_countries%in%x[,country_iso3]]})
    suppressWarnings({
      avail_geosets <- data.table::data.table(geosets=avail_geosets)[,country_iso3:=avail_countries[match(geosets,names(avail_countries))]][,desc:=paste0(geosets," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(geosets," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    })
    n_topics <- sum(!grepl("^Geoset",names(avail_subset)))
    avail_topics <- names(avail_subset)[!grepl("^Geoset",names(avail_subset))]
    avail_topics <- data.table::data.table(topics=avail_topics)[,country_iso3:=country_by_topic[match(topics,names(country_by_topic))]][,desc:=paste0(topics," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(topics," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    # Output
    summary_message <- c("-----------------------------------------","Summary of data available through SUNGEO","-----------------------------------------","",paste0("Found data for ",n_countries," countries, ",n_topics," topics"),paste0(n_units," spatial units, ",n_geosets," geographic boundary sets"),"","-------","Topics","-------",avail_topics[,(desc)],"","-------","Geosets","-------",avail_geosets[,(desc)],"","-------","See 'topics' and 'geosets' slots for details on spatial units and years.","-------" )
    print(summary_message)
    out <- list(summary=summary_message,topics=avail_subset[!grepl("^Geoset",names(avail_subset))],geosets=lapply(SUNGEO::available_data[avail_geosets[,geosets]],function(x){x[country_iso3%in%unique(unlist(avail_countries)),]}))
  }

  # By country
  if(!is.null(country_iso3s)&is.null(topics)){
    # Take subset
    avail_subset <- sapply(SUNGEO::available_data,function(x){x[country_iso3%in%c(country_iso3s)]})
    avail_subset <- avail_subset[sapply(avail_subset,nrow)>0]
    # Summarize
    country_by_topic <- sapply(sapply(avail_subset,function(x){x[,country_iso3]}),function(x){paste0(x,collapse=", ")})
    units_by_topic <- sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){paste0(x,collapse=", ")})
    years_by_topic <- sapply(sapply(sapply(avail_subset[!grepl("^Geoset",names(avail_subset))],function(x){x[,year_range]}),unique),function(x){paste0(x,collapse=", ")})
    years_by_geoset <- sapply(sapply(avail_subset[grepl("^Geoset",names(avail_subset))],function(x){x[,geoset_years]}),unique)
    years <- paste0(range(unlist(stringr::str_split(years_by_topic,", |-"))),collapse="-")
    n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
    n_units <- length(unique(unlist(sapply(stringr::str_split(units_by_topic,", "),length))))
    n_geosets <- length(grep("^Geoset",names(avail_subset)))
    avail_geosets <- names(avail_subset)[grep("^Geoset",names(avail_subset))]
    avail_geosets <- data.table::data.table(geosets=avail_geosets)[,country_iso3:=country_by_topic[match(geosets,names(country_by_topic))]][,desc:=paste0(geosets," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(geosets," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    n_topics <- sum(!grepl("^Geoset",names(avail_subset)))
    avail_topics <- names(avail_subset)[!grepl("^Geoset",names(avail_subset))]
    avail_topics <- data.table::data.table(topics=avail_topics)[,country_iso3:=country_by_topic[match(topics,names(country_by_topic))]][,desc:=paste0(topics," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(topics," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    # Output
    summary_message <- c("-----------------------------------------","Summary of data available through SUNGEO","-----------------------------------------","",paste0("Found data for ",n_countries," countries, ",n_topics," topics"),paste0(n_units," spatial units, ",n_geosets," geographic boundary sets"),"","-------","Topics","-------",avail_topics[,(desc)],"","-------","Geosets","-------",avail_geosets[,(desc)],"","-------","See 'topics' and 'geosets' slots for details on spatial units and years.","-------" )
    print(summary_message)
    out <- list(summary=summary_message,topics=avail_subset[!grepl("^Geoset",names(avail_subset))],geosets=avail_subset[grepl("^Geoset",names(avail_subset))])
  }

  # By topic and country
  if(!is.null(country_iso3s)&!is.null(topics)){
    # Take subset
    avail_subset <- sapply(SUNGEO::available_data,function(x){x[country_iso3%in%c(country_iso3s)]})
    avail_subset <- avail_subset[topics]
    avail_subset <- avail_subset[sapply(avail_subset,nrow)>0]
    # Summarize
    if(length(topics)==1){
      country_by_topic <- paste0(sapply(avail_subset,function(x){x[,country_iso3]}),collapse=", ")
      names(country_by_topic) <- topics
      units_by_topic <- paste0(sort(unique(unlist(sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){stringr::str_split(x,", ")})))),collapse=", ")
      names(units_by_topic) <- topics
      years_by_topic <- paste0(range(unlist(sapply(sapply(avail_subset,function(x){x[,year_range]}),function(x){stringr::str_split(x,", |-")}))),collapse="-")
      names(years_by_topic) <- topics
      years <- paste0(range(unlist(stringr::str_split(years_by_topic,", |-"))),collapse="-")
      n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
      n_units <- length(unique(unlist(unique(stringr::str_split(units_by_topic,", ")))))
      n_geosets <- length(unique(unlist(sapply(sapply(sapply(avail_subset,function(x){x[,geosets]}),function(x){stringr::str_split(as.character(x),"-|, ")}),unlist))))
    }
    if(length(topics)>1){
      country_by_topic <- sapply(sapply(avail_subset,function(x){x[,country_iso3]}),function(x){paste0(x,collapse=", ")})
      names(country_by_topic) <- topics
      units_by_topic <- sapply(sapply(sapply(sapply(sapply(sapply(avail_subset,function(x){x[,space_units]}),function(x){stringr::str_split(x,", ")}),unlist),unique),sort),function(x){paste0(x,collapse=", ")})
      years_by_topic <- apply(sapply(sapply(sapply(avail_subset,function(x){x[,year_range]}),function(x){stringr::str_split(as.character(x),"-|, ")}),range),2,function(x){paste0(unique(x),collapse="-")})
      years <- paste0(range(unlist(stringr::str_split(years_by_topic,", |-"))),collapse="-")
      n_countries <- max(sapply(stringr::str_split(country_by_topic,", "),length))
      n_units <- length(unique(unlist(unique(stringr::str_split(units_by_topic,", ")))))
      n_geosets <- length(unique(unlist(sapply(sapply(sapply(avail_subset,function(x){x[,geosets]}),function(x){stringr::str_split(as.character(x),"-|, ")}),unlist))))
    }
    avail_countries <- unique(unlist(stringr::str_split(country_by_topic,", ")))
    avail_geosets <- paste0("Geoset:",unique(unlist(stringr::str_split(unlist(sapply(avail_subset,function(x){x[,geosets]})),", "))))
    avail_countries <- sapply(sapply(SUNGEO::available_data[avail_geosets],function(x){avail_countries[avail_countries%in%x[,country_iso3]]}),function(x){paste0(x,collapse=", ")})
    suppressWarnings({
      avail_geosets <- data.table::data.table(geosets=avail_geosets)[,country_iso3:=avail_countries[match(geosets,names(avail_countries))]][,desc:=paste0(geosets," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(geosets," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    })
    n_topics <- sum(!grepl("^Geoset",names(avail_subset)))
    avail_topics <- names(avail_subset)[!grepl("^Geoset",names(avail_subset))]
    avail_topics <- data.table::data.table(topics=avail_topics)[,country_iso3:=country_by_topic[match(topics,names(country_by_topic))]][,desc:=paste0(topics," (",country_iso3,")")][sapply(stringr::str_split(country_iso3,","),length)>3,desc:=paste0(topics," (",sapply(stringr::str_split(country_iso3,","),length)," countries)")]
    # Output
    summary_message <- c("-----------------------------------------","Summary of data available through SUNGEO","-----------------------------------------","",paste0("Found data for ",n_countries," countries, ",n_topics," topics"),paste0(n_units," spatial units, ",n_geosets," geographic boundary sets"),"","-------","Topics","-------",avail_topics[,(desc)],"","-------","Geosets","-------",avail_geosets[,(desc)],"","-------","See 'topics' and 'geosets' slots for details on spatial units and years.","-------" )
    print(summary_message)
    out <- list(summary=summary_message,topics=avail_subset[!grepl("^Geoset",names(avail_subset))],geosets=lapply(SUNGEO::available_data[avail_geosets[,geosets]],function(x){x[country_iso3%in%unique(unlist(avail_countries)),]}))
  }



  return(out)
}


