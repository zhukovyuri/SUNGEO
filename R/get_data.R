#' Download data from SUNGEO server
#'
#' Function to download data files through the SUNGEO API. Function produces a data.table object, corresponding to the user's choice of countries, topics, sources, and spatial and temporal units.
#'
#' @param country_names Country name(s). Character string (single country) or vector of character strings (multiple countries).
#' @param country_iso3 Country code (ISO 3166-1 alpha-3). Character string (single country) or vector of character strings (multiple countries).
#' @param geoset Name of geographic boundary set. Can be one of \code{"GADM"} (Database of Global Administrative Areas), \code{"GAUL"} (Global Administrative Unit Layers), \code{"geoBoundaries"}, \code{"GRED"} (GeoReferenced Electoral Districts Datasets), \code{"HEXGRID"} (SUNGEO Hexagonal Grid), \code{"MPIDR"} (Max Planck Institute for Demographic Research Population History GIS Collection), \code{"NHGIS"} (National Historical Geographic Information System), \code{"PRIOGRID"} (PRIO-GRID 2.0), \code{"SHGIS"} (SUNGEO Historical GIS). Default is \code{"geoBoundaries"}. Character string.
#' @param geoset_yr Year of geographic boundaries. See \code{get_info()['geosets']} for availability. Default is \code{2020}. Integer.
#' @param space_unit Geographic level of analysis. Can be one of \code{"adm0"} (country), \code{"adm1"} (province), \code{"adm2"} (district), \code{"cst"} (GRED electoral constituency), \code{"hex05"} (SUNGEO Hexagonal Grid cell), \code{"prio"} (PRIO-GRID cell). See \code{get_info()['geosets']} for availability by geoset, country and topic. Default is \code{"adm1"}. Character string.
#' @param time_unit Temporal level of analysis. Can be one of \code{"year"}, \code{"month"}, \code{"week"}. See \code{get_info()['topics']} for availability by topic. Default is \code{"year"}. Character string.
#' @param topics Data topics. See \code{get_info()['summary']} for full list. Character string (single topic) or vector of character strings (multiple topics).
#' @param year_min Time range of requested data: start year. See \code{get_info()['topics']} for availability by topic. Default is \code{1990}. Integer.
#' @param year_max Time range of requested data: end year. See \code{get_info()['topics']} for availability by topic. Default is \code{2017}. Integer.
#' @param print_url Print url string of requested data to console? Default is \code{TRUE}. Logical.
#' @param print_time Print processing time for API query to console? Default is \code{TRUE}. Logical.
#' @param error_stop Error handling. If \code{TRUE}, function terminates request if an error is encountered. If \code{FALSE}, error is skipped and error message is recorded in a new \code{message} column. Default is \code{FALSE}. Logical.
#' @param by_topic Break query down by topic and country? If \code{TRUE}, a separate request is sent to the API for each country and topic, and the results are combined on the client side. This ensures that data that are available for some, but not all countries are returned, rather than resulting in a failed request. If \code{FALSE}, a single request is sent to the API for all countries and topics, and the results are combined on the server side. Only data that are available for all countries are returned. Default is \code{TRUE}. Logical.
#' @param skip_missing Skip missing data topics? If \code{TRUE}, missing data topics are skipped, columns are populated with NAs, and corresponding error message is recorded in a new \code{message} column. If \code{FALSE}, returns NULL results for missing topics. Default is \code{TRUE}. Logical.
#' @param cache_param Store cached query on server? This can speed up processing for repeated queries. Default is \code{FALSE}. Logical.
#' @param short_message Shorten error messages? If \code{TRUE}, a short, informative error message is recorded in the \code{message} column. If \code{FALSE}, full error message is recorded. Default is \code{TRUE}. Logical.
#' @return data.table object, with requested data from SUNGEO API.
#' @importFrom data.table data.table as.data.table setnames rbindlist copy
#' @importFrom jsonlite fromJSON
#' @importFrom RCurl getURL
#' @importFrom stringr str_split str_extract_all
#' @importFrom dplyr last
#' @seealso \code{\link{get_info}}
#' @examples
#' # Single country, single topic
#' \dontrun{
#' out_1 <- get_data(country_name="Afghanistan",topics="Demographics:Population:GHS")
#' out_1
#' }
#'
#  # Multiple countries, multiple topics
#' \dontrun{
#' out_2 <- get_data(
#' 	country_name=c("Afghanistan","Moldova"),
#' 	topics=c("Demographics:Ethnicity:EPR","Demographics:Population:GHS"))
#' out_2
#' }
#'
#' # Other boundary sets, spatial and time units
#' \dontrun{
#' out_3 <- get_data(
#' 	country_name="Albania",
#' 	topics="Weather:AirTemperatureAndPrecipitation:NOAA",
#' 	geoset="GAUL",geoset_yr=1990,space_unit="adm2",time_unit="month",
#' 	year_min=1990,year_max=1991)
#' out_3
#' }
#' @export
get_data <- function(
    country_names=NULL,
    country_iso3=NULL,
    geoset="geoBoundaries",
    geoset_yr=2020,
    space_unit="adm1",
    time_unit="year",
    topics=NULL,
    year_min=1990,
    year_max=2017,
    print_url=TRUE,
    print_time=TRUE,
    error_stop=FALSE,
    by_topic=TRUE,
    skip_missing =TRUE,
    cache_param=FALSE,
    short_message=TRUE){

  ####
  # Error messages
  ####

  if(is.null(country_names)&&is.null(country_iso3)){
    stop("ERROR: Please provide at least one country_names or country_iso3")
  }

  if(is.null(topics)){
    stop("ERROR: Please specify at least one topic")
  }


  ####
  # Define internal functions
  ####

  # Bindings
  .SD <- .N <- DATE <- MID <-TID <- WID <- YEAR <- YRMO <- YRWK <- error <- topic <- country_name_alt <- NULL

  # Time matrix
  make_ticker <- function(
    date_min=15000101,
    date_max=as.numeric(gsub("-","",as.Date(Sys.Date())))
  ){
    date_max_ <- paste0(substr(date_max,1,4),"-",substr(date_max,5,6),"-",substr(date_max,7,8))
    ticker <- data.table::data.table(DATE =  gsub("-","",seq(as.Date("1900-01-01"), as.Date(date_max_), by="days")))[,TID:=1:.N][,WID:=rep(1:.N,each=7)[1:.N]][,YRMO:=substr(gsub("-","",DATE),1,6)][,MID:=as.numeric(as.factor(YRMO))][,YEAR:=substr(gsub("-","",DATE),1,4)]
    ticker2 <- data.table::data.table(DATE =  gsub("-","",seq(as.Date("1500-01-01"), as.Date("1899-12-31"), by="days")))[,TID:=1-(.N:1)][,WID:=1-rev(rep(1:.N,each=7)[1:.N])][,YRMO:=substr(gsub("-","",DATE),1,6)][,MID:=1-as.numeric(factor(YRMO,levels=rev(unique(YRMO))))][,YEAR:=substr(gsub("-","",DATE),1,4)]
    ticker <- data.table::as.data.table(rbind(ticker2,ticker))[,lapply(.SD,as.numeric)]
    ticker <- ticker[DATE>=date_min&DATE<=date_max]
    rm(ticker2)
    ticker <- merge(ticker,data.table::setnames(ticker[,unique(YEAR),by=WID],"V1","YEAR")[!duplicated(WID,fromLast=TRUE)][, YRWK := match(WID, unique(WID)) + (YEAR*1e3), by=YEAR][,YEAR := NULL],by="WID")
    return(ticker)
  }

  # Merge list items on common variables
  merge_list <- function(lst){
    while(length(lst) > 1) {
      idxlst <- seq(from=1, to=length(lst), by=2)
      lst <- lapply(idxlst, function(i) {
        if(i==length(lst)) { return(lst[[i]]) }
        cmnvar <- intersect(names(lst[[i]]),names(lst[[i+1]]))
        out <- data.table::as.data.table(merge(lst[[i]], lst[[i+1]],by=cmnvar,all=TRUE,allow.cartesian=TRUE))
        out <- out[,.SD,.SDcols=unique(names(out))]
        return(out)
      })
    }
    lst[[1]]
  }


  ####
  # Download and integrate data
  ####

  # Country codes
  if (length(country_iso3) == 0 & length(country_names) > 0) {
    country_iso3 <- unique(SUNGEO::cc_dict[match(country_names,SUNGEO::cc_dict[,country_name_alt]),country_iso3])
    country_iso3 <- country_iso3[!is.na(country_iso3)]
  }

  if(by_topic==FALSE){
    # Create query
    url_string <- paste0("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/",paste0(country_iso3,collapse=","),"/",geoset,"/",geoset_yr,"/",space_unit,"/",time_unit,"?include=",paste0(topics,collapse=","),"&date=",year_min,"-",year_max,"&cacheEnabled=",cache_param)
    if(print_url==TRUE){
      print(paste0("Fetching: ",gsub("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/","",url_string,fixed=TRUE)))
    }
    # Download, parse JSON, convert to dt
    t1 <- Sys.time()
    json_return <- jsonlite::fromJSON(gsub("\\}\\{","},{",RCurl::getURL(url_string)))
    if("error"%in%names(json_return)&error_stop==TRUE){
      stop(paste0("ERROR ",json_return$error$code,". ",json_return$error$message," (",gsub("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/","",url_string,fixed=TRUE),")"))
    } else if("error"%in%names(json_return)&error_stop==FALSE){
      print(paste0("ERROR ",json_return$error$code," (",gsub("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/","",url_string,fixed=TRUE),")"))
      sungeo.data <- data.table::data.table(error=json_return$error$code,message=json_return$error$message)
    } else if(!"error"%in%names(json_return)){
      sungeo.data <- data.table::as.data.table(data.table::rbindlist(json_return[["data"]],fill=TRUE))
    }

    t2 <- Sys.time(); if(print_time==TRUE){print(t2-t1)}
  } else {
    # Time precision table
    time_invariant <- c("Demographics:Ethnicity:GREG","Infrastructure:Roads:gRoads","Terrain:Elevation:ETOPO1","Terrain:LandCover:GLCC")
    topic_time <- data.table::data.table(topic=c("Demographics:Ethnicity:EPR","Demographics:Ethnicity:GREG","Demographics:Population:GHS","Elections:LowerHouse:CLEA","Events:PoliticalViolence:ABADarfur","Events:PoliticalViolence:ACLED","Events:PoliticalViolence:BeissingerProtest","Events:PoliticalViolence:BeissingerRiot","Events:PoliticalViolence:BeissingerUkraine","Events:PoliticalViolence:COCACW","Events:PoliticalViolence:ESOCAfghanistanWITS","Events:PoliticalViolence:ESOCIraqSIGACT","Events:PoliticalViolence:ESOCIraqWITS","Events:PoliticalViolence:ESOCMexicoDrugRelatedMurders","Events:PoliticalViolence:ESOCMexicoHomicide","Events:PoliticalViolence:ESOCPakistanBFRS","Events:PoliticalViolence:ESOCPakistanWITS","Events:PoliticalViolence:GED","Events:PoliticalViolence:Lankina","Events:PoliticalViolence:NIRI","Events:PoliticalViolence:NVMS","Events:PoliticalViolence:PITF","Events:PoliticalViolence:SCAD","Events:PoliticalViolence:yzCaucasus2000","Events:PoliticalViolence:yzChechnya","Events:PoliticalViolence:yzLibya","Events:PoliticalViolence:yzUkraine2014","Infrastructure:Roads:gRoads","Infrastructure:NightLights:DMSP","PublicHealth:Covid19:JHUCSSEC19","Terrain:Elevation:ETOPO1","Terrain:LandCover:GLCC","Weather:AirTemperatureAndPrecipitation:NOAA"),
                                         time_unit=c("YEAR","YEAR","YEAR","YEAR|MONTH","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR","YEAR","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR|MONTH|WEEK","YEAR","YEAR","YEAR|MONTH|WEEK","YEAR","YEAR","YEAR|MONTH")
    )
    t1 <- Sys.time()
    # tp0 <- topics[1]
    # is0 <- country_iso3[1]
    data_list <- lapply(country_iso3,function(is0){
      topic_list <- lapply(topics,function(tp0){
        # Adjust time precision
        if(!grepl(topic_time[topic%in%tp0,time_unit],time_unit)){
          time_unit_ <- topic_time[topic%in%tp0,sapply(stringr::str_split(time_unit,"\\|"),dplyr::last)]
        } else {time_unit_ <- time_unit}
        if(tp0%in%time_invariant){
          url_string <- paste0("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/",is0,"/",geoset,"/",geoset_yr,"/",space_unit,"/",time_unit_,"?include=",tp0,"&cacheEnabled=",cache_param)
        } else {
          url_string <- paste0("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/",is0,"/",geoset,"/",geoset_yr,"/",space_unit,"/",time_unit_,"?include=",tp0,"&date=",year_min,"-",year_max,"&cacheEnabled=",cache_param)
        }
        if(print_url==TRUE){
          print(paste0("Fetching: ",gsub("https://api-sungeo-org-sungeo-api.apps.gnosis.lsa.umich.edu/data/","",url_string,fixed=TRUE)))
        }
        # Download, parse JSON, convert to dt
        tryCatch({
          # csv_return_ <- data.table::fread(RCurl::getURL(paste0(url_string,"&csv=true")))
          # names(csv_return_)[1]<-"POLYGON_ID"
          json_return_ <- RCurl::getURL(url_string)
          if(grepl("Fatal error",json_return_)){
            json_return <- list()
            json_return$error <- data.table::data.table(code=500,message=gsub("<.*?>|\\n", "", json_return_))
          }else{
            json_return <-  jsonlite::fromJSON(gsub("\\}\\{","},{",json_return_))
          }
        },error=function(e){
          print(paste0("ERROR: ",url_string," ",e))
          json_return <<- list(); json_return$error <<- data.table::data.table(code=500,message=paste0(e))
        })

        # }
        if("error"%in%names(json_return)&error_stop==TRUE){
          if(short_message == TRUE & grepl("not all of the provided topics",json_return$error$message)){
            errm <- unlist(stringr::str_extract_all(json_return$error$message,"\\[.+?\\]"))
            errm <- paste0(grep("\\:",errm,value=TRUE)," missing for ",errm[!grepl("\\:",errm)])
          } else {
            errm <- json_return$error$message
          }
          stop(paste0("ERROR ",json_return$error$code,". ",errm))
        } else if("error"%in%names(json_return)&error_stop==FALSE){
          if(short_message == TRUE & grepl("not all of the provided topics",json_return$error$message)){
            errm <- unlist(stringr::str_extract_all(json_return$error$message,"\\[.+?\\]"))
            errm <- paste0(grep("\\:",errm,value=TRUE)," missing for ",errm[!grepl("\\:",errm)])
          } else {
            errm <- json_return$error$message
          }
          print(paste0("ERROR ",json_return$error$code,", ",errm))
          sungeo.data <- data.table::data.table(error=json_return$error$code,message=errm)
        } else if(!"error"%in%names(json_return)){
          # Parse JSON
          sungeo.data <- data.table::as.data.table(data.table::rbindlist(json_return,fill=TRUE))
          # Add date fields if applicable
          if(any(names(sungeo.data)%in%c("WID","DATE","TID","YRMO","MID","YEAR","YRWK"))){
            ticker <- make_ticker(date_min=paste0(year_min,"0101"),date_max=min(paste0(year_max,"1231"),as.numeric(gsub("-","",as.Date(Sys.Date())))))[,c("YEAR","YRMO","YRWK"):=lapply(.SD,as.character),.SDcols=c("YEAR","YRMO","YRWK")]
            common_var <- intersect(names(ticker),names(sungeo.data))
            if(length(common_var)==1){
              if(grepl("YEAR",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YEAR)&get(common_var)%in%sungeo.data[,get(common_var)]]
              }
              if(grepl("MONTH",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YRMO)&get(common_var)%in%sungeo.data[,get(common_var)]]
              }
              if(grepl("WEEK",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YRWK)&get(common_var)%in%sungeo.data[,get(common_var)]]
              }
              sungeo.data <- merge(data.table::copy(ticker)[,eval(common_var):=as.character(get(common_var))],data.table::copy(sungeo.data)[,eval(common_var):=as.character(get(common_var))],by=common_var,all=TRUE,allow.cartesian=TRUE)
            }
            if(length(common_var)>1){
              if(grepl("YEAR",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YEAR)&get(common_var)%in%sungeo.data[,get(common_var)]]
              }
              if(grepl("MONTH",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YRMO)&YRMO%in%sungeo.data[,YRMO]]
              }
              if(grepl("WEEK",time_unit,ignore.case=TRUE)){
                ticker <- ticker[!duplicated(YRWK)&WID%in%sungeo.data[,WID]]
              }
              sungeo.data <- merge(data.table::copy(ticker)[,eval(common_var):=lapply(.SD,as.character),.SDcols=common_var],data.table::copy(sungeo.data)[,eval(common_var):=lapply(.SD,as.character),.SDcols=common_var],by=common_var,all=TRUE,allow.cartesian=TRUE)
            }
          }
        }

        sungeo.data
      })
      # Skip missing topics?
      good_results <- which(sapply(topic_list,function(x){nrow(x)>0&!"error"%in%names(x)}))
      errorz <- which(sapply(topic_list,function(x){"error"%in%names(x)}))
      zeroz <- which(sapply(topic_list,function(x){nrow(x)==0}))
      if(skip_missing == TRUE){
        if(length(good_results)>0){
          if(length(errorz)>0&length(zeroz)>0){
            error_msgz <- unlist(sapply(topic_list[errorz],function(x){x[,substring(sapply(stringr::str_split(paste0(error,". ",message)," in \\'on clause"),"[",1),1,200)]}))
            topic_out <- merge_list(topic_list[good_results])[,message:=paste0("ERRORS: ",paste0(paste0(topics[errorz],". Error code ",error_msgz),collapse=", "),".   ZERO ROWS: ",paste0(topics[zeroz],collapse=", "),".")]
          }
          if(length(errorz)>0&length(zeroz)==0){
            error_msgz <- unlist(sapply(topic_list[errorz],function(x){x[,substring(sapply(stringr::str_split(paste0(error,". ",message)," in \\'on clause"),"[",1),1,200)]}))
            topic_out <- merge_list(topic_list[good_results])[,message:=paste0("ERRORS: ",paste0(paste0(topics[errorz],". Error code ",error_msgz),collapse=", "),".")]
          }
          if(length(errorz)==0&length(zeroz)>0){
            topic_out <- merge_list(topic_list[good_results])[,message:=paste0("ZERO ROWS: ",paste0(topics[zeroz],collapse=", "),".")]
          }
          if(length(errorz)==0&length(zeroz)==0){
            topic_out <- merge_list(topic_list[good_results])
          }
        } else {
          topic_out <- NULL
        }
      }
      topic_out
    })
    if(print_url==TRUE){
      print(paste0("Combining..."))
    }
    sungeo.data <- data.table::rbindlist(data_list,fill=TRUE)
    # Remove redundancies
    suppressWarnings({
      if(grepl("YEAR",time_unit,ignore.case=TRUE)){
        sungeo.data[,c("YRMO","MID","YRWK","WID","DATE","TID") := NULL]
      }
      if(grepl("MONTH",time_unit,ignore.case=TRUE)){
        sungeo.data[,c("YRWK","WID","DATE","TID") := NULL]
      }
      if(grepl("WEEK",time_unit,ignore.case=TRUE)){
        sungeo.data[,c("DATE","TID") := NULL]
      }
    })
    # Remove duplicates
    sungeo.data <- unique(sungeo.data)
    t2 <- Sys.time(); if(print_time==TRUE){print(t2-t1)}
  }

  # Return dt
  return(sungeo.data)
}
