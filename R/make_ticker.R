#' Make date ticker
#'
#' Function to create a table of consecutive dates, in SUNGEO-compliant format.
#'
#' @param date_min Start date, in YYYYMMDD format. Default is \code{19000101}. Integer.
#' @param date_max End date, in YYYYMMDD format. Default is today. Integer.
#' @return \code{data.table} object, with seven columns:
#' \itemize{
##'  \item \code{DATE}. Date in YYYYMMDD format. Integer.
##'  \item \code{DATE_ALT}. Date in \code{Date} (YYYY-MM-DD) format. Date.
##'  \item \code{TID}. Date ID, in consecutive integer format. Integer.
##'  \item \code{YRWK}. Week in YYYYWW format. Integer.
##'  \item \code{WID}. Weed ID, in consecutive integer format. Integer.
##'  \item \code{YRMO}. Month in YYYYMM format. Integer.
##'  \item \code{MID}. Month ID, in consecutive integer format. Integer.
##'  \item \code{YEAR}. Year in YYYY format. Integer.
##' }
#' @importFrom data.table data.table as.data.table
#' @examples
#' # All dates from January 1, 1900 to today
#' \dontrun{
#' out_1 <- make_ticker()
#' out_1
#' }
#' 
#' # All dates from January 1, 1200 to today
#' \dontrun{
#' out_2 <- make_ticker(date_min=12000101)
#' out_2
#' }
#'  
#' # All dates from January 1, 1500 to December 31, 1899
#' \dontrun{
#' out_3 <- make_ticker(date_min=15000101, date_max=18991231)
#' out_3
#' }
#' @export


make_ticker <- function(
  date_min=19000101,
  date_max=as.integer(gsub("-","",as.Date(Sys.Date())))
  ){
  # Convert to "Date" format
  date_min_ <- paste0(substr(date_min,1,4),"-",substr(date_min,5,6),"-",substr(date_min,7,8))
  date_max_ <- paste0(substr(date_max,1,4),"-",substr(date_max,5,6),"-",substr(date_max,7,8))
  # Create time table
  if(date_max>=19000101){
    ticker <- data.table::data.table(DATE_ALT =  seq(as.Date("1900-01-01"), as.Date(date_max_), by="days"))[,DATE:=as.integer(gsub("-","",DATE_ALT))][,TID:=1:.N][,WID:=rep(1:.N,each=7)[1:.N]][,YRMO:=as.integer(substr(gsub("-","",DATE),1,6))][,MID:=as.integer(as.factor(YRMO))][,YEAR:=as.integer(substr(gsub("-","",DATE),1,4))]
  }
  # Combine with pre-1900 dates (if applicable)
  if(date_min<19000101&date_max>=19000101){
    ticker2 <- data.table(DATE_ALT =  seq(as.Date(date_min_), as.Date("1899-12-31"), by="days"))[,DATE:=as.integer(gsub("-","",DATE_ALT))][,TID:=1-(.N:1)][,WID:=as.integer(1-rev(rep(1:.N,each=7)[1:.N]))][,YRMO:=as.integer(substr(gsub("-","",DATE),1,6))][,MID:=as.integer(1-as.integer(factor(YRMO,levels=rev(unique(YRMO)))))][,YEAR:=as.integer(substr(gsub("-","",DATE),1,4))]
    ticker <- data.table::as.data.table(rbind(ticker2,ticker))[DATE>=date_min&DATE<=date_max]
    rm(ticker2)
  }
  if(date_min<19000101&date_max<19000101){
    ticker <- data.table(DATE_ALT =  seq(as.Date(date_min_), as.Date("1899-12-31"), by="days"))[,DATE:=as.integer(gsub("-","",DATE_ALT))][,TID:=1-(.N:1)][,WID:=as.integer(1-rev(rep(1:.N,each=7)[1:.N]))][,YRMO:=as.integer(substr(gsub("-","",DATE),1,6))][,MID:=as.integer(1-as.integer(factor(YRMO,levels=rev(unique(YRMO)))))][,YEAR:=as.integer(substr(gsub("-","",DATE),1,4))][DATE>=date_min&DATE<=date_max]
  }
  # Add week codes, re-order columns
  ticker <- data.table::setkey(merge(ticker,data.table::setnames(ticker[,unique(YEAR),by=WID],"V1","YEAR")[!duplicated(WID,fromLast=TRUE)][, YRWK := match(WID, unique(WID)) + (YEAR*1e3), by=YEAR][,YEAR := NULL] ,by="WID")[DATE>=date_min&DATE<=date_max][,.(DATE,DATE_ALT,TID,YRWK,WID,YRMO,MID,YEAR)],TID)
  return(ticker)
}
