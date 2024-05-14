#' Merge list of tables on common variable(s)
#'
#' Function that finds a set of common columns in a list of tables, and merges the tables on these columns.
#'
#' @param lst List of tables to be merged. List object.
#' @return \code{data.table} object
#' @importFrom data.table data.table as.data.table
#' @examples
#' # Merge list of three tables with different common variables
#' \dontrun{
#' A <- data.table::data.table(month=month.name,year=rep(1991:1992,each=12),A=rnorm(24))
#' B <- data.table::data.table(year=c(1991,1992),B=rbeta(2,1,1))
#' C <- data.table::data.table(month=month.name,C=runif(12))
#' 
#' out_1 <- merge_list(list(A,B,C))
#' out_1
#' }
#' 
#' @export

merge_list <- function(lst){
  while(length(lst) > 1) {
    # Create index of list objects
    idxlst <- seq(from=1, to=length(lst), by=2)
    # Iterate over pairs of list objects
    lst <- lapply(idxlst, function(i){
      if(i==length(lst)){ return(lst[[i]]) }
      # Identify common variables
      cmnvar <- intersect(names(lst[[i]]),names(lst[[i+1]]))
      # Merge tables
      out <- data.table::as.data.table(merge(lst[[i]], lst[[i+1]],by=cmnvar,all=TRUE,allow.cartesian=TRUE))
      # Drop duplicate columns
      out <- out[,.SD,.SDcols=unique(names(out))]
      return(out)
    })
  }
  return(lst[[1]])
}


