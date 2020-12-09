#' @importFrom httr handle config status_code headers GET warn_for_status
#' @importFrom rlang abort

## html_sessionSunGeo Function
## Function that replaces rvest::html_session()

html_sessionSunGeo <- function(url, ...) {
  session <- structure(
    list(
      handle   = httr::handle(url),
      config   = c(..., httr::config(autoreferer = 1L)),
      url      = NULL,
      back     = character(),
      forward  = character(),
      response = NULL,
      html     = new.env(parent = emptyenv(), hash = FALSE)
    ),
    class = "session"
  )
  request_GET(session, url)
}

print.session <- function(x, ...) {
  cat("<session> ", x$url, "\n", sep = "")
  cat("  Status: ", httr::status_code(x$response), "\n", sep = "")
  cat("  Type:   ", httr::headers(x)$`Content-Type`, "\n", sep = "")
  cat("  Size:   ", length(x$response$content), "\n", sep = "")
}

request_GET <- function(x, url, ...) {
  x$response <- httr::GET(url, x$config, ..., handle = x$handle)
  x$html <- new.env(parent = emptyenv(), hash = FALSE)
  x$url <- x$response$url

  httr::warn_for_status(x$response)

  x
}


# ## udconvertSunGeo Function
# ## Function that replaces udunits2::udconvert()
# udconvertSunGeo <-
#   function(x, u1, u2) {
#     if (! ud.are.convertible(u1, u2)) {
#       stop(paste("Units", u1, "and", u2, "are not convertible"))
#     }
#     ## Filter out NA's before passing them to the C function
#     ## since it can't handle them
#     rv <- rep(NA, length(x))
#     i <- which(! is.na(x))
#
#     len <- length(i)
#     c.rv <- .C('R_ut_convert',
#                as.double(x)[i],
#                as.integer(len),
#                as.character(u1),
#                as.character(u2),
#                converted=double(len)
#     )
#     rv[i] <- c.rv$converted
#     ## If it's a matrix/vector or anything else, convert it back to it's original type
#     attributes(rv) <- attributes(x)
#     return(rv)
#   }
#
# .onLoad <- function(libname, pkgname) {
#   ## By default, configure udunits with path set (presumably) by the
#   ## user through the UDUNITS2_XML_PATH environment variable
#   .C('R_ut_init', as.integer(0))
#   if (!ud.have.unit.system()) {
#     ## Failing that, override it with the in-package XML file
#     p0 <- system.file("share/udunits2.xml", package="udunits2")
#     Sys.setenv(UDUNITS2_XML_PATH=p0)
#     .C('R_ut_init', as.integer(1))
#     ## If *that* fails, give the user some instructions for how to remedy
#     ## the problem
#     if (!ud.have.unit.system()) {
#       warning(
#         "Failed to read udunits system database: udunits2 will not work properly.\n Please set the UDUNITS2_XML_PATH environment variable before attempting to load the package")
#     }
#   }
# }
#
# .onAttach <- function(libname, pkgname) {
#   msg <- "udunits system database read"
#   p0 <- Sys.getenv("UDUNITS2_XML_PATH")
#   if (p0 != "") {
#     msg <- paste(msg, "from", p0)
#   }
#   packageStartupMessage(msg)
# }
#
# ud.are.convertible <-
#   function(u1, u2) {
#     if (! (ud.is.parseable(u1) && ud.is.parseable(u2))) {
#       return(FALSE)
#     }
#     rv <- .C('R_ut_are_convertible',
#              as.character(u1),
#              as.character(u2),
#              convertible=logical(1))
#     return(rv$convertible)
#   }
#
# ud.get.name <-
#   function(unit.string) {
#     stopifnot(ud.is.parseable(unit.string))
#     rv <- .C('R_ut_get_name',
#              as.character(unit.string),
#              ud.name=character(length=1))
#     return(rv$ud.name)
#   }
#
# ud.get.symbol <-
#   function(unit.string) {
#     stopifnot(ud.is.parseable(unit.string))
#     rv <- .C('R_ut_get_symbol',
#              as.character(unit.string),
#              ud.symbol=character(length=1))
#     return(rv$ud.symbol)
#   }
#
# ud.is.parseable <-
#   function(unit.string) {
#     rv <- .C('R_ut_is_parseable',
#              as.character(unit.string),
#              parseable=logical(1))
#     return(rv$parseable)
#   }
#
# ud.set.encoding <-
#   function(enc.string) {
#     .C('R_ut_set_encoding',
#        as.character(enc.string))
#     return()
#   }
#
# ud.have.unit.system <-
#   function() {
#     rv <- .C('R_ut_has_system',
#              exists=logical(1))
#     return(rv$exists)
#   }



## all_ofSunGeo Function
## Function that replaces tidyselect::all_of()
## We utilize the function from tidyselec::all_of(), which calls vctrs::vec_as_location(), and ellipsis::check_dots_empty()
all_ofSunGeo <- function (x, ...)
{
  if (is.function(x)) {
    if (!missing(...))
      dots <- dots(parent.frame())
    if (length(dots) == 0) {
      return()
    }
    action_dots(action = abort, message = "`...` is not empty.",
                dot_names = names(dots), note = "These dots only exist to allow future extensions and should be empty.",
                .subclass = "rlib_error_dots_nonempty")
    .Call(vctrs_as_location, i = x, n = 0L, names = NULL, loc_negative = "invert",
          loc_oob = "error", loc_zero = "remove", missing =  c("propagate", "error"),
          arg = NULL)
    rlang::abort("Internal error: `all_of()` should have failed sooner")
  }
  x
}



# #
# ## fromJSONSunGeo Function
# ## Function that replaces jsonlite::fromJSON()

# fromJSONSunGeo <- function(txt, simplifyVector = TRUE, simplifyDataFrame = simplifyVector,
#                            simplifyMatrix = simplifyVector, flatten = FALSE, ...) {

#   # check type
#   if (!is.character(txt) && !inherits(txt, "connection")) {
#     stop("Argument 'txt' must be a JSON string, URL or file.")
#   }

#   # overload for URL or path
#   if (is.character(txt) && length(txt) == 1 && nchar(txt, type="bytes") < 2084 && !validate(txt)) {
#     if (grepl("^https?://", txt, useBytes=TRUE)) {
#       loadpkg("curl")
#       h <- curl::new_handle(useragent = paste("jsonlite /", R.version.string))
#       curl::handle_setheaders(h, Accept = "application/json, text/*, */*")
#       txt <- curl::curl(txt, handle = h)
#     } else if (file.exists(txt)) {
#       # With files we can never know for sure the encoding. Lets try UTF8 first.
#       # txt <- raw_to_json(readBin(txt, raw(), file.info(txt)$size));
#       txt <- file(txt)
#     }
#   }

#   # call the actual function (with deprecated arguments)
#   parse_and_simplify(txt = txt, simplifyVector = simplifyVector, simplifyDataFrame = simplifyDataFrame,
#                      simplifyMatrix = simplifyMatrix, flatten = flatten, ...)
# }

# parse_and_simplify <- function(txt, simplifyVector = TRUE, simplifyDataFrame = simplifyVector,
#                                simplifyMatrix = simplifyVector, flatten = FALSE, unicode = TRUE, validate = TRUE, bigint_as_char = FALSE, ...){

#   if(!missing(unicode)){
#     message("Argument unicode has been deprecated. YAJL always parses unicode.")
#   }

#   if(!missing(validate)){
#     message("Argument validate has been deprecated. YAJL automatically validates json while parsing.")
#   }

#   # parse
#   obj <- parseJSON(txt, bigint_as_char)

#   # post processing
#   if (any(isTRUE(simplifyVector), isTRUE(simplifyDataFrame), isTRUE(simplifyMatrix))) {
#     return(simplify(obj, simplifyVector = simplifyVector, simplifyDataFrame = simplifyDataFrame,
#                     simplifyMatrix = simplifyMatrix, flatten = flatten, ...))
#   } else {
#     return(obj)
#   }
# }


# # Backward compatiblity
# fromJSON_string <- parse_and_simplify





# parseJSON <- function(txt, bigint_as_char = FALSE) {
#   if(inherits(txt, "connection")){
#     parse_con(txt, bigint_as_char)
#   } else {
#     parse_string(txt, bigint_as_char)
#   }
# }

# #' @useDynLib SUNGEO R_parse
# parse_string <- function(txt, bigint_as_char){
#   if (length(txt) > 1) {
#     txt <- paste(txt, collapse = "\n")
#   }
#   .Call("R_parse", txt, bigint_as_char)
# }


# #simplify (from JSON)
# simplify <- function(x, simplifyVector = TRUE, simplifyDataFrame = TRUE, simplifyMatrix = TRUE,
#                      simplifyDate = simplifyVector, homoList = TRUE, flatten = FALSE, columnmajor = FALSE,
#                      simplifySubMatrix = simplifyMatrix) {

#   #This includes '[]' and '{}')
#   if (!is.list(x) || !length(x)) {
#     return(x)
#   }

#   # list can be a dataframe recordlist
#   if (isTRUE(simplifyDataFrame) && is.recordlist(x)) {
#     mydf <- simplifyDataFrame(x, flatten = flatten, simplifyMatrix = simplifySubMatrix)
#     if(isTRUE(simplifyDate) && is.data.frame(mydf) && is.datelist(mydf)){
#       return(parse_date(mydf[["$date"]]))
#     }
#     return(mydf)
#   }

#   # or a scalar list (atomic vector)
#   if (isTRUE(simplifyVector) && is.null(names(x)) && is.scalarlist(x)) {
#     return(list_to_vec(x))
#   }

#   # apply recursively
#   out <- lapply(x, simplify, simplifyVector = simplifyVector, simplifyDataFrame = simplifyDataFrame,
#                 simplifyMatrix = simplifySubMatrix, columnmajor = columnmajor, flatten = flatten)

#   # fix for mongo style dates turning into scalars *after* simplifying
#   # only happens when simplifyDataframe=FALSE
#   if(isTRUE(simplifyVector) && is.scalarlist(out) && all(vapply(out, inherits, logical(1), "POSIXt"))){
#     return(structure(list_to_vec(out), class=c("POSIXct", "POSIXt")))
#   }

#   # test for matrix. Note that we have to take another look at x (before
#   # list_to_vec on its elements) to differentiate between matrix and vector.
#   if (isTRUE(simplifyMatrix) && isTRUE(simplifyVector) && is.matrixlist(out) && all(unlist(vapply(x, is.scalarlist, logical(1))))) {
#     if(isTRUE(columnmajor)){
#       return(do.call(cbind, out))
#     } else {
#       #this is currently the default
#       return(do.call(rbind, out))
#     }
#   }

#   # Simplify higher arrays
#   if (isTRUE(simplifyMatrix) && is.arraylist(out)){
#     if(isTRUE(columnmajor)){
#       return(array(
#         data = do.call(cbind, out),
#         dim = c(dim(out[[1]]), length(out))
#       ));
#     } else {
#       #this is currently the default
#       return(array(
#         data = do.call(rbind, lapply(out, as.vector)),
#         dim = c(length(out), dim(out[[1]]))
#       ));
#     }
#   }

#   # try to enfoce homoList on unnamed lists
#   if (isTRUE(homoList) && is.null(names(out))) {
#     # coerse empty lists, caused by the ambiguous fromJSON('[]')
#     isemptylist <- vapply(out, identical, logical(1), list())
#     if (any(isemptylist) & !all(isemptylist)) {
#       # if all the others look like data frames, coerse to data frames!
#       if (all(vapply(out[!isemptylist], is.data.frame, logical(1)))) {
#         for (i in which(isemptylist)) {
#           out[[i]] <- data.frame()
#         }
#         return(out)
#       }

#       # if all others look like atomic vectors, unlist all
#       if (all(vapply(out[!isemptylist], function(z) {
#         isTRUE(is.vector(z) && is.atomic(z))
#       }, logical(1)))) {
#         for (i in which(isemptylist)) {
#           out[[i]] <- vector(mode = typeof(out[[which(!isemptylist)[1]]]))
#         }
#         return(out)
#       }
#     }
#   }

#   # convert date object
#   if( isTRUE(simplifyDate) && is.datelist(out) ){
#     return(parse_date(out[["$date"]]))
#   }

#   # return object
#   return(out)
# }

# is.matrixlist <- function(x) {
#   isTRUE(is.list(x)
#          && length(x)
#          && is.null(names(x))
#          && all(vapply(x, is.atomic, logical(1)))
#          && all.identical(vapply(x, length, integer(1)))
#          #&& all.identical(vapply(x, mode, character(1))) #this fails for: [ [ 1, 2 ], [ "NA", "NA" ] ]
#   );
# }

# is.arraylist <- function(x) {
#   isTRUE(is.list(x)
#          && length(x)
#          && is.null(names(x))
#          && all(vapply(x, is.array, logical(1)))
#          && all.identical(vapply(x, function(y){paste(dim(y), collapse="-")}, character(1)))
#   );
# }

# is.datelist <- function(x){
#   isTRUE(is.list(x)
#          && identical(names(x), "$date")
#          && (is.numeric(x[["$date"]]) || is.character(x[["$date"]]))
#   );
# }

# parse_date <- function(x){
#   if(is.numeric(x)){
#     return(structure(x/1000, class=c("POSIXct", "POSIXt")))
#   } else if(is.character(x)) {
#     #tz is not vectorized, so assume all() are the same
#     is_utc <- ifelse(all(grepl("Z$", x)), "UTC", "")
#     return(as.POSIXct(strptime(x, format = '%Y-%m-%dT%H:%M:%OS', tz = is_utc)))
#   } else {
#     return(x)
#   }
# }

# all.identical <- function(x){
#   if(!length(x)) return(FALSE)
#   for(i in x){
#     if(x[1] != i) return(FALSE)
#   }
#   return(TRUE)
# }


# #is.recordlist (from JSON)
# is_recordlist_c <- function(x){
#   .Call("C_is_recordlist", x)
# }

# is_recordlist_r <- function(x) {
#   if (!(is.unnamedlist(x) && length(x))) {
#     return(FALSE)
#   }
#   at_least_one_object = FALSE
#   for(i in x){
#     if(!(is.namedlist(i) || is.null(i))) return(FALSE)
#     if(!at_least_one_object && is.namedlist(i)) at_least_one_object <- TRUE
#   }
#   return(at_least_one_object)
# }

# is.recordlist <- is_recordlist_c;

# is.namedlist <- function(x) {
#   isTRUE(is.list(x) && !is.null(names(x)))
# }

# is.unnamedlist <- function(x) {
#   isTRUE(is.list(x) && is.null(names(x)))
# }


# #simplifyDataFrame (from JSON)
# simplifyDataFrame <- function(recordlist, columns, flatten, simplifyMatrix) {

#   # no records at all
#   if (!length(recordlist)) {
#     if (!missing(columns)) {
#       return(as.data.frame(matrix(ncol = length(columns), nrow = 0, dimnames = list(NULL,
#                                                                                     columns))))
#     } else {
#       return(data.frame())
#     }
#   }

#   # only empty records and unknown columns
#   if (!any(vapply(recordlist, length, integer(1), USE.NAMES = FALSE)) && missing(columns)) {
#     return(data.frame(matrix(nrow = length(recordlist), ncol = 0)))
#   }

#   # find columns if not specified
#   if (missing(columns)) {
#     columns <- unique(unlist(lapply(recordlist, names), recursive = FALSE, use.names = FALSE))
#   }

#   # Convert row lists to column lists. This is the heavy lifting
#   # columnlist <- lapply(columns, function(x) lapply(recordlist, "[[", x))
#   # Now slighlty optimized
#   columnlist <- transpose_list(recordlist, columns)

#   # simplify vectors and nested data frames
#   columnlist <- lapply(columnlist, simplify, simplifyVector = TRUE, simplifyDataFrame = TRUE,
#                        simplifyMatrix = FALSE, simplifySubMatrix = simplifyMatrix, flatten = flatten)

#   # check that all elements have equal length
#   columnlengths <- unlist(vapply(columnlist, function(z) {
#     ifelse(length(dim(z)) > 1, nrow(z), length(z))
#   }, integer(1)))
#   n <- unique(columnlengths)
#   if (length(n) > 1) {
#     stop("Elements not of equal length: ", paste(columnlengths, collapse = " "))
#   }

#   # add the column names before flattening
#   names(columnlist) <- columns

#   # flatten nested data frames
#   if(isTRUE(flatten)) {
#     dfcolumns <- vapply(columnlist, is.data.frame, logical(1))
#     if(any(dfcolumns)){
#       columnlist <- c(columnlist[!dfcolumns], do.call(c, columnlist[dfcolumns]))
#     }
#   }

#   # make into data frame
#   class(columnlist) <- "data.frame"

#   # set row names
#   if("_row" %in% names(columnlist)) {
#     rn <- columnlist[["_row"]];
#     columnlist["_row"] <- NULL;

#     # row.names() casts double to character which is undesired.
#     if(is.double(rn)) {
#       rn <- as.integer(rn);
#     }

#     # Replace missing values with numbers
#     rn_na <- is.na(rn)
#     if(sum(rn_na) > 0){
#       rn[rn_na] <- paste0("NA_", seq_len(sum(rn_na)))
#     }

#     # data frames MUST have row names
#     if(any(duplicated(rn))){
#       warning('Duplicate names in "_row" field. Data frames must have unique row names.', call. = FALSE)
#       if(is.character(rn)) {
#         row.names(columnlist)  <- make.unique(rn)
#       } else {
#         row.names(columnlist) <- seq_len(n)
#       }
#     } else {
#       row.names(columnlist) <- rn;
#     }
#   } else {
#     row.names(columnlist) <- seq_len(n)
#   }

#   return(columnlist)
# }


# #is.scalarlist (from JSON)
# is_scalarlist_r <- function(x) {
#   if(!is.list(x)) return(FALSE)
#   for(i in x){
#     if(!is.atomic(i) || length(i) > 1) return(FALSE)
#   }
#   return(TRUE)
# }
# is.scalarlist <- is_scalarlist_r;


# #' #' @useDynLib jsonlite C_is_scalarlist
# #' is_scalarlist_c <- function(x){
# #'   .Call("C_is_scalarlist", x)
# #' }
# #' is.scalarlist <- is_scalarlist_c;


# #transpose_list (from JSON)
# transpose_list <- function(x, names) {
#   .Call("C_transpose_list", x, names)
# }


# #list_to_vec (from JSON)
# list_to_vec <- function(x) {
#   isdates <- is_datelist(x)
#   out <- unlist(null_to_na(x), recursive = FALSE, use.names = FALSE)
#   if(isdates && is.numeric(out)){
#     structure(out, class = c("POSIXct", "POSIXt"))
#   } else{
#     out
#   }
# }


# #is_datelist (from JSON)
# is_datelist <- function(x){
#   .Call("C_is_datelist", x)
# }


# #null_to_na (from JSON)
# null_to_na <- function(x) {
#   .Call("C_null_to_na", x)
# }
