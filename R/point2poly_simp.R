#' Point-to-polygon interpolation, simple overlay method
#'
#' Function for assigning values from a source point layer to a destination polygon layer, using simple point-in-polygon overlays
#'
#' @param pointz Source points layer. \code{sf} object.
#' @param polyz Destination polygon layer. Must have identical CRS to \code{pointz}. \code{sf} object.
#' @param varz Names of variable(s) to be assigned from source polygon layer to destination polygons. Character string or vector of character strings.
#' @param char_varz Names of character string variable(s) in \code{varz}. Character string or vector of character strings.
#' @param funz Aggregation function to be applied to variables specified in \code{varz}. Must take as an input a vector \code{x}. Function or list of functions.
#' @param na_val Value to be assigned to missing values. Defaul is \code{NA}. Logical or list.
#' @param drop_na_cols Drop columns with completely missing values. Defaul is \code{FALSE}. Logical.
#' @return Returns a \code{sf} polygon object, with variables from \code{pointz} assigned to the geometries of \code{polyz}.
#' @details Assignment procedures are the same for numeric and character string variables. All variables supplied in \code{varz} are passed directly to the function specified in \code{funz}. If different sets of variables are to be aggregated with different functions, both \code{varz} and \code{funz} should be specified as lists (see examples below).
#' @import sf
#' @importFrom data.table data.table rbindlist as.data.table
#' @importFrom methods as
#' @importFrom dplyr full_join left_join
#' @importFrom purrr reduce
#' @examples
#' # Assignment of a single variable (sums)
#' \dontrun{
#' data(hex_05_deu)
#' data(clea_deu2009_pt)
#' out_1 <- point2poly_simp(pointz=clea_deu2009_pt,
#'                          polyz=hex_05_deu,
#'                          varz="vv1")
#' plot(out_1["vv1"])
#' }
#'
#' # Replace NA's with 0's
#' \dontrun{
#' out_2 <- point2poly_simp(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = "vv1",
#'                          na_val = 0)
#' plot(out_2["vv1"])
#' }
#'
#' # Multiple variables, with different assignment functions
#' \dontrun{
#' out_3 <- point2poly_simp(pointz = clea_deu2009_pt,
#'                          polyz = hex_05_deu,
#'                          varz = list(
#'                            c("to1","pvs1_margin"),
#'                            c("vv1"),
#'                            c("incumb_pty_n","win1_pty_n")),
#'                          funz = list(
#'                            function(x){mean(x,na.rm=TRUE)},
#'                            function(x){sum(x,na.rm=TRUE)},
#'                            function(x){paste0(unique(na.omit(x)),collapse=" | ") }),
#'                          na_val = list(NA_real_,0,NA_character_))
#' }
#' @export
#'

point2poly_simp <- function(pointz,
                            polyz,
                            varz,
                            char_varz=NULL,
                            funz=list(function(x){sum(x,na.rm=TRUE)}),
                            na_val=NA,
                            drop_na_cols=FALSE){
  ################################################
  #Part 1 -  Put variables and functions into list
  ################################################
  #Part i -
  if(length(varz) != length(funz) && is.list(varz) && is.list(funz)){
    stop("ERROR: Length of variable list does not equal length of function list. Pleasse ensure that each set of variables correspsonds to a specific function argument.")
  }

  #Part ii -
  if(class(varz)=="character"){varz <- list(varz)}

  #Part iii -
  if(class(funz)!="list"){funz <- list(funz)}
  funz <- lapply(funz, function(sub_iter){
    if(class(sub_iter)%in%"function" == FALSE && class(sub_iter)%in%'character'){
      funzInt <- get(sub_iter, mode = 'function')
      return(funzInt)
    } else {
      return(sub_iter)
    }
  })

  #Part iv -
  if(class(na_val)=="logical"){na_val <- list(na_val)}

  ################################################
  #Part 2:
  ################################################
  if(sf::st_crs(pointz) != sf::st_crs(polyz)){
    #Part i -
    polyz <- sf::st_transform(polyz, crs = sf::st_crs(pointz))
  }

  ############################
  #Part 3 - Empty points layer
  ############################
  #Part i -
  pointz_dt <- data.table::as.data.table(pointz) #Convert to data table object

  #Part ii -
  pointz_dt_NACols <- sapply(pointz_dt, function(x){  #Generate a Vector for Identifying all NA columns
    all(is.na(x))
  })

  #Part iii - drop NA columns (optional)
  if(drop_na_cols==TRUE){
    pointz_dt <- pointz_dt[,pointz_dt_NACols == FALSE, with = FALSE] #Subset columns
  }
  ###################################
  #Part 4 -  Match points to polygons
  ###################################
  #Part i -
  polyz$polygon_id <- 1:nrow(polyz)

  #Part ii -
  suppressWarnings({
    suppressMessages(
      Point_Overlay <- sf::st_within(pointz,polyz)
    )
  })

  #Part iii -
  ErrorCheck <- sapply(Point_Overlay, function(x) length(x))

  ##########################################
  #Part 5 - Add polygon index to point layer
  ##########################################
  #Part i -
  pointz_dt$polygon_id <- NA
  pointz_dt$polygon_id[ErrorCheck == 1] <- polyz$polygon_id[unlist(Point_Overlay[ErrorCheck == 1])]

  #Part ii -
  pointz$polygon_id <- NA
  pointz$polygon_id[ErrorCheck == 1] <- polyz$polygon_id[unlist(Point_Overlay[ErrorCheck == 1])]

  #######################################
  #Part 6 - Aggregate (numeric variables)
  #######################################
  #Part i -
  ClassTypes_Variables <- data.table::data.table('Varz' = varz, 'Funz' = funz, 'NAval' = na_val)

  #Part ii -
  ClassTypes_Variables <- lapply(1:length(ClassTypes_Variables$Varz), function(sub_iter) {
    Varz_Int <- unlist(ClassTypes_Variables$Varz[sub_iter])

    Funz_Int <- unlist(ClassTypes_Variables$Funz[sub_iter], recursive = TRUE)

    NAValue <- unlist(ClassTypes_Variables$NAval[sub_iter], recursive = TRUE)

    OutputMatrix <- data.table::data.table('Varz' = Varz_Int, 'Funz' = Funz_Int, 'NAval' = NAValue)

    return(OutputMatrix)

  })
  ClassTypes_Variables <- data.table::rbindlist(ClassTypes_Variables)

  #Part iii -
  # sub_iter <- 1
  Aggregation_Matrix <- lapply(1:nrow(ClassTypes_Variables), function(sub_iter){
    #Part a -
    Select_Columns <- pointz_dt[,names(pointz_dt)%in%c(ClassTypes_Variables$Varz[sub_iter], 'polygon_id'), with = FALSE]

    #Part b -
    names(Select_Columns) <- c('var', 'polygon_id')

    #Part c -
    IntermediateFunction <- ClassTypes_Variables$Funz[sub_iter][[1]]

    #Part d -
    suppressWarnings({
      Output <- by(data=Select_Columns,INDICES=Select_Columns$polygon_id,FUN=function(x){IntermediateFunction(x$var)},simplify=TRUE)
      Output <- data.table::data.table(polygon_id=as.numeric(names(Output)),V1=c(Output))
    })

    #Part e -
    names(Output)[2] <- ClassTypes_Variables$Varz[sub_iter]

    #RETURN
    return(Output)
  })

  #Part iv -
  Empty_Sets <- sapply(Aggregation_Matrix, function(x) is.null(x) || grepl('Error', x))
  if(TRUE%in%Empty_Sets){
      Aggregation_Matrix <- Aggregation_Matrix[Empty_Sets == FALSE]
  }

  #Part v -
  Combined_Matrix <- purrr::reduce(Aggregation_Matrix,dplyr::full_join, by = 'polygon_id')

  #Part vi - Remove NA Rows (Failed Matches)
  if(any(is.na(Combined_Matrix$polygon_id))){
      Combined_Matrix <- Combined_Matrix[is.na(Combined_Matrix$polygon_id) == FALSE,]
  }

  #######################################
  #Part 7 -
  #######################################
  #Part i -
  Coordinates <- sf::st_geometry(polyz)
  sf::st_geometry(polyz) <- NULL

  #Part ii -
  polyz_ <- dplyr::left_join(polyz, Combined_Matrix, by = 'polygon_id')

  #Part iii -
  for(iter in 1:length(ClassTypes_Variables$Varz)){
    VectorSubset <- polyz_[,names(polyz_)%in%ClassTypes_Variables$Varz[iter]]

    VectorSubset[is.na(VectorSubset)] <- ClassTypes_Variables$NAval[iter]

    polyz_[,names(polyz_)%in%ClassTypes_Variables$Varz[iter]] <- VectorSubset
  }

  #Part iv -
  sf::st_geometry(polyz_) <- Coordinates

  # Ensure classes are same as in source file (+ coerce char_varz to character)
  classez <- sapply(as.data.frame(pointz,stringsAsFactors=FALSE)[,unique(unlist(varz))],class)
  if(length(char_varz)>0){classez[names(classez)%in%char_varz] <- "character"}
  polyz_ng <- sf::st_drop_geometry(polyz_)
  for(cl0 in 1:length(classez)){polyz_ng[,names(classez[cl0])] <- methods::as(polyz_ng[,names(classez[cl0])],classez[cl0])}
  polyz_ <- sf::st_set_geometry(polyz_ng,polyz_$geometry); rm(polyz_ng)

  #Output
  return(polyz_)
}
