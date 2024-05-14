#' Smart numerical rounding function
#'
#' Function to round numerical values with minimal information loss (e.g. to avoid "0.000" values in tables).
#'
#' @param x Vector of values to be rounded. Numeric.
#' @param rnd Requested number of decimal places. Default is 0. Non-negative integer.
#' @param return_char Return rounded values as character string? Default is TRUE. Logical.
#' @return If \code{return_char=TRUE}, returns a character string of same length as \code{x}. If \code{return_char=FALSE}, returns a numerical vector of same length as \code{x}.
#' @details Rounds the values in its first argument to the specified number of decimal places (default 0). If brute-force rounding produces zero values (e.g. "0.00"), the number of decimal places is expanded to include the first significant digit.
#' @examples
#' # Round a vector of numbers, character string output (best for tables) 
#' \dontrun{
#' out_1 <- smart_round(c(.0013,2.3,-1,pi),rnd=2)
#' out_1
#' }
#' 
#' # Round a vector of numbers, numerical output 
#' \dontrun{
#' out_2 <- smart_round(c(.0013,2.3,-1,pi),rnd=2,return_char=FALSE)
#' out_2
#' }
#' @export

smart_round <- function(x,rnd=0,return_char=TRUE){
  x_ <- x
  x_[abs(round(x,rnd))<(1*10^(0-rnd))] <- signif(x[abs(round(x,rnd))<(1*10^(0-rnd))],1)
  x_[abs(round(x,rnd))>=(1*10^(0-rnd))] <- round(x[abs(round(x,rnd))>=(1*10^(0-rnd))],rnd)
  x_[abs(round(x_,rnd+1))>=(1*10^(0-(rnd+1)))] <- round(x_[abs(round(x_,rnd+1))>=(1*10^(0-(rnd+1)))],rnd+1)

  # Change class (if applicable)
  if(return_char==TRUE){
    x_ <- as.character(x_)
  }
  return(x_)
}

