#'
#' @name LongitudinalDVine_class
#' @rdname LongitudinalDVine_class
#' @title Class to create an instance of a longitudinal D-Vine copula structure.
#' @description Class to create an instance of a longitudinal D-Vine copula structure.
#'
#' @import methods
#' 
#' @export
methods::setClass("LongitudinalDVine", contains = "list")


#' @rdname LongitudinalDVine_class
#' 
#' @param x (list) an embedded list describing the D-vine structure.
#' Each element of the list corresponds to the main characteristics of the
#' `bicop_dist()` function from the `rvinecopulib` package, which are :
#' * `family` (characters), 
#' * `rotation` (An integer between 0 and 365),
#' * `parameters` (A named vector).
#' 
#' See Examples for more details.
#' 
#' @return An S4 object that contains the D-vine copula structure.
#' 
#' @examples 
#' LongitudinalDVine_dist(x = list(
#'    list(
#'       family = "gaussian",
#'       rotation = 0,
#'       parameters = 0.85
#'    ),
#'    list(
#'       family = "indep",
#'       rotation = 0,
#'       parameters = numeric(0)
#'    )))
#' 
#' @export
LongitudinalDVine_dist <- function(x) {
   assertthat::assert_that(inherits(x, "list"))
   
   x <- lapply(x, .family_to_lower)
   
   return(methods::new("LongitudinalDVine", x))
}


.family_to_lower <- function(y) {
   y$family <- tolower(y$family)
   return(y)
}
