#' Conditional cumulative density function
#' 
#' @name build_u_tree
#' 
#' @description Computing the conditionnal cdf of an observation, knowing the 
#' surrounding ones.
#' 
#' @import utils
#' @import rvinecopulib
#' 
#' 
#' @param u (numeric) a vector of marginal uniforms
#' @param longitudinalDVine An S4 object of class LongitudinalDVine
#' 
#' @return (numeric) The triangle of conditionnal uniforms
#' 
#' @examples 
#' u <- pnorm(stats::arima.sim(list(order = c(1,0,0), ar = 0.85), n = 100))
#' trained_DVine <- fit_LongitudinalDVine(u)
#' build_u_tree(u[1:5], trained_DVine)
#' @export
build_u_tree <- function(u, longitudinalDVine)
{
   assertthat::assert_that(inherits(longitudinalDVine, "LongitudinalDVine"))
   
   .T <- length(u)
   
   max_lag <- min(.T, length(longitudinalDVine))
   
   uu <- array(
      NaN,
      dim = c(max_lag, .T, 2), 
      dimnames = list(NULL, NULL, c("forward", "backward"))
   )
   uu[1,, 1] <- uu[1,, 2] <- head(u, .T)
   
   for (l in 1:(max_lag - 1))
   {
      .pairs <- .pairing_sequence(uu, l)
      
      .c <- longitudinalDVine[[l]]
      
      # forward
      uu[l + 1,,"forward"] <- c(
         hbicop(.pairs, 2, .c$family, .c$rotation, .c$parameters),
         rep(NaN, l)
      )
      
      # Backward
      uu[l + 1,,"backward"] <- c(
         hbicop(.pairs, 1, .c$family, .c$rotation, .c$parameters),
         rep(NaN, l)
      )
      
   }
   
   return(uu)
}
