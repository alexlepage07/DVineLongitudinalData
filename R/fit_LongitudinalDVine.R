#' Fitting a D-Vine copula on a longitudinal data set.
#' 
#' @import methods
#' @import utils
#' @import rvinecopulib
#' 
#' @param x (numeric) a vector of longitudinal observations 
#' @param pmargin (function) the cdf of the marginal distribution
#' @param ... arguments of pmargin
#' @param max_lag (integer) the maximum lag the dependence structure can support
#' 
#' @return an S4 object containing a list of all fitted bivariate copulas. 
#' As an attribute, you can retrieve a version of the `pmargin` value with
#' its parameters nested inside the function.
#' 
#' @examples 
#' # Simulating an AR(1) time series
#' x <- stats::arima.sim(list(order = c(1,0,0), ar = 0.85), n = 500L)
#' x <- pnorm(x)
#' fit_LongitudinalDVine(x)
#' 
#' 
#' # Fit a longitudinal Dvine on an archimedian copula
#' Dvine_copula <- DVineSD::LongitudinalDVine_dist(x = list(
#'    list(family = "frank",  rotation = 0, parameters = 3)
#' ))
#' x <- DVineSD::simul_serial_dependence(Dvine_copula, nseq = 500L)
#' fit_LongitudinalDVine(x, max_lag = 2L)
#' 
#' @export

fit_LongitudinalDVine <- function(x, pmargin, ..., max_lag)
{
   if (missing(pmargin)) {
      pmargin <- stats::ecdf(x)
   } else {
      pmargin <- function(x) pmargin(x, ...)
   }
   
   u <- pmargin(x)
   
   .T <- length(u)
   
   max_lag <- if (missing(max_lag)) .T else max_lag
   max_lag <- min(.T, max_lag)

   uu <- array(
      NaN,
      dim = c(max_lag + 1, .T, 2), 
      dimnames = list(NULL, NULL, c("forward", "backward"))
   )

   trained_bicop <- list()
   
   uu[1,, 1] <- uu[1,, 2] <- u
   
   pb <- txtProgressBar(max = max_lag, style = 3)
   setTxtProgressBar(pb, 0)
   
   for (l in 1:max_lag)
   {
      .pairs <- .pairing_sequence(uu, l)
      
      .c <- rvinecopulib::bicop(
         .pairs, 
         family_set = c("onepar", "indep"), 
         par_method = "itau", 
         selcrit = "aic",
         cores = 4L
      )
      
      # Test of independence of Mantel-Haenszel:
      # Null hypothesis: Independence
      if(.MantelHaenszel_test(.pairs)$pvalue > 0.05) {
         .c <- rvinecopulib::bicop_dist(family = "indep")
      }
      
      trained_bicop[[l]] <- list(
            family = .c$family,
            rotation = .c$rotation,
            parameters = .c$parameters
      )
      
      # forward
      uu[l + 1,,"forward"] <- c(
         rvinecopulib::hbicop(.pairs, 2, .c$family, .c$rotation, .c$parameters),
         rep(NaN, l)
      )

      # Backward
      uu[l + 1,,"backward"] <- c(
         rvinecopulib::hbicop(.pairs, 1, .c$family, .c$rotation, .c$parameters),
         rep(NaN, l)
      )
      
      
      if(.c$family == "indep") {
         setTxtProgressBar(pb, max_lag)
         break
      }
      
      setTxtProgressBar(pb, l)
   }
   
   close(pb)
   
   
   trained_Dvine <- LongitudinalDVine_dist(trained_bicop)
   attr(trained_Dvine, "pmargin") <- pmargin
   
   return(trained_Dvine)
}


.pairing_sequence <- function(uu, .lag=1) 
{
   .l <- 1:sum(!is.na(uu[.lag, , 1]))
   
   .pairs <- cbind(
      forward  = uu[.lag, utils::head(.l, -1), "forward"],
      backward = uu[.lag, utils::tail(.l, -1), "backward"]
   )
   
   return(as.matrix(.pairs))
}


.MantelHaenszel_test <- function(.pairs) {
   out <- list()
   out$stat <- (nrow(.pairs) - 1) * cor(.pairs)[2]^2
   out$pvalue <- 1-pchisq(out$stat, 1)
   return(out)
}
