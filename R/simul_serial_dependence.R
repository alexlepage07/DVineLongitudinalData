#' Simulate realisations of a time series
#'
#' @name simul_serial_dependence
#' 
#' @description Simulating time-dependant sequence
#' 
#' @import utils
#' @import parallel
#' @import rvinecopulib
#' @import stats
#' 
#' @include internals.R
#' 
#' 
#' @param longitudinalDVine A S3 object of class LongitudinalDVine
#' @param x (numeric) a vector of marginal uniforms on which we condition
#' @param qmargin (function) the marginal quantile function of x
#' @param nsim (integer) the desired number of simulations
#' @param nseq (integer) the length of the simulated sequence
#' @param seed (integer) a simulation seed
#' @param cores (integer) the number of available threads to do parallel computation
#' @param ... further arguments to qmargin
#' 
#' @return (numeric) a matrix of dimension nsim x n_seq of simulated values
#' 
#' @examples 
#' Dvine_copula <- DVineSD::LongitudinalDVine_dist(x = list(
#'    list(
#'       family = "clayton",
#'       rotation = 0,
#'       parameters = 0.9
#'    ),
#'    list(
#'       family = "indep",
#'       rotation = 0,
#'       parameters = numeric(0)
#'    )
#' ))
#' x <- DVineSD::simul_serial_dependence(Dvine_copula, nseq = 10L)

#' @export
simul_serial_dependence <- function(
   longitudinalDVine, x, qmargin, nsim=1L, nseq=100L, seed=32L, cores=4L, ...)
{
   # nsim = 10L; nseq=5L; seed=32L; cores=4L
   
   assertthat::assert_that(inherits(longitudinalDVine, "LongitudinalDVine"))
   
   qmargin <- if(missing(qmargin)) function(y) y else qmargin
   
   if(!missing(x)) {
      pmargin <- attr(longitudinalDVine, "pmargin")
      u <- pmargin(tail(x, length(longitudinalDVine)))
      uu <- build_u_tree(u, longitudinalDVine)
      uu <- .extend_array(uu, dim(uu) + c(nseq, nseq, 0))
      
   } else {
      uu <- array(NaN, dim = c(nseq, nseq, 2))
      dimnames(uu) <- list(NULL, NULL, c("backward", "forward")) 
      x <- numeric()
   }
   
   set.seed(seed)
   
   on.exit(parallel::stopCluster(cl))
   cl <- parallel::makeCluster(cores)
   parallel::clusterEvalQ(cl, require(rvinecopulib))
   
   simulations <- t(parallel::parSapply(
      cl, seq(to=nsim, by=1), .ssd_single_simul, uu, longitudinalDVine, nseq))
   
   
   return(structure(
      qmargin(simulations, ...),
      past_observations = x,
      class = c("sequence_simulation", "matrix")
   ))
} 


setOldClass(c("sequence_simulation", "matrix"))


.ssd_single_simul <- function(i, uu, longitudinalDVine, nseq) 
{
   simul_id <- tail(1:dim(uu)[1], nseq)
   uu[simul_id, 1, "backward"] <- stats::runif(nseq)
   
   for (.t in simul_id) 
   {
      if (.t == 1) 
      {
         uu[1,1,"forward"] <- uu[1,1,"backward"]
         next
      }
      for (j in 1:(.t-1)) 
      {
         lag <- .t-j
         
         if(lag > length(longitudinalDVine)) {
            .c <- list(family = "indep", rotation = 0, parameters = numeric())
         } else {
            .c <- longitudinalDVine[[lag]]
         }
         
         
         # Computing backward utj using the inverse h function
         u1 <- uu[lag, j, "forward"]
         u2 <- uu[lag+1, j, "backward"]
         utj <- rvinecopulib::hbicop(
            c(u1, u2), 1, .c[[1]], .c[[2]], .c[[3]], TRUE)
         uu[lag, j+1, "backward"] <- utj
         
         if(lag == 1) {
            uu[lag, j+1, "forward"] <- uu[lag, j+1, "backward"]
         }
         
         # Compute forward ujt using the h function and utj
         u1 <- uu[lag, j, "forward"]
         u2 <- uu[lag, j+1, "backward"]
         ujt <- rvinecopulib::hbicop(
            c(u1, u2), 2, .c[[1]], .c[[2]], .c[[3]], FALSE)
         uu[lag+1, j, "forward"] <- ujt
      }
   }
   
   return(as.vector(uu[1, simul_id, 1]))
}
