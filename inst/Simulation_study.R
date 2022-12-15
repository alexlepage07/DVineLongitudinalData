# Simulation_study.R

# A script to try to simulate a time series and then to retrive the parameters


# Packages ---------------------------------------------------------------------


library(data.table)
library(DVineSD)
library(rvinecopulib)
library(xtable)


# Inputs -----------------------------------------------------------------------


seq_length <- 1e+3
max_lag <- 4L
seed <- 20221128


# Functions --------------------------------------------------------------------


fit_longitudinal_DVine_rvinecopulib <- function(
      x, 
      max_lag = 4L, 
      par_method = "mle",
      selcrit = "mbic",
      .margins_controls=list(xmin=0, xmax=1)
   ) {
   
   par_method <- match.arg(par_method, c("itau", "mle"))
   selcrit <- match.arg(selcrit, c("aic", "bic", "mbic"))
   
   
   data <- .format_autocorrelated_data(x, max_lag)
   
   dvine_longitudinal <- rvinecopulib::vine(
      data, 
      margins_controls = .margins_controls,
      copula_controls = list(
         family_set = c("onepar", "indep"), 
         par_method = par_method,
         selcrit = selcrit,
         structure = rvinecopulib::dvine_structure(ncol(data))
   ))
   
   dt <- summary(dvine_longitudinal)$copula[,c(1:2, 6:8, 10)]
   dt$parameters <- as.numeric(dt$parameters)
   
   print(xtable::xtable(dt, align = c("lcc|lccr")), row.names=FALSE)
   
   return(invisible(dvine_longitudinal))
}


print_results_DVineSD <- function(trained_Dvine) {
   
   .trained_Dvine <- trained_Dvine
   out <- matrix(nrow = length(trained_Dvine), ncol =3)
   
   for(l in seq_along(trained_Dvine)) {
      if (.trained_Dvine[[l]]$family == "indep"){
         .trained_Dvine[[l]]$parameters <- NA
      }
      .trained_Dvine[[l]]$parameters <- as.numeric(.trained_Dvine[[l]]$parameters)
      out[l,] <- rbind(unlist(.trained_Dvine[[l]]))
   }
   
   colnames(out) <- names(trained_Dvine[[1]])
   out <- as.data.frame(out)

   out$parameters <- as.numeric(out$parameters)
   
   print(xtable::xtable(out, digits = 2L, align = "rlcr"), row.names=FALSE)
   
   return(invisible(NULL))
}
   

.format_autocorrelated_data <- function(x, max_lag) {
   
   data <- x <- as.numeric(x)
   for (l in 1:(max_lag)) {
      data <- cbind(head(data, -1), tail(x, -l))
   }
   
   return(data)
}


# 1 Testing with an independence scenario --------------------------------------


set.seed(seed)
x <- runif(seq_length)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)


# 2 Test of the training algorithm on an AR(2) time series dependence ----------


set.seed(seed)
x <- stats::arima.sim(list(ar = c(0.8, -0.6)), n = seq_length)
x <- pnorm(x)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)


# 3 Test of the training algorithm on an AR(3) time series dependence ----------


set.seed(seed)
x <- stats::arima.sim(list(ar = c(0.8, -0.6, 0.4)), n = seq_length)
x <- pnorm(x)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)


# 4 Test of the training algorithm on an AR(4) time series dependence ----------


set.seed(seed)
x <- stats::arima.sim(list(ar = c(0.8, -0.6, 0.4, -0.3)), n = seq_length)
x <- pnorm(x)

# Our function:
trained_Dvine <- DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag)
print_results_DVineSD(trained_Dvine)

# rvinecopulib:
fitted_rvinecopulib <- fit_longitudinal_DVine_rvinecopulib(x, max_lag)


# Verifying that we can retrieve the dependence structure from the initial model.
y <- DVineSD::simul_serial_dependence(trained_Dvine, nseq=seq_length, seed=seed)
z <- rvine(seq_length, fitted_rvinecopulib, cores = 4L)

xtable(cor(.format_autocorrelated_data(x, max_lag)), digits = 4) # Original
xtable(cor(.format_autocorrelated_data(y, max_lag)), digits = 4) # Simulated
xtable(cor(z), digits = 4) # Simulated with rvinecopulib


# 5 Testing the simul_serial_dependence with 4 Gaussian trees ------------------


Dvine_copula <- DVineSD::LongitudinalDVine_dist(x = list(
   list(family = "gaussian", rotation = 0, parameters =  0.8),
   list(family = "gaussian", rotation = 0, parameters = -0.6),
   list(family = "gaussian", rotation = 0, parameters =  0.4),
   list(family = "gaussian", rotation = 0, parameters = -0.3),
))
x <- DVineSD::simul_serial_dependence(Dvine_copula, nseq=seq_length, seed=seed)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)


# 6 Testing the simul_serial_dependence with 2 archimedian trees ---------------


Dvine_copula <- DVineSD::LongitudinalDVine_dist(x = list(
   list(family = "clayton", rotation = 0, parameters =  5),
   list(family = "frank",   rotation = 0, parameters = -4)
))
x <- DVineSD::simul_serial_dependence(Dvine_copula, nseq=seq_length, seed=seed)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)


plot(bicop_dist("clayton", 0, 9))
plot(bicop_dist("frank",   0, 28.09))


# 7 Testing the simul_serial_dependence with 3 archimedian trees ---------------


Dvine_copula <- DVineSD::LongitudinalDVine_dist(x = list(
   list(family = "clayton", rotation = 0, parameters =  5),
   list(family = "frank",   rotation = 0, parameters = -4),
   list(family = "gumbel",  rotation = 0, parameters =  3)
))
x <- DVineSD::simul_serial_dependence(Dvine_copula, nseq=seq_length, seed=seed)

# Our function:
print_results_DVineSD(DVineSD::fit_LongitudinalDVine(x, max_lag = max_lag))

# rvinecopulib:
fit_longitudinal_DVine_rvinecopulib(x, max_lag)
