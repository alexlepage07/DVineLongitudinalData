# Modeling Longitudinal Data Using a Pair-Copula Decomposition of Serial Dependence

## An implementation of [[Smith and al. 2010](https://doi.org/10.1198/jasa.2010.tm09572)]

By [Alexandre Lepage](alexlepage07@hotmail.com), Laval University

Last update: 2022-12-14


**Disclaimer** This package has not been revised. Some bugs might come up and
it is not meant to be published on CRAN. This is the result of a university 
project.

----

### Context

In [[Smith and al. 2010](https://doi.org/10.1198/jasa.2010.tm09572)],
they explain how to use pair-copulas to model the dependence structure between
longitudinal data. However, in their paper they make the assumption that all
copulas are of the same family, but with different parameters. In our work, we
relax this assumption by assuming that, for a single lag, all copulas are
identically distributed, but across different lags, we can have a variety of
copula families. In addition, we use the assumption that, when at a given lag
the dependence is negligible, the independence is assumed for all further lags.

As a result, the Algorithm 1 from the paper works really well. However, the
implementation of the method in this package is not robust as it hardly retrieve
the original copula structure when we simulate time series. For more details, see
the simulation study here: `inst/Simulation_study.R`.

--- 

### Package installation

All you need to install the package is this file : `DVineSD_0.0.1.zip`
by using this command once you have pulled the repo:
 
```
install.packages("~/DVineSD_0.0.1.zip", repos = NULL, type = "win.binary")
```

### Package content

- `DVineSD.Rproj`: You should always work from an R project. Before opening an R file, open this one.
- `R`: This is the main folder where all the R code is saved
- `man`: All the package documentation is written here
- `inst`: You will find the dataset that [Smith and al. 2010] used for their reel case example.


### Content of the `R` folder

- `conditionnal_cdf.R`: The implementation of Algorithm 1 from [Smith and al. 2010]
- `fit_LongitudinalDVine.R`: The code to fit a D-Vine copula for longitudinal data
- `simul_serial_dependance.R`: The code to simulate a sequence of time-dependant observations.
- `LongitudinalDVine_class.R`: Where the formal S4 LongitudinalDVine class is declared, so we could transport information from the `fit_LongitudinalDVine.R` method to the `simul_serial_dependance.R`. Otherwise, it is possible to simply declare an 
instance of the LongitudinalDVine class, and then simulate from it.
- `internals.R`: Internal functions used in the other R files.

### Content of the `inst` folder

- `extdata`: The dataset that [Smith and al. 2010] used for their reel case example.
- `Simulation_study.R`: This is the simulation study that I made to test the package


---

**Before sharing, remember that this package is a prototype and is not meant for production !**
